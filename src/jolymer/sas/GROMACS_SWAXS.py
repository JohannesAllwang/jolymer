from scipy import optimize
import numpy as np
import re
from os.path import join
from matplotlib.colors import LogNorm
from pylab import cm
import matplotlib.pyplot as plt
import pandas as pd

import shutil
import os

import sasmodels
from sasmodels import data as sasmodels_data
import pyFAI
import fabio
from pathlib import Path

from .. import os_utility as osu
from .SAXS_Measurement import SAXS_Measurement

from dataclasses import dataclass

def _colorbar(mappable):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib.pyplot as plt
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax)
    plt.sca(last_axes)
    return cbar

@dataclass
class GROMACS_SWAXS(SAXS_Measurement):

    instrument = 'gromacs swaxs'

    rawpath: str= ''
    path: str = '~'
    filename: str='waxs_final.xvg'
    scectra_name: str='waxs_spectra.xvg'
    name: str='unnamed'
    angular_unit: str='$\\mathrm{nm^{-1}}$'
    gromacs_path: str='not defined'
    maxq: float=np.inf
    color: str=None
    marker: str=''
    linestyle: str='-'
    mdpath: str=None
    spectra_filename: str=None
    proj_filename: str=None
    rmsf_filename: str=None
    eigenval_filename: str=None
    swaxspot_filename: str=None

    def get_eigenval(self, mdpath=None, **kwargs):
        if mdpath is None:
            mdpath = self.mdpath
        # if eigenval_filename is None:
        eigenval_filename = self.eigenval_filename
        if eigenval_filename is None:
            return None
        # print(mdpath, eigenval_filename)
        filename = Path(mdpath, eigenval_filename)
        with open(filename) as file:
            out = {}
            i = 1
            data = []
            for line in file:
                if not line.startswith(('@', '#', '&')):
                    data.append([float(x) for x in line.split()])
        df = pd.DataFrame({'pca': [int(row[0]) for row in data],
                           'RMSF': [row[1] for row in data]}).set_index('pca')
        return df

    def get_swaxspot(self, path=None, filename=None):
        if path is None:
            path = self.mdpath
        if filename is None:
            filename = self.swaxspot_filename
        with open(Path(path,filename)) as f:
            outdict = {'time': []}
            qs = []
            for line in f:
                if len(line.split('legend')) > 1:
                    qeqq = line.split('legend')[1]
                    q = float(line.split('=')[1].replace('"', ''))
                    qs.append(q)
                    outdict[q] = []
                elif len(line.split()) > 3 and not line[0] in ['#', '@']:
                    numbers = [float(n) for n in line.split()]
                    outdict['time'].append(numbers.pop(0))
                    for number, q in zip(numbers, qs):
                        outdict[q].append(number)
        df = pd.DataFrame(outdict)
        df.time = df.time / 1000
        return df

    def get_filename(self):
        return join(self.path, self.filename)

    def get_spectra_filename(self):
        return join(self.path, self.spectra_name)

    def get_rg(self, path=None, filename='rg.xvg'):
        if path is None:
            path = self.mdpath
        full_path = Path(path, filename)
        gdict = self.read_xvg(full_path, include1=False, every_n=1, max_out=np.inf,
                 names=['step', 'Rg', 'Rx', 'Ry', 'Rz'],
                 onlyfirst=False)
        df = gdict['dfs'][0]
        df['time'] = df.step/1000
        return df

    def get_data(self, engine='pandas', **kwargs):
        """
        Reads the data file at self.path / self.filename
        Alternatively, a path can be provided as a keyword argument.
        """
        out =  self.read_xvg(self.get_filename(), names=['q', 'I', 'err_I'],
                             include1=True, onlyfirst=True)
        return out['dfs'][0]

    def get_gromacs(self, path=None, filename='waxs_final.xvg'):
        gromacs_path = path
        if path is None:
            gromacs_path = self.gromacs_path
        filename = join(gromacs_path, filename)
        with open(filename) as f:
            i = 0
            for line in f:
                if line == '@type xydy\n':
                    skiprows = i+1
                if line == '&\n':
                    print(line)
                    nrows = i - skiprows
                    break
                i += 1
        print('skiprows', skiprows)
        print('nrows', nrows)
        df = pd.read_csv(filename,
                 sep=r'\s+', names=['q', 'I', 'errI'], skiprows=skiprows,
                 nrows = nrows,
                 dtype=np.float64)
        return df

    def get_proj(self, mdpath=None, filename='proj.xvg', pc='pc', **kwargs):
        eigenval_df = self.get_eigenval(mdpath=mdpath)
        # N_eigenval = len(eigenval_df)
        if mdpath is None:
            mdpath = self.mdpath
        filename = Path(mdpath, filename)
        with open(filename) as file:
            out = {}
            i = 1
            data = []
            for line in file:
                if line.startswith(('&')):
                    df = pd.DataFrame({'x': [row[0] for row in data],
                                       'y': [row[1] for row in data]})
                    df['time'] = df.x/1000
                    out['time'] = df.time
                    # RMSFi = float(eigenval_df.loc[i].iloc[0])
                    out[f'{pc}{i}'] = df.y
                    data = []
                    i += 1
                if not line.startswith(('@', '#', '&')):
                    data.append([float(x) for x in line.split()])
        outdf = pd.DataFrame(out)
        return outdf

    # def get_spectra():
    #     file = self.get_spectra_filename()
    #     out = m.get_waxs_spectra(file, every_n=every_n, max_out=max_out)

    ## Plot ##

    def plot_waxspot3D(self, cmap='viridis', **kwargs):
        df = self.get_swaxspot(**kwargs)
        from mpl_toolkits.mplot3d import Axes3D  # not strictly necessary, but keeps IDEs happy
        T, Q = np.meshgrid(df.time, df.columns[1:].astype(float), indexing="ij")  # shape: (n_timepoints, n_qvalues)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(T, Q, df.iloc[:, 1:].values, cmap=cmap)
        ax.set_xlabel("Time [ns]")
        ax.set_ylabel("$q$ [nm$^{-1}$]")
        ax.set_zlabel("$E_{saxs}$ [kJ/mol]")
        return ax

    def plot_spectra(self, filename=None, path=None,
                     plot=True, get_Rg=False, maxqfit=np.inf,
                     rerun_filename=None,
                      every_n=1, max_out=1000, index_list=None):
        if path is None:
            path = self.mdpath
        if filename is None:
            filename = self.spectra_filename
        outdict = {'time': [],
                   'df': [],
                   'chi2': [],
                   'Rg': [],
                   'err_Rg': [],
                   'I0': [],
                   'err_I0': []}
        file = Path(path, filename)
        if filename is None:
            file = self.get_spectra_filename()
        df_data = self.get_data()
        if plot:
            fig = plt.figure(figsize=(3,3.5))
            gs = fig.add_gridspec(2, 1, height_ratios=[3, 1])  # 3:1 ratio
            ax = fig.add_subplot(gs[0])
            ax_res = fig.add_subplot(gs[1], sharex=ax)
            ax_res.plot([df_data.q.min(), df_data.q.max()], [0,0])
            ax = self.plot_data(ax=ax, label=f'{self.name} data', marker='o', linestyle='')
        out = self.get_waxs_spectra(file, every_n=every_n, max_out=max_out)
        df_data['I'] = df_data.I
        df_data['err_I'] = df_data.err_I
        dfs = out['dfs'][1::]
        times = out['times'][1::]
        if not rerun_filename is None:
            rerun_out = self.get_waxs_spectra(Path(path, rerun_filename), every_n=every_n, max_out=max_out)
            rerun_dfs = rerun_out['dfs'][1::]

        if not index_list is None:
            dfs = [out['dfs'][i] for i in index_list]
            times = [out['times'][i] for i in index_list]
        for df, time in zip(dfs, times):
            df['fit'] = df.I
            df = df[df.q>0.15]
            scale = None
            offset = None
            if maxqfit < df.q.max():
                if rerun_filename is None:
                    df_maxqfit = df.copy()
                else:
                    df_maxqfit =
                df_maxqfit = df_maxqfit[df_maxqfit.q < maxqfit]
                odict = self.scale_and_offset_fit(df_data, df_maxqfit)
                scale = odict['scale']
                offset = odict['offset']
            odict = self.scale_and_offset_fit(df_data, df, scale=scale, offset=offset)
                 # p0=None, scale=None, offset=None, bounds=(-np.inf, np.inf))
            df = odict['df']
            outdict['df'].append(df)
            chi2 = odict['chi2']
            outdict['chi2'].append(chi2)
            time_ns = time * 2e-6
            outdict['time'].append(time_ns)
            if plot:
                ax.errorbar(df.q, df.I, fmt='-',
                            label=f'$t = {time_ns:.0f}$ ns; $\\chi^2={chi2:.1f}$')
                ax_res.errorbar(df.q, df.res, fmt='o-',
                            label='$\\chi^2={:.0f}$'.format(chi2))
            else:
                ax=None
            if get_Rg:

                Rgdict = SAXS_Measurement.get_rg(self, df=df, plot=plot, ax=ax,
                                                 qmin=0.1, qmax=1)
                # print('rgdict', Rgdict)
                for key in Rgdict.keys():
                    # print(Rgdict[key])
                    outdict[key].append(Rgdict[key])

        if plot:
            ax.legend(fontsize='small')
            ax_res.set_xlabel(ax.get_xlabel())
            ax_res.set_ylabel('Residuals')
            ax_res.legend()
            ax_res.set_xscale('log')
        # ax.set_ylim(1e3, 1e7)
        return outdict

    def get_all(self, mdpath=None,
                plot_spectra=False,
                spectra_filename=None,
                proj_filename=None,
                rmsf_filename=None,
                waxspot_name=None):
        if mdpath is None:
            mdpath = self.mdpath
        if spectra_filename is None:
            spectra_filename = self.spectra_filename
        if proj_filename is None:
            proj_filename = self.proj_filename
        if rmsf_filename is None:
            rmsf_filename = self.rmsf_filename
        if not spectra_filename is None:
            spectra_dict = self.plot_spectra(filename=spectra_filename,
                                         path=mdpath,
                                         every_n=1,
                                         max_out=10000,
                                         plot=plot_spectra,
                                         get_Rg=True)
            spectra_df = pd.DataFrame({'time': spectra_dict['time'],
                                   'chi2': spectra_dict['chi2'],
                                   'I0': spectra_dict['I0'],
                                   'err_I0': spectra_dict['err_I0'],
                                   'Rg': spectra_dict['Rg'],
                                   'err_Rg': spectra_dict['err_Rg']})
        if not proj_filename is None:
            proj_df = self.get_proj(mdpath=mdpath, filename=proj_filename)
        if not rmsf_filename is None:
            try:
                rmsf_df = self.get_proj(mdpath=mdpath, filename=rmsf_filename, pc='rmsf')
            except:
                print(f'no file called {rmsf_filename}')
        # print('proj_df', proj_df)
        # print('spectra_df', spectra_df)
        if proj_filename is None:
            return spectra_df
        elif spectra_filename is None:
            return proj_df
        else:
            merged_df = pd.merge_asof(proj_df, spectra_df[['time', 'chi2', 'Rg']], on='time', direction='forward')
            return merged_df

    def plot_rg(self, ax=None, rgkwargs={}, **kwargs):
        out = {}
        if ax is None:
            fig, ax = plt.subplots(figsize=(3,1))
        df = self.get_rg(**rgkwargs)
        out['df'] = df
        annotation = '{}, $R_g = {:.2f} \\pm {:.2f}$ nm'.format(
                        self.name, df.Rg.mean(), df.Rg.std())
        ax.errorbar(df.time, df.Rg,
                    marker='',
                    linestyle='-',
                    **kwargs)
        ax.annotate(annotation, xy=(0.01, 0.8), xycoords='axes fraction',
                    fontsize='small', color='#e5e7fa', bbox=dict(facecolor='#6e79ad'))

        ax.set_ylabel('$R_G$ [nm]')
        ax.set_xlabel('time [ns]')
        out['ax'] = ax
        return out


    @staticmethod
    def read_xvg(filename, include1=False, every_n=1, max_out=np.inf,
                 names=['q', 'I', 'errI', 'I2', 'I3', 'I4'],
                 onlyfirst=False):
        out = {'dfs': [],
               'times': []}
        with open(filename) as f:
            start_row = 0
            time = 1
            plot = 0 if include1 else 1
            for i, line in enumerate(f):
                if line[0] == '@':
                    start_row = i+1
                elif line[0] == '#':
                    start_row = i+1
                elif line == '\n':
                    start_row = i+1
                elif line == '&\n':
                    end_row = i-2
                    if start_row>0 and (plot%every_n)==0 and len(out['dfs'])<max_out:
                        df = pd.read_csv(filename,
                                         skiprows=start_row,
                                         nrows=end_row-start_row, sep=r'\s+',
                                         names=names)
                        if onlyfirst:
                            return df
                        out['dfs'].append(df)
                        out['times'].append(time)
                    time += 1
                    plot += 1
            df = pd.read_csv(filename, skiprows=start_row, sep=r'\s+',
                            names=names)
            out['dfs'].append(df)
            out['times'].append(time)
        return out



    @staticmethod
    def get_crysol(path):
        df = pd.read_csv(path, skiprows=2, sep=r'\s+',
                         names=['q', 'I', 'errI', 'fit'])
        return df

    # @staticmethod
    # def get_chi2(df, df_fit):
    #     if len(df) != len(df_fit):
    #         print('Interpolating since df and df_fit donnot have the same len')
    #         from scipy.interpolate import interp1d
    #         df_fit = df_fit[df_fit.q<df.q.max()]
    #         print(df.q.max(), df_fit.q.max())
    #         interp_function = interp1d(df.q, df.I)
    #         df_fit['data'] = interp_function(df_fit.q)
    #         err_interp_function = interp1d(df.q, df.err_I)
    #         df_fit['err_data'] = err_interp_function(df_fit.q)
    #         df_fit = df_fit[df_fit.err_data>0]
    #     else:
    #         df_fit['data'] = df.I
    #         df_fit['err_data'] = df.err_I
    #     chi2 = np.sum(((df_fit.data - df_fit.fit) / df_fit.err_data)**2 /
    #                          (len(df_fit) - 1))
    #     return df_fit, chi2
