"""
"""

from scipy import optimize
import numpy as np
import re
from os.path import join
from matplotlib.colors import LogNorm
from pylab import cm
import matplotlib.pyplot as plt
import pandas as pd
import subprocess

import shutil
import os

import sasmodels
from sasmodels import data as sasmodels_data
import pyFAI
import fabio

from .. import database_operations as dbo
from ..Measurement import Measurement
from .. import os_utility as osu

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
class SAXS_Measurement(Measurement):

    instrument = 'no instrument'

    rawpath: str= ''
    path: str = '~'
    filename: str='merge_001.dat'
    name: str='unnamed'
    short_name: str='unnamed'
    angular_unit: str='A' # or 'nm'
    gromacs_path: str='not defined'
    qmax: float=np.inf
    qmin: float=0
    color: str=None
    marker: str=None
    linestyle: str=''
    shift: float=1

    imagepath: str=''
    waxs_imagepath: str=''

    def __post_init__(self):
        if self.angular_unit=='A':
            self.str_angular_unit ='$\\mathrm{\\AA^{-1}}$'
        if self.angular_unit=='nm':
            self.str_angular_unit ='$\\mathrm{nm^{-1}}$'

    def get_filename(self):
        return join(self.path, self.filename)

    def get_data(self, engine='pandas', **kwargs):
        """
        Reads the data file at self.path / self.filename
        Alternatively, a path can be provided as a keyword argument.
        """
        qmax = np.inf
        if 'qmax' in kwargs:
            qmax = kwargs.pop('qmax')
        qmin = 0
        if 'qmin' in kwargs:
            qmax = kwargs.pop('qmin')
        path = self.get_filename()
        if 'path' in kwargs:
            path = kwargs.pop('path')
        scale = 1
        if 'scale' in kwargs:
            scale = kwargs.pop('scale')
        self.data1d = sasmodels_data.load_data(path, **kwargs)
        if engine == 'sasview':
            return self.data1d
        elif engine == 'pandas':
            outdict = {'q': self.data1d.x,
                       'I': self.data1d.y*scale,
                       'err_I': self.data1d.dy*scale}
        df = pd.DataFrame(outdict)
        df = df[df.q<self.qmax]
        df = df[df.q<qmax]
        df = df[df.q>self.qmin]
        df = df[df.q>qmin]
        return df

    def save_data(self, filepath, df=None, **kwargs):
        if df is None:
            df = self.get_data(engine='pandas', **kwargs)
        df.to_csv(filepath, sep='\t', float_format='{:.7e}'.format, index=False)

    def get_filename_sasImage(self):
        with open(self.get_filename(), 'r') as f:
            lines = f.readlines()
        matched_lines = [line for line in lines if re.search('Sample filename', line)]
        filename = matched_lines[0].split()[0]
        return join(self.rawpath, filename)
        return self.imagepath

    def get_sasImage(self, filename=None, data=True, frame=None):
        if filename is None:
            if data:
                print('test', self.get_filename_sasImage())
                return fabio.open(self.get_filename_sasImage(), frame=frame).data
            else:
                return fabio.open(self.get_filename_sasImage(frame=frame))
        else:
            if data:
                return fabio.open(filename, frame=frame).data
            else:
                return fabio.open(filename, frame=frame)

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
                 sep=r'\s+', names=['q', 'I', 'err_I'], skiprows=skiprows,
                 nrows = nrows,
                 dtype=np.float64)
        return df

    def pyfai_integrate1d(self):
        masked_image = self.get_masked()
        nbins=200.
        sdd = self.sasImage.detector_distance[0]
        centerx, centery = masked_image.center
        pixelsizex, pixelsizey = masked_image.pixel_size
        ai = pyFAI.azimuthalIntegrator.AzimuthalIntegrator(dist=sdd,
                                                           poni1=centerx*pixelsizex,
                                                           poni2=centery*pixelsizey,
                                                           detector='pilatus300k')
        # ai.setFit2D(sdd, centerx, centery)
        ai.wavelength = masked_image.wavelength[0] * 10**-10

        q, I, err_I = ai.integrate1d(data=masked_image.data, npt=nbins,
                                     unit="q_nm^-1", mask=masked_image.mask,
                                     error_model='poisson', correctSolidAngle=True)
        dict={'q':q, 'I':I/self.exposure_time, 'err_I':err_I/self.exposure_time}
        df = pd.DataFrame(dict)
        return df

    def show_sasImage(self, filename=None, ax=None, vmin=1,
                      vmax=1000, frame=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots()
        img = self.get_sasImage(filename=filename, frame=frame)
        # img = np.ma.masked_less(img, 200)
        img = np.ma.masked_greater(img, 10000)
        cmap = cm.inferno
        cmap.set_bad(color='red')  # Set masked values to red
        im = ax.matshow(img, cmap=cmap, origin='lower',
                        norm=LogNorm(vmin=vmin, vmax=vmax), **kwargs)
        ax.grid(False)
        return _colorbar(im)

    def smooth_data(self, file_path=None, out_path=None):
        if file_path is None:
            file_path = self.get_filename()
        if out_path is None:
            out_path = self.get_workdir()
        script_path = '/home/johannes/pCloudDrive/gromacs/odna_cufix/ac6/smooth-saxs-curve.sh'
        result = subprocess.run(['bash', script_path, '-f', file_path],
                                cwd=out_path)
        return result

    def get_waxs_spectra(self, filename, include1=False, every_n=1, max_out=np.inf):
        out = {'dfs': [],
               'times': []}
        with open(filename) as f:
            start_row = 0
            time = 0
            plot = 0 if include1 else 1
            for i, line in enumerate(f):
                if len(line.split('simulation step')) > 1:
                    time = float(line.split('step')[1])
                if line == '@type xydy\n':
                    start_row = i+1
                if line == '&\n':
                    end_row = i-2
                    if start_row>0 and (plot%every_n)==0 and len(out['dfs'])<max_out:
                        df = pd.read_csv(filename, skiprows=start_row, nrows=end_row-start_row, sep=r'\s+',
                                        names=['q', 'I', 'err_I', 'I2', 'I3', 'I4'])
                        out['dfs'].append(df)
                        out['times'].append(time)
                    # time += 1
                    plot += 1
            df = pd.read_csv(filename, skiprows=end_row+5, sep=r'\s+',
                            names=['q', 'I', 'err_I', 'I2', 'I3', 'I4'])
            out['dfs'].append(df)
            out['times'].append(time)
        try:
            df_final = self.get_gromacs(path='', filename=filename.replace('spectra.xvg', 'final.xvg'))
            out['df_final'] = df_final
            print('added final')
        except:
            out['df_final'] = None
            print('no _final found')
        return out


    def plot_data(self, label=None, scale=1, buf=False,
                  linear_shift=0, **kwargs):
        df = self.get_data()
        if buf==True:
            df = self.get_averaged(buf=True)
        if 'figure' in kwargs:
            fig, ax = kwargs['figure']
            kwargs.pop('figure')
        if 'ax' in kwargs:
            ax = kwargs['ax']
            kwargs.pop('ax')
        else:
            fig, ax = plt.subplots()
            ax.set_xlabel('$q$ [1/nm]')
            ax.set_ylabel('$I$ [1/cm]')
        if 'shift' in kwargs:
            scale = kwargs['shift']
            kwargs.pop('shift')
        if 'marker' in kwargs:
            marker=kwargs['marker']
            kwargs.pop('marker')
        else:
            marker = self.marker
        if 'linestyle' in kwargs:
            linestyle=kwargs['linestyle']
            kwargs.pop('linestyle')
        else:
            linestyle = self.linestyle
        if 'color' not in kwargs:
            kwargs['color'] = self.color
        unit = self.angular_unit
        if 'unit' in kwargs:
            unit = kwargs.pop('unit')
        str_unit = 'not assigned'
        if unit == 'nm':
            str_unit ='$\\mathrm{nm^{-1}}$'
            if self.angular_unit == 'A':
                df.q = df.q * 10
        if unit == 'A':
            str_unit ='$\\mathrm{\\AA^{-1}}$'
            if self.angular_unit == 'nm':
                df.q = df.q / 10
        every_n = 1
        if 'every_n' in kwargs:
            every_n = kwargs.pop('every_n')
        df = df[::every_n]

        markers, caps, bars = ax.errorbar(df.q, (df.I+linear_shift)*scale, df.err_I*scale,
                                          marker=marker, linestyle=linestyle,
                                          label=label, elinewidth=0.2,  **kwargs)
        # [bar.set_alpha(0.2) for bar in bars]

        ax.set_xlabel(f'q [{str_unit}]')
        ax.set_ylabel('Intensity [cm$^{-1}$]')
        ax.set_xscale('log')
        ax.set_yscale('log')
        return ax

    def plot_kratky(self, label=None, scale=1, buf=False, **kwargs):
        df = self.get_data()
        if buf==True:
            df = self.get_averaged(buf=True)
        if 'figure' in kwargs:
            fig, ax = kwargs['figure']
            kwargs.pop('figure')
        if 'ax' in kwargs:
            ax = kwargs['ax']
            kwargs.pop('ax')
        else:
            fig, ax = plt.subplots()
            ax.set_xlabel('$q$ [1/nm]')
            ax.set_ylabel('$I$ [1/cm]')
        if 'shift' in kwargs:
            scale = kwargs['shift']
            kwargs.pop('shift')
        else:
            shift=self.shift
        if 'marker' in kwargs:
            marker=kwargs['marker']
            kwargs.pop('marker')
        else:
            marker = self.marker
        if 'linestyle' in kwargs:
            linestyle=kwargs['linestyle']
            kwargs.pop('linestyle')
        else:
            linestyle=''
        if 'color' not in kwargs:
            kwargs['color'] = self.color

        markers, caps, bars = ax.errorbar(df.q, df.I*df.q*df.q*scale, df.err_I*scale,
                                          marker=marker, linestyle=linestyle,
                                          label=label, elinewidth=0.2,  **kwargs)
        # [bar.set_alpha(0.2) for bar in bars]

        ax.set_xlabel(f'q [{self.str_angular_unit}]')
        ax.set_ylabel('Intensity $\\cdot q^2$ [A.U.]')
        return ax

    def get_workdir(self):
        dirname = f'./workdirs/{self.name}'
        return dirname

    def make_workdir(self, dirname=None):
        if dirname is None:
            dirname = self.get_workdir()
        osu.create_path('workdirs')
        osu.create_path(dirname)
        shutil.copyfile(self.get_filename(), join(dirname, 'data.dat'))

    @staticmethod
    def get_crysol(path):
        df = pd.read_csv(path, skiprows=2, sep=r'\s+',
                         names=['q', 'I', 'err_I', 'fit'])
        return df

    @staticmethod
    def get_chi2(df, df_fit):
        if 'fit' in df_fit.keys():
            fit = df_fit.fit
        else:
            fit = df_fit.I
            # print('df_fit.I used because no df_fit.fit was found')
        data = df.I
        err_data = df.err_I
        chi2 = np.sum(((data - fit) / err_data)**2 /
                             (len(df) - 1))
        return df_fit, chi2

    def scale_and_offset_fit(self, df_ref, df_scale, p0=None,
                             scale=None, offset=None,
                             bounds=(-np.inf, np.inf)):
        """Fit scale and offset of y_original to match y_intact."""
        import numpy as np
        from scipy import optimize
        from scipy.interpolate import interp1d

        outdict = {}
        # Interpolate y_original onto the x-values of y_intact
        df_ref_out = df_ref[df_ref.q < np.max(df_scale.q)]
        df_ref_out = df_ref_out[df_ref_out.q > np.min(df_scale.q)]
        f = interp1d(df_scale.q, df_scale.I, kind='linear')
        df_scale_out = pd.DataFrame({'q': df_ref_out.q,
                                     'I': f(df_ref_out.q)})
        if 'err_I' in df_scale.keys():
            err_f = interp1d(df_scale.q, df_scale.err_I, kind='linear')
            df_scale_out = pd.DataFrame({'q': df_ref_out.q,
                                         'I': f(df_ref_out.q),
                                         'err_I': err_f(df_ref_out.q)})
        # Fit scale and offset
        def linear_model(y, scale, offset):
            """Linear model."""
            return scale * y + offset
        if scale is None and offset is None:
            popt, pcov = optimize.curve_fit(linear_model, df_scale_out.I, df_ref_out.I,
                                         sigma=df_ref_out.err_I, p0=p0, bounds=bounds)
            scale, offset = popt
            err_scale, err_offset = np.sqrt(np.diag(pcov))
        else:
            err_scale, err_offset = [None, None]
            print('using input scale and offset parameters')
        df_scale_out.I = scale * df_scale_out.I + offset
        if 'err_I' in df_scale_out.keys():
            df_scale_out.err_I = scale * df_scale_out.err_I
        df_scale_out['res'] = (df_ref_out['I'] - df_scale_out['I']) / df_ref_out.err_I
        outdict['df'] = df_scale_out
        outdict['chi2'] = self.get_chi2(df_ref_out, df_scale_out)[1]
        outdict['scale'] = scale
        outdict['offset'] = offset
        outdict['err_scale'] = err_scale
        outdict['err_offset'] = err_offset
        return outdict

    def bin_data(self, df, qs):
        from scipy.interpolate import interp1d
        f = interp1d(df.q, df.I, kind='linear')
        binned_df = pd.DataFrame({'q': qs,
                                  'I': f(qs)})
        if 'err_I' in df.keys():
            err_f = interp1d(df.q, df.err_I, kind='linear')
            binned_df = pd.DataFrame({'q': qs,
                                     'I': f(qs),
                                     'err_I': err_f(qs)})
        return binned_df


    def get_rg(self, df=None, qmin=0, qmax=10,
               plot=False, Rg0=1, I00=None, ax=None,
               bounds=((0,0), (np.inf, np.inf)),
               plot_kwargs={}):
        outdict = {}
        if I00 is None:
            I00 = df.I.max()
        if df is None:
            df = self.get_data()
        df = df[df.q > qmin]
        df = df[df.q < qmax]
        df['lnI'] = np.log(df.I)
        def guinier(q, Rg, I0):
            I = I0 * np.exp(-(Rg*q)**2/3)
            return I
        popt, pcov = optimize.curve_fit(guinier, df.q, df.I, sigma=df.err_I, p0=[Rg0, I00], bounds=bounds)
        outdict['Rg'] = popt[0]
        outdict['err_Rg'] = np.sqrt(pcov[0, 0])
        outdict['I0'] = popt[1]
        outdict['err_I0'] = np.sqrt(pcov[1, 1])
        if plot:
            if ax is None:
                ax = self.plot_data(**plot_kwargs)
            ax.errorbar(df.q, guinier(df.q, *popt), marker='', linestyle='-.', label='fit')
        return outdict

def gen_guinier_fitfunc(alpha):
    def inner(q, Rg, A):
        if alpha == 0:
            pre = 1
        elif alpha == 1 or alpha == 2:
            pre = alpha * np.pi * q ** -alpha
        else:
            raise TypeError('alpha needs to be in 0,1,2')

        I = pre * A * np.exp(-Rg ** 2 * q ** 2 / (3 - alpha))
        return I
    return inner

def guinier_porod_3D(q, Rg1, s1, Rg2, s2, G2, dd):
    q = np.atleast_1d(q)

    # define parameters for smooth transitions
    Q1 = (1 / Rg1) * ((dd - s1) * (3 - s1) / 2) ** 0.5
    Q2 = ((s1 - s2) / (2 / (3 - s2) * Rg2 ** 2 - 2 / (3 - s1) * Rg1 ** 2)) ** 0.5
    G1 = G2 / (np.exp(-Q2 ** 2 * (Rg1 ** 2 / (3 - s1) -
                                  Rg2 ** 2 / (3 - s2))) * Q2 ** (s2 - s1))
    D = G1 * np.exp(-Q1 ** 2 * Rg1 ** 2 / (3 - s1)) * Q1 ** (dd - s1)

    # define functions in different regions
    def _I1_3regions(q):
        res = G2 / q ** s2 * np.exp(-q ** 2 * Rg2 ** 2 / (3 - s2))
        return res

    def _I2_3regions(q):
        res = G1 / q ** s1 * np.exp(-q ** 2 * Rg1 ** 2 / (3 - s1))
        return res

    def _I3_3regions(q):
        res = D / q ** dd
        return res

    I = np.piecewise(q, [q < Q2, (Q2 <= q) & (q < Q1), q >= Q1],
                     [_I1_3regions, _I2_3regions, _I3_3regions])
    return I

