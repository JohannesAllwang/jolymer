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

    def get_filename(self):
        return join(self.path, self.filename)

    def get_spectra_filename(self):
        return join(self.path, self.spectra_name)

    def get_data(self, engine='pandas', **kwargs):
        """
        Reads the data file at self.path / self.filename
        Alternatively, a path can be provided as a keyword argument.
        """
        return self.read_xvg(self.get_filename(), names=['q', 'I', 'err_I'],
                             include1=True, onlyfirst=True)

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
                 sep='\s+', names=['q', 'I', 'errI'], skiprows=skiprows,
                 nrows = nrows,
                 dtype=np.float64)
        return df

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
                if line == '@type xydy\n':
                    start_row = i+1
                if line == '&\n':
                    end_row = i-2
                    if start_row>0 and (plot%every_n)==0 and len(out['dfs'])<max_out:
                        df = pd.read_csv(filename, skiprows=start_row, nrows=end_row-start_row, sep='\s+',
                                        names=names)
                        if onlyfirst:
                            return df
                        out['dfs'].append(df)
                        out['times'].append(time)
                    time += 1
                    plot += 1
            df = pd.read_csv(filename, skiprows=end_row+5, sep='\s+',
                            names=names)
            out['dfs'].append(df)
            out['times'].append(time)
        return out


    @staticmethod
    def get_crysol(path):
        df = pd.read_csv(path, skiprows=2, sep='\s+',
                         names=['q', 'I', 'errI', 'fit'])
        return df

    @staticmethod
    def get_chi2(df, df_fit):
        if len(df) != len(df_fit):
            print('Interpolating since df and df_fit donnot have the same len')
            from scipy.interpolate import interp1d
            df_fit = df_fit[df_fit.q<df.q.max()]
            print(df.q.max(), df_fit.q.max())
            interp_function = interp1d(df.q, df.I)
            df_fit['data'] = interp_function(df_fit.q)
            err_interp_function = interp1d(df.q, df.err_I)
            df_fit['err_data'] = err_interp_function(df_fit.q)
            df_fit = df_fit[df_fit.err_data>0]
        else:
            df_fit['data'] = df.I
            df_fit['err_data'] = df.err_I
        chi2 = np.sum(((df_fit.data - df_fit.fit) / df_fit.err_data)**2 /
                             (len(df_fit) - 1))
        return df_fit, chi2
