"""
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from os.path import join
from dataclasses import dataclass

from .. import database_operations as dbo
from .. import Sample as Sample
from .. import os_utility as osu

from .SAXS_Measurement import SAXS_Measurement

@dataclass
class TPS13A(SAXS_Measurement):

    instrument = 'tps13a'

    def get_parameter(self, parstring):
        with open(self.get_filename()) as f:
            for line in f:
                if line.split(':')[0] == parstring:
                    return line.split(':')[1]

    def get_TC(self):
        """
        gets the temperature in Celsius from the processed_subtracted file.
        """
        return float(self.get_parameter('SC Target Temperature'))

    def get_notparent_fullpaths(self):
        path = self.absolute_path
        parents = os.listdir(path)
        nums = []
        for parent in parents:
            num = parent.split('_')[1]
            if not num in nums:
                nums.append(num)
        path = self.origpath
        filenames = []
        for num in nums:
            filenames += [fil for fil in os.listdir(self.origpath) if fil.split('_')[1]==num]
        filenames = [os.path.join(self.origpath, fil) for fil in filenames if not fil in parents]
        return filenames

    def get_notparent_dfs(self):
        out = []
        files = []
        for path in self.get_notparent_fullpaths():
            df = pd.read_csv(path, skiprows=2, header=None, names=['q', 'I', 'err_I'],
                     delimiter='\s+', nrows = 2652)
            out.append(df)
            files.append(path)
        return out, files

    def get_absolute_fullpaths(self, buf=False):
        path=self.absolute_path
        if buf:
            path = self.buffer_absolute_path
        files = [join(path, fil) for fil in os.listdir(path) ]
        return files

    def get_frame_dfs(self, buf=False):
        out=[]
        path=self.frames_path
        if buf:
            path = self.buffer_frames_path
        for fil in os.listdir(path):
            df = pd.read_csv(os.path.join(path, fil), skiprows=2, header=None, names=['q', 'I', 'err_I'],
                     delimiter='\s+', nrows = 2652)
            out.append(df)
        return out

    def get_absolute_dfs(self, buf=False):
        out=[]
        path=self.absolute_path
        if buf:
            path = self.buffer_absolute_path
        files = os.listdir(path)
        for file in files:
            df = pd.read_csv(os.path.join(path, file), skiprows=2, header=None, names=['q', 'I', 'err_I'],
                     delimiter='\s+', nrows = 2652)
            out.append(df)
        return out, files

    def get_averaged_fullpath(self, buf=False):
        path=self.averaged_path
        files = os.listdir(path)
        sorb = 'buffer' if buf else 'sample'
        for fil in files:
            if len(fil.split(sorb)) >1:
                fullpath = join(path, fil)
        return fullpath

    def get_averaged(self, buf=False):
        path=self.averaged_path
        files = os.listdir(path)
        sorb = 'buffer' if buf else 'sample'
        for fil in files:
            if len(fil.split(sorb)) >1:
                df = pd.read_csv(os.path.join(path, fil), skiprows=3, header=None, names=['q', 'I', 'err_I'],
                     delimiter='\s+', nrows = 2652)
        return df

    def get_filename(self):
        return os.path.join(self.path, 'processed_subtracted.dat')

    def get_data(self, cout=True, altername='No', nrows=2653):
        "get data from data.csv and apply some filter."
        # path = os.path.join(self.path, 'data.csv')
        if altername=='No':
            path = self.get_filename()
        else:
            path = os.path.join(self.path, f'{altername}.dat')
        df = pd.read_csv(path, sep='\s+', header=None, skiprows=3, nrows=nrows, names=['q', 'I', 'err_I'])
        # df = pd.read_csv(filename, sep='\t', skiprows=16+xmin, nrows=xmax-xmin,
        #                  header=None, names=['t', 'g2'], engine='python')
        len_before = len(df)
        # df = df[df.I>0]
        len_after = len(df)
        if cout:
            print(f'{len_before-len_after} negative I values!')
        # df = df[df.I>df.err_I]
        # df = df[df.q<7.]
        len_after = len(df)
        if cout:
            print(f'{len_before-len_after} excluded I values!')
        return df

    def get_bumpsmodel(self, modelname):
        out = ''
        return out

    def plot_data(self, label=None, scale=1, buf=False, **kwargs):
        "plots data"
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
            marker='.'
        if 'linestyle' in kwargs:
            linestyle=kwargs['linestyle']
            kwargs.pop('linestyle')
        else:
            linestyle=''

        markers, caps, bars = ax.errorbar(df.q, df.I * scale, df.err_I * scale, marker = marker,
                    linestyle=linestyle, label = label, elinewidth=0.2,  **kwargs)
        # [bar.set_alpha(0.2) for bar in bars]

        ax.set_xscale('log')
        ax.set_yscale('log')
        return ax

    def plot_kratky(self, label=None, scale=1, buf=False, **kwargs):
        "plots kratky"
        df = self.get_data()
        df = df[df.q<2]
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
            ax.set_ylabel('$q^2I$')
        if 'shift' in kwargs:
            scale = kwargs['shift']
            kwargs.pop('shift')
        if 'marker' in kwargs:
            marker=kwargs['marker']
            kwargs.pop('marker')
        else:
            marker='.'
        if 'linestyle' in kwargs:
            linestyle=kwargs['linestyle']
            kwargs.pop('linestyle')
        else:
            linestyle=''

        xdata = df.q
        ydata = df.q * df.q * df.I
        erry = ydata * df.err_I/df.I
        ax.errorbar(xdata, ydata, yerr=erry, marker=marker, linestyle=linestyle, **kwargs)
        return ax
