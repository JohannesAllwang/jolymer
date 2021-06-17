#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 15:21:55 2020

@author: johannes
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

from .. import database_operations as dbo
from .. import Sample as Sample
from .. import os_utility as osu

from .SAXS_Measurement import SAXS_Measurement

class Desy(SAXS_Measurement):

    @classmethod
    def setup(cls, values):
        """
        This method inserts the measurement into the database and creates the desired folderstructure.
        """
        # dbo.insert_values('desy_measurements', values)
        did = values[0]
        m = cls.__init__(m, did, issetup=False, count=False)
        osu.create_path(os.path.join(m.path,'buffer'))
        osu.create_path(m.frames_path)
        osu.create_path(m.absolute_path)
        osu.create_path(m.averaged_path)
        osu.create_path(m.buffer_frames_path)
        osu.create_path(m.buffer_absolute_path)
        unsorted_path = os.path.join(m.rawdatapath, f'desy{m.datestring}')
        osu.create_path(unsorted_path)
    
    def __init__(self, desy_id, issetup=False, cout=True, **kwargs):
        with dbo.dbopen() as c:
            c.execute("SELECT * FROM desy_measurements WHERE id=?;", (desy_id,))
            self.id, self.datestring, self.given_name, self.samplestring, _,\
                self.comment, self.timestring = list(c.fetchone())
        if cout:
            print(self.samplestring)
        if self.samplestring == None:
            self.sample = None
        else:
            sample_type, sample_id = self.samplestring.split('_')
            self.sample = Sample.get_sample(sample_type, sample_id)
        self.path = os.path.join(self.rawdatapath, 'desy{0:03}/'.format(self.id))
        self.processed_subtracted_file = os.path.join(self.path, 'processed_subtracted.dat')
        self.processed_file = os.path.join(self.path, 'processed.dat')
        self.frames_path = os.path.join(self.path, 'frames')
        self.absolute_path = os.path.join(self.path, 'absolute')
        self.averaged_path = os.path.join(self.path, 'averaged')
        self.buffer_frames_path = os.path.join(self.path,'buffer', 'frames')
        self.buffer_absolute_path = os.path.join(self.path,'buffer', 'absolute')

    def get_parameter(self, parstring):
        with open(self.get_filename()) as f:
            for line in f:
                if line.split(':')[0] == parstring:
                    return line.split(':')[1] 

    def get_TC(self):
        """
        gets the temperature in Celsius from the processed_subtracted file.
        """
        return get_parameter('SC Target Temperature')

    def get_frame_dfs(self, buf=False):
        out=[]
        path=self.frames_path
        if buf:
            path = self.buffer_frames_path
        for file in os.listdir(path):
            df = pd.read_csv(os.path.join(path, file), skiprows=2, header=None, names=['q', 'I', 'err_I'],
                     delimiter='   ', nrows = 2652)
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
                     delimiter='   ', nrows = 2652)
            out.append(df)
        return out, files

    def get_averaged(self, buf=False):
        path=self.averaged_path
        files = os.listdir(path)
        sorb = 'buffer' if buf else 'sample'
        for fil in files:
            if len(fil.split(sorb)) >1:
                df = pd.read_csv(os.path.join(path, fil), skiprows=3, header=None, names=['q', 'I', 'err_I'],
                     delimiter='   ', nrows = 2652)
        return df

    def get_filename(self):
        return os.path.join(self.path, 'processed_subtracted.dat')
        
    def get_data(self, cout=True):
        "get data from data.csv and apply some filter."
        # path = os.path.join(self.path, 'data.csv')
        path = self.get_filename()
        df = pd.read_csv(path, sep='  ', header=None, skiprows=3, nrows=2653, names=['q', 'I', 'err_I'])
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
    
    def plot_data(self, label=None, scale=1, **kwargs):
        "plots data"
        df = self.get_data()
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
    
