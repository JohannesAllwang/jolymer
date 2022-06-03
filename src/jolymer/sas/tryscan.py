#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 15:21:55 2020

@author: johannes
"""

import numpy as np
import pandas as pd
import os
from os.path import join

from .. import database_operations as dbo
from .. import Sample as Sample
from .. import os_utility as osu
from .. import plot_utility as plu

from .SAXS_Measurement import SAXS_Measurement
from .ms import Ms
from .desy import Desy

class TryscanSample:

    def __init__(self, PS, NaOH, pH, tt):
        self.PS = PS
        self.NaOH = NaOH * 0.5e-3 # Molar
        self.pH = pH
        self.tt = tt

sampledict = {
    'hau' : ['HA', 0, 4.17, False],
    'hat0' : ['HA', 0, 4.17, True],
    'hat1' : ['HA', 11, 4.85, True],
    'hat2' : ['HA', 22, 6.57, True],
    'hat3' : ['HA', 33, 7.09, True],
    'hat4' : ['HA', 44, 10.07, True],
    'hat5' : ['HA', 55, 10.81, True],
    'hat6' : ['HA', 66, 10.95, True],
    'hat7' : ['HA', 77, 11.11, True],
    'hat8' : ['HA', 88, 11.36, True],
    'csu' : ['CS', 0, 4.30, False],
    'cst0' : ['CS', 0, 4.30, True],
    'cst1' : ['CS', 11, 5.24, True],
    'cst2' : ['CS', 22, 6.38, True],
    'cst3' : ['CS', 33, 7.96, True],
    'cst4' : ['CS', 44, 10.06, True],
    'cst5' : ['CS', 55, 10.88, True],
    'cst6' : ['CS', 66, 11.17, True],
    'cst7' : ['CS', 77, 11.19, True],
    'cst8' : ['CS', 88, 11.25, True],
    'tru' : [None, 0, 4.17, False]
}


class Tryscan(Desy):

    def __init__(self, tryscanid, issetup=False, cout=True, **kwargs):
        with dbo.dbopen() as c:
            query = f"""SELECT id, filename, sample, comment
            FROM tryscan_measurements WHERE id=?"""
            c.execute(query, (tryscanid,))
            self.datestring = '20220525'
            self.id, self.filenames, self.samplestring, self.comment, = list(c.fetchone())
            self.filenames = self.filenames.split(';')
            self.filename = self.filenames[0]
        if cout:
            print(self.samplestring)
        self.sample = TryscanSample(*sampledict[self.samplestring])
        # self.path = os.path.join(self.rawdatapath, f'desy{self.datestring}/')
        self.path = join(self.rawdatapath, 'tryscan', self.filename)
        self.processed_subtracted_file = os.path.join(self.path, 'processed_subtracted.dat')
        self.processed_file = os.path.join(self.path, 'processed.dat')
        self.analysis_path = join(self.path, 'analysis')
        self.averaged_path = join(self.path, 'averaged')
        self.subtracted_path = join(self.path, 'subtracted')
        self.frames_path = os.path.join(self.path, 'frames')
        self.absolute_path = os.path.join(self.path, 'absolute')
        self.averaged_path = os.path.join(self.path, 'averaged')
        self.buffer_frames_path = os.path.join(self.path,'buffer', 'frames')
        self.buffer_absolute_path = os.path.join(self.path,'buffer', 'absolute')
        self.origpath = os.path.join(self.rawdatapath, f'desy{self.datestring}', 'datacollection', 'data', 'absolute')
        # get the path of the


    def get_filename(self):
        return os.path.join(self.path, 'processed_subtracted.dat')

    def get_data(self, cout=True, altername='No', nrows=2640):
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


class TryscanJumpPart(Desy):

    def __init__(self, jumpname, T, time):
        self.datestring = '20220525'
        self.id = f'{jumpname}-{T}-{time}'
        self.samplestring = f'{jumpname[0:2]}u'
        self.path = join(self.rawdatapath, 'tryscan', jumpname)
        self.time = time
        self.filename = join(T, 'subtracted', f't{time}.dat')
        self.sample = TryscanSample(*sampledict[self.samplestring])
        # self.path = os.path.join(self.rawdatapath, f'desy{self.datestring}/')
        self.processed_subtracted_file = os.path.join(self.path, 'processed_subtracted.dat')
        self.processed_file = os.path.join(self.path, 'processed.dat')
        self.analysis_path = join(self.path, 'analysis')
        self.averaged_path = join(self.path, 'averaged')
        self.subtracted_path = join(self.path, 'subtracted')
        self.frames_path = os.path.join(self.path, 'frames')
        self.absolute_path = os.path.join(self.path, 'absolute')
        self.averaged_path = os.path.join(self.path, 'averaged')
        self.buffer_frames_path = os.path.join(self.path,'buffer', 'frames')
        self.buffer_absolute_path = os.path.join(self.path,'buffer', 'absolute')
        self.origpath = os.path.join(self.rawdatapath, f'desy{self.datestring}', 'datacollection', 'data', 'absolute')

    def get_filename(self):
        return os.path.join(self.path, self.filename)

    def get_data(self, cout=True, altername='No', nrows=2400):
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

times = 75 + 45*np.linspace(0, 30, num=31)
csT20 = [TryscanJumpPart('cstryjump', '20C', int(time)) for time in times][0:-2]
csT60 = [TryscanJumpPart('cstryjump', '60C', int(time)) for time in times][0:-2]
haT20 = [TryscanJumpPart('hatryjump', '20C', int(time)) for time in times][0:-2]
haT60 = [TryscanJumpPart('hatryjump', '60C', int(time)) for time in times][0:-2]
tryT20 = [TryscanJumpPart('tryjump', '20C', int(time)) for time in times][0:-2]
tryT60 = [TryscanJumpPart('tryjump', '60C', int(time)) for time in times][0:-2]
colors = list(plu.colorgradient('colors', [plu.tum_lblue, plu.tum_dred], 31))
colorbar = plu.colorgradient('colors', [plu.tum_dblue, 'white', plu.tum_dred], 31)
colors = list(colorbar)

for csT in [csT20, csT60, haT20, haT60, tryT20, tryT60]:
    for i, m in enumerate(csT):
        m.model = None
        m.iqmin = 0
        m.iqmax = 2300
        m.dataqmin = 0
        m.dataqmax = 2300
        m.color = colors[i]
        m.label = f'{m.time} sec'
        m.marker = 'o'
        m.model = None
        m.bounds = None
        m.p0 = None
        m.fixed_pars = None

csT20 = Ms(csT20)
csT60 = Ms(csT60)
haT20 = Ms(haT20)
haT60 = Ms(haT60)
tryT20 = Ms(tryT20)
tryT60 = Ms(tryT60)
