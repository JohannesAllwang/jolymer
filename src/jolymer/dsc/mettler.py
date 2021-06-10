# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 12:45:08 2021

@author: xcill
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import shutil

from .. import database_operations as dbo
from .. import Sample as Sample
from .. import os_utility as osu
from ..Measurement import Measurement

class _isorun:

    def __init__(self, df, typen, hrate):
        self.df = df
        self.typen = typen
        self.type = typen[0:-1]
        self.n = typen[-1]
        self.hrate = hrate
        print('K per sec:', hrate)
        print('K per min:', hrate*60)

    def get_cp(self, gpm):
        df = self.df
        df['cp'] = df.hf / self.hrate
        df['cpmol'] = df.cp * gpm
        return df
        
    def plot_t(self, figure=None, pars=['hf', 'Tr'], **kwargs):
        if figure==None:
            fig, ax = plt.subplots()
        else:
            fig, ax = figure
        ax_Tr = ax.sharex()

    def plot_T(self, ax, **kwargs):
        df = self.df
        ax.plot(df.Tr, df.hf, **kwargs)
        return ax

    def plot_cpmol(self, ax, gpm, **kwargs):
        df = self.get_cp(gpm)
        ax.plot(df.Tr, df.cpmol, **kwargs)
        return ax

        
class Mettler(Measurement):

    dump_path = os.path.join(Measurement.rawdatapath, 'mettler_dump')

    def __init__(self, mettler_id, cout=True):
            
        with dbo.dbopen() as c:
            c.execute("SELECT * FROM mettler_measurements WHERE id=?;", (mettler_id,))
            self.id, self.datestring, self.samplestring,self.comment, self.timestring = list(c.fetchone())
                
        if cout:
            print(self.samplestring)
        if self.samplestring == None:
            self.sample = None
        else:
            sample_type, sample_id = self.samplestring.split('_')
            self.sample = Sample.get_sample(sample_type, sample_id)
        self.path = os.path.join(self.rawdatapath, 'mettler', 'mettler{0:03}/'.format(self.id))
        filenames = os.listdir(self.path) 
        for filename in filenames:
            if filename[-7::]=='{0:03}.txt'.format(self.id):
                self.filename = self.path + filename

    def get_data(self, fromto=(0,-1)):
        df = pd.read_csv(self.filename, skiprows=2, index_col=0, skipfooter=1,\
                 delimiter=r"\s+", header=None, names=['index', 't', 'hf', 'Ts', 'Tr'])[fromto[0]:fromto[1]]
        return df
    
    def split_data(self):
        """
        return df for each isothermal and dynamic scann
        """
        df = self.get_data()
        minT = np.min(df.Tr)
        maxT = np.max(df.Tr)
        Tchange = np.diff(df.Tr)
        ka = np.diff(Tchange)
        # Get values where the third derivative is non-zero
        # to find the start/stop of isothermals/dynamic scans
        ind = np.where(abs(ka)>=0.001)
        ind1 = list(ind[0])
        ind0 = [0] + ind1[0:-1]
        out = [[], [], []]
        isocount, hrcount, crcount = [1, 1, 1]
        # for loop over the inflection points:
        # i0 indexes the start of a isothermal/dynamic scan
        # i1 indexes the end or the beginning of the next one
        for i0, i1 in zip(ind0, ind1):
            out[0].append(df[i0:i1])
            # check the slope in the currently inspected trange
            if abs(Tchange[i0+5])<0.0001:
                out[1].append(f'iso{isocount}')
                isocount+=1
            elif Tchange[i0+5]>0.0001:
                out[1].append(f'hr{hrcount}')
                hrcount+=1
            elif Tchange[i0+5]< -0.0001:
                out[1].append(f'cr{crcount}')
                crcount += 1
            out[2].append(Tchange[i0+5])
        return out
    
    def get_iso_or_run(self, typen):
        """
        input sth like 'hr1' for the first heating run.
        return: isorun object
        """
        dfs, types, hrate = self.split_data()
        index = types.index(f'{typen}')
        out = _isorun(dfs[index], typen, hrate[index])
        return out

    def get_iso(self, n):
        return self.get_iso_or_run(f'iso{n}')  
 
    def get_cr(self, n):
        return self.get_iso_or_run(f'cr{n}')  
    
    def get_hr(self, n):
        return self.get_iso_or_run(f'hr{n}' )

    def plot_vst(self, fromto=(0,-1), ax=None, hfkwargs={'color' : 'tab:pink'}, Trkwargs={'color':'tab:orange'}):
        """lol I am lol"""
        df = self.get_data(fromto = fromto)
        if ax==None:
            fig, ax_hf = plt.subplots()
        else:
            ax_hf = ax
        ax_Tr = ax_hf.twinx()
        hfp, = ax_hf.plot(df.t, df.hf, **hfkwargs)
        Trp, = ax_Tr.plot(df.t, df.Tr, **Trkwargs)
        ax_hf.set_ylabel('Heat flow [Wg^-1]')
        ax_Tr.set_ylabel('Temperature [$^\\circ$C]')
        ax_hf.set_xlabel('Time [s]')
        ax_Tr.yaxis.label.set_color(Trp.get_color())
        ax_hf.yaxis.label.set_color(hfp.get_color())
        return ax_hf, ax_Tr

    def plot_T(self, ax=None, **kwargs):
        if ax==None:
            fig, ax = plt.subplots()
        else:
            ax=ax
        df = self.get_data()
        ax.plot(df.Tr, df.hf, **kwargs)
        ax.set_xlabel('$T$ [$^\circ$C]')
        ax.set_ylabel('Heat flow [Wg^-1]')
        return ax

    def plot_cpmol(self, ax=None, **kwargs):
        molpg = self.sample.get_molpg()
        hr = self.get_hr(2)
        df = hr.get_cp(1 / molpg)
        print(df)
        print('Calculates the hr to kcal/mol K')
        df.cpmol = df.cpmol * 0.000239006
        if ax==None:
            fig, ax = plt.subplots()
        else:
            ax=ax
        ax.plot(df.Tr, df.cpmol, **kwargs)
        ax.set_xlabel('$T$ [$^\circ$C]')
        ax.set_ylabel('$c_p$ [kcal/mol K]')
        return ax

def distribute(path):
    fnfs=os.listdir(path)
    for fnf in fnfs:
        if os.path.isdir(os.path.join(path, fnf)):
            distribute(os.path.join(path, fnf))
        else:
            mid = fnf[-7:-4]
            dist_dir = os.path.join(Mettler.rawdatapath, 'mettler', 'mettler'+mid)
            osu.create_path(dist_dir)
            shutil.move(os.path.join(path, fnf), os.path.join(dist_dir, fnf))

def distribute_dump():
    distribute(Mettler.dump_path)

def distribute_usb():
    usb_path = 'E:/DSC'
    for folder in os.listdir(usb_path):
        path = os.path.join(usb_path, folder)
        for f in os.listdir(path):
            mid = f[-7:-4]
            dist_dir = os.path.join(Mettler.rawdatapath, 'mettler', 'mettler'+mid)
            try:
                os.mkdir(dist_dir)
                shutil.copy(os.path.join(path, f), os.path.join(dist_dir, f))
                print('succesfully created', dist_dir)
            except:
                pass
