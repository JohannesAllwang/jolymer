#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 16:52:47 2020

@author: johannes
"""

from .. import database_operations as dbo
import pandas as pd
import numpy as np
import datetime as dt
import os
import matplotlib.pyplot as plt
from scipy import optimize, constants
from .. import Sample
from ..Measurement import Measurement


class lsi(Measurement):

    
    def __init__(self, lsi_id, instrument = 'lsi', evaluate_script = True):
        self.instrument = instrument
        with dbo.dbopen() as c:
            c.execute(f"SELECT * FROM {self.instrument}_measurements WHERE id=?", (lsi_id,))
            self.id, self.datestring, self.filename, self.samplestring, self.mode, self.comment = list(c.fetchone())
            c.execute(f"SELECT * FROM {self.instrument}_exceptions WHERE {self.instrument}_id=?", (lsi_id,))
            exceptions = list(c.fetchall())
            self.exceptions = [e[1] for e in exceptions]

        self.figures_path = os.path.join(Measurement.figures_path, f'{self.instrument}_{self.id}')

        if self.samplestring == None:
            self.sample = None
        else:
            sample_type, sample_id = self.samplestring.split('_')
            self.sample = Sample.get_sample(sample_type, sample_id)

        self.path =  os.path.join(self.rawdatapath, 'dls', f'{self.datestring}')
        if evaluate_script:
            script_df = pd.read_csv(f'{self.path}/{self.filename}.lsiscript', sep='\t', names=['index', 'start_angle', 'end_angle'
                                                                ,'step_angle', 'per_angle', 'sec', 'TODO1',
                                                                'TODO2', 'TODO3', 'TC'], index_col=0)
            Tlist = [x for x in set(script_df.TC)]
            Tlist.sort()
            self.TCs = Tlist
            self.smin = 1
            self.Tdict = {}

            seq_number = self.smin
            print('Use function self.get_scriptrow(i) with the followning:')
            for index, row in script_df.iterrows():
                print(self.samplestring)
                angles = range(int(row.start_angle), int(row.end_angle+1), int(row.step_angle))
                angles = np.array(angles)
                seq_numbers = list(range(seq_number, int(row.per_angle)*len(angles)+1 + seq_number))
                # print(seq_number, row.per_angle)
                # print(seq_numbers)
                oneT = _oneT(row.TC, angles, seq_numbers, row.per_angle, self.id, index, instrument=self.instrument)
                self.Tdict[index] = oneT
                seq_number += int(row.per_angle)*len(angles)
                print(f'({index}) : {row.TC}')
            self.smax = seq_number

            self.seq_numbers = range(int(self.smin), int(self.smax)+1)
            # self.exceptions = [float(e) for e in exceptions]

        self.rmin = 2e-9
        self.rmax = 2_000e-9
        self.filename += '_#.dat'
        if os.path.exists(self.rawdata_filename('0')):
            self.seq_numbers = [x-1 for x in self.seq_numbers]
        self.name = f'{self.instrument}{self.id}'
        # return self.Tdict

    def get_filename(self, seq_number):
        return self.filename.replace('#', str(seq_number))

    def get_fitpath(self, model):
        out = os.path.join(self.path, f'{model.name}_{self.instrument}_{self.id}')
        return out

    def get_fitfile(self, model, seq_number):
        out = os.path.join(self.get_fitpath(model),\
                self.get_filename(seq_number))
        return out

    def get_scriptrow(self, i):
        return self.Tdict[i]

    def get_metadata(self, seq, name):
        with open(self.rawdata_filename(seq)) as f:
            for line in f:
                if line.split(':')[0] == name:
                    return line.split(':')[1]

    def get_TK(self):
        return self.TC + 273.15

    def get_wl(self):
        wl = self.get_metadata(1, "Wavelength (nm)")
        wl = float(wl) * 10**-9
        return wl

    def rawdata_filename(self, seq_number):
        out = os.path.join(self.path,  self.get_filename(seq_number))
        return out

    def get_exceptions(self):
        query = f"""
        SELECT seq_number FROM {self.instrument}_exceptions
        WHERE {self.instrument}_id = {self.id}
        """
        out = dbo.execute(query)
        return out

    def add_exception(self, seq_number, reason=''):
        query = f"""
        INSERT INTO {self.instrument}_exceptions
        VALUES ({self.id}, {seq_number}, '{reason}');
        """
        self.exceptions.append(seq_number)
        out = dbo.execute(query)
        return out

    def add_exceptions(self, seq_numbers, reason=''):
        for s in seq_numbers:
            self.add_exception(s, reason=reason)
        return dbo.get_table('{self.instrument}_exceptions')

    def get_data(self, seq_number, xmin=15, xmax=221):
        if self.mode=='mod3d':
            xmax=170
        filename = self.rawdata_filename(seq_number)
        df = pd.read_csv(filename, sep='\s+', skiprows=16+xmin, nrows=xmax-xmin,
                         header=None, names=['t', 'g2'], engine='python')
        return df

    def get_summary(self):
        filename = self.rawdata_filename('')
        # Scattering angle   Mean_CR * sin(angle) / Laser Intensity (kHz/mW)    
        # g2(t=0)-1    CR CHA (kHz)    CR CHB (kHz)    Temperature (K)    Laser intensity (mW) Hydrodynamic Radius (nm) Width (nm)
        df = pd.read_csv(filename, sep='\t', skiprows=1, header=None, names=['angle', 'I', 'g20', 'CRA', 'CRB', 'TK', 'I0', 'rh', 'width'])
        # Sequence number is not in the file so we find it from the index:
        if os.path.exists(self.rawdata_filename('0')):
            df['seq_number'] = range(len(df))
        else:
            df['seq_number'] = range(1, len(df)+1)
        # select only rows that correspond to self
        mask = df.seq_number.apply(lambda x: any(item for item in self.seq_numbers if item == x))
        df = df[mask]
        return df

    def get_fitpar(self, fit, parameter):
        pass

    def get_gammastar(self):
        pass
    
    def phirange(self, phi):
        "returns all the seq_numbers of the argument angle. Exceptions excepted."
        index = np.where(self.angles == phi)[0]
        first = int(self.per_angle * index + self.seq_numbers[0])
        delta = int(self.per_angle)
        out = []
        for n in range(first, first + delta):
            if n in self.exceptions:
                continue
            out.append(n)
        return out

    def seqinphi(self, seq_number, phi):
        out = seq_number in self.phirange(phi)
        return out

    def phifromseq(self, seq_number):
        for phi in self.angles:
            if self.seqinphi(seq_number, phi):
                return phi
    
    def RfromD(self, D):
        TK = self.get_TK()
        visc = self.sample.get_viscosity(TK)
        out = constants.Boltzmann*TK / (6*np.pi*visc*D)
        return out

    def RfromG(self, Gamma, qq):
        D = Gamma / qq
        out = self.RfromD(D)
        return out

    def Rfromt(self, t, qq):
        D = 1 / (qq * t)
        out = self.RfromD(D)
        return out
    
    def DfromR(self, rh):
        TK = self.get_TK()
        out = constants.Boltzmann*TK / \
                (6*np.pi*self.sample.get_viscosity(TK)*rh)
        return out

    def q(self, phi):
        "calculates the scattering vector q[m^-1] from the scattering angle 2\Theta."
        wl = self.get_wl()
        n = self.sample.get_n(self.get_TK(), wl)
        out = n *4* np.pi *np.sin((phi*np.pi/360))/wl
        return out
    
    def qq(self, phi):
        "calculates the square of scattering vector q^2[m^-2] from the scattering angle 2\Theta."
        return self.q(phi)**2

    def get_visc(self):
        return self.sample.get_viscosity(self.get_TK())
    
    def plot_data(self, seq_number, ax=None, **kwargs):
        if ax==None:
            fig, ax = plt.subplots()
            ax.set_xlabel('$\\tau$ [s]')
            ax.set_ylabel('$g_2-1$')
        if 'marker' not in kwargs:
            kwargs['marker']='.'
        if 'linestyle' not in kwargs:
            kwargs['linestyle']=''
        df = self.get_data(seq_number)

        # Autocorrelation Plots:
        ax.plot(df.t, df.g2 , **kwargs)
        ax.set_xscale('log')
        ax.set_ylim(0,1)
        if self.mode == '3dcross':
            ax.set_ylim(0, 0.25)
        if self.mode == 'mod3d':
            ax.set_ylim(0, 0.85)
        return ax
    
    def plot_fit(self, seq_number, fit, fitcolor='lightsalmon', ax = None, **kwargs):
        ax = self.plot_data(seq_number, ax=ax, **kwargs)
        axres = ax.twinx()
        fitdf = fit.get_fit(self, seq_number)
        fitdf = fitdf[fitdf.fit<1]
        ax.plot(fitdf.t, fitdf.fit, '-', color=fitcolor)
        rescolor=fitcolor
        if 'color' in kwargs:
            rescolor=kwargs['color']
        axres.plot(fitdf.t, fitdf.res, alpha=0.2, color=rescolor)
        axres.set_ylabel('Relative Residuals')
        axres.set_ylim(-0.05, 0.05)
        return ax

    def plot_dist(self, seq_number, fit, ax=None, xspace='t', **kwargs):
        "plot distribution, if exists, in equal area representation."
        if ax==None:
            fig, ax = plt.subplots()
        dfd = fit.get_dist(self, seq_number)
        if xspace=='rh':
            # rapp in nm as xaxis
            ax.plot(dfd.rh * 10**9, dfd.dist, **kwargs)
        elif xspace=='t':
            ax.plot(dfd.t, dfd.dist, **kwargs)
        ax.set_xscale('log')
        return ax
    
    def qplot(self, fit, parameters, function=lambda x:x, 
              err_function=lambda x,err_x:err_x, ax=None, **kwargs):
        if ax==None:
            fig, ax = plt.subplots()
            ax.set_ylabel(fit.pardict[parameters[0]][0])
            ax.set_xlabel("$q^2$ $\\mathrm{[\\mu m^{-2}]}$")
        if 'marker' not in kwargs:
            kwargs['marker']='o'
        if 'linestyle' not in kwargs:
            kwargs['linestyle']=''
    
        df = fit.get_phitable(self)
        xdata=df.qq
        args = [df[par] for par in parameters]
        ydata = function(*args)
        err_args = [df[f'err_{par}'] for par in parameters]
        err_ydata = err_function(*args, *err_args)
        ax.errorbar(xdata, ydata, err_ydata, **kwargs)
        ax.set_xlim(left=0)
        return ax

class lsi3d(lsi):

    instrument = 'lsi3d'
    
    def __init__(self, lsi_id, evaluate_script=True):
        lsi.__init__(self, lsi_id, instrument='lsi3d', evaluate_script=evaluate_script)


class _oneT(lsi):

    def __init__(self, TC, angles, seq_numbers, per_angle, lsi_id, script_row, instrument='lsi'):
        self.seq_numbers= seq_numbers
        if instrument == 'lsi':
            lsi.__init__(self, lsi_id, evaluate_script=False)
        elif instrument == 'lsi3d':
            lsi3d.__init__(self, lsi_id, evaluate_script=False)
        self.TC = TC
        self.angles = angles
        self.per_angle = per_angle
        self.script_row = script_row
        self.name = f'{self.instrument}{self.id}_{int(script_row)}'
