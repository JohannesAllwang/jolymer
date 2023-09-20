#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 16:52:47 2020

@author: johannes
"""

from .. import database_operations as dbo
import pandas as pd
import numpy as np
# import datetime as dt
import os
from os.path import join
import matplotlib.pyplot as plt
from scipy import constants
from .. import Sample
from ..Measurement import Measurement
from .. import os_utility as osu
from . import repes_utility as ru
from dataclasses import dataclass, field


@dataclass
class lsi(Measurement):

    color: str = '#000000'
    exceptions: list[int] = field(default_factory=list)

    def __post_init__(self):
        self.lol = 'lol'

    def __init__(self, path, filename,
                 TC, allangles, seq_numbers,
                 per_angle, sample=None,
                 lsid='X', samplestring=None,
                 datestring='XX.XX.XXXX',
                 exceptions=[],
                 badangles=[], mode='X',
                 instrument='lsi'):

        self.filename = filename
        self.path = path
        self.per_angle = per_angle
        self.TC = TC
        self.allangles = allangles
        self.angles = allangles
        self.seq_numbers = seq_numbers
        self.instrument = instrument
        self.id = lsid
        self.datestring = datestring
        self.samplestring = samplestring
        self.mode = mode
        self.comment = ''
        self.exceptions = exceptions
        self.badangles = badangles

        if sample is not None:
            self.sample = sample
        elif self.samplestring is None:
            self.sample = None
        else:
            print(samplestring)
            sample_type, sample_id = self.samplestring.split('_')
            self.sample = Sample.get_sample(sample_type, sample_id)

        self.phidls_path = join(self.path, 'phidls')
        osu.create_path(self.phidls_path)

        self.rmin = 2e-9
        self.rmax = 2_000e-9
        self.filename += '_#.dat'
        if os.path.exists(self.rawdata_filename('0')):
            self.seq_numbers = [x-1 for x in self.seq_numbers]
        self.name = f'{self.instrument}{self.id}'

    @classmethod
    def from_scriptfile(cls, scriptpath, scriptfile, scriptrow,
                        path=None, filename=None, smin=1,
                        **kwargs):
        if scriptfile.split('.')[-1] != 'lsiscript':
            scriptfile = f'{scriptfile}.lsiscript'
        if path is None:
            path = scriptpath
        if filename is None:
            filename = scriptfile[:-10]
        script_filepath = join(scriptpath, scriptfile)
        script_df = pd.read_csv(script_filepath, sep='\t',
                                names=[
                                    'index', 'start_angle',
                                    'end_angle', 'step_angle',
                                    'per_angle', 'sec', 'TODO1',
                                    'TODO2', 'TODO3', 'TC'],
                                index_col=0)
        seq_number = smin
        for index, row in script_df.iterrows():
            index = index-1
            allangles = range(int(row.start_angle),
                              int(row.end_angle+1),
                              int(row.step_angle))
            allangles = np.array(allangles)
            firstseq = seq_number
            lastseq = int(row.per_angle)*len(allangles)\
                + 1 + seq_number
            seq_numbers = list(range(firstseq, lastseq))
            seq_number += int(row.per_angle)*len(allangles)
            TC = row.TC
            per_angle = row.per_angle
            if index == scriptrow:
                break
        m = cls(path, filename,
                TC, allangles, seq_numbers,
                per_angle, **kwargs)
        return m

    @classmethod
    def from_db(cls, lsi_id, instrument='lsi',
                scriptrow=0, **kwargs):
        # self.instrument = instrument
        with dbo.dbopen() as c:
            query = f"""
            SELECT id, mdate, filename, sample, mode, comment
            FROM {instrument}_measurements WHERE id=?"""
            c.execute(query, (lsi_id,))
            id, datestring, filename, samplestring, mode,\
                comment = list(c.fetchone())
            c.execute(f"""
                       SELECT * FROM {instrument}_exceptions
                       WHERE {instrument}_id=?""",
                      (lsi_id,))
            exceptions = list(c.fetchall())
            exceptions = [e[1] for e in exceptions]
            c.execute(f"""
                       SELECT * FROM {instrument}_badangles
                       WHERE {instrument}_id=?""",
                      (lsi_id,))
            badangles = list(c.fetchall())
            badangles = [phi[1] for phi in badangles]

        path = os.path.join(cls.rawdatapath, 'dls', f'{datestring}')
        scriptfile = f'{path}/{filename}.lsiscript'
        m = cls.from_scriptfile(path, scriptfile, scriptrow,
                                lsid=lsi_id, mode=mode,
                                samplestring=samplestring,
                                datestring=datestring,
                                exceptions=exceptions,
                                badangles=badangles,
                                instrument=instrument, **kwargs)
        m.id = lsi_id
        m.datestring = datestring
        m.phidls_path = join(m.rawdatapath, 'phidls',
                             f'{datestring}_{filename}')
        m.figures_path = join(Measurement.figures_path,
                              f'{instrument}_{id}')
        osu.create_path(m.figures_path)
        osu.create_path(m.phidls_path)
        m.name = f'{m.instrument}{m.id}_{int(scriptrow)}'
        return m

    def get_filename(self, seq_number):
        return self.filename.replace('#', str(seq_number))

    def get_fitpath(self, model):
        out = os.path.join(self.path,
                           f'{model.name}_{self.instrument}_{self.id}')
        return out

    def get_fitfile(self, model, seq_number):
        out = os.path.join(self.get_fitpath(model),
                           self.get_filename(seq_number))
        return out

    def get_scriptrow(self, i):
        return self.Tdict[i]

    def get_metadata(self, seq, name):
        with open(self.rawdata_filename(seq)) as f:
            for line in f:
                if line.split(':')[0] == name:
                    return line.split(':')[1]

    def get_visc(self):
        visc = self.get_metadata(self.seq_numbers[0], 'Viscosity (mPas)')
        visc = float(visc)/1000
        return visc

    def get_n(self):
        n = self.get_metadata(self.seq_numbers[0], 'Refractive index')
        n = float(n)
        return n

    def get_TK(self):
        return self.get_TC() + 273.15

    def get_TC(self):
        return self.TC

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
        return dbo.get_table(f'{self.instrument}_exceptions')

    def get_badangles(self):
        query = f"""
        SELECT angle FROM {self.instrument}_badangles
        WHERE {self.instrument}_id = {self.id}
        """
        out = dbo.execute(query)
        return out

    def add_badangle(self, angle, reason=''):
        query = f"""
        INSERT INTO {self.instrument}_badangles
        VALUES ({self.id}, {angle}, '{reason}');
        """
        self.badangles.append(angle)
        out = dbo.execute(query)
        return out

    def add_badangles(self, badangles, reason=''):
        for phi in badangles:
            self.add_badangle(phi, reason=reason)
        return dbo.get_table(f'{self.instrument}_badangles')

    def get_data(self, seq_number, **kwargs):
        if 'xmin' in kwargs:
            xmin = kwargs.pop('xmin')
        else:
            xmin = 7
        if 'xmax' in kwargs:
            xmax = kwargs.pop('xmax')
        else:
            xmax = 170 if self.mode == 'mod3d' else 221
        filename = self.rawdata_filename(seq_number)
        df = pd.read_csv(filename, sep='\s+',
                         skiprows=16+xmin, nrows=xmax-xmin,
                         header=None, names=['t', 'g2'], engine='python')
        return df

    def get_summary(self):
        """
        gets the summary file as a pandas dataframe
        """
        filename = self.rawdata_filename('')
        # Scattering angle   Mean_CR * sin(angle) / Laser Intensity (kHz/mW)
        # g2(t=0)-1    CR CHA (kHz)    CR CHB (kHz)    Temperature (K)
        # Laser intensity (mW) Hydrodynamic Radius (nm) Width (nm)
        df = pd.read_csv(filename, sep='\t', skiprows=1, header=None,
                         names=['angle', 'I', 'g20', 'CRA', 'CRB',
                                'TK', 'I0', 'rh', 'width'])
        # Sequence number is not in the file so we find it from the index:
        if os.path.exists(self.rawdata_filename('0')):
            df['seq_number'] = range(len(df))
        else:
            df['seq_number'] = range(1, len(df)+1)
        # select only rows that correspond to self
        mask = df.seq_number.apply(
                lambda x: any(item for item in self.seq_numbers if item == x))
        df = df[mask]
        return df

    def get_sls(self, buf=None, toluene=None):
        df_summary = self.get_summary()
        dict_sls = {'angle': self.angles,
                    'q': [self.q(phi) for phi in self.angles],
                    # 'g20' : [],
                    'CRA': [],
                    'CRB': [],
                    'I0': [],
                    'Isample': [],
                    # 'err_g20' : [],
                    'err_CRA': [],
                    'err_CRB': [],
                    'err_I0': [],
                    'err_Isample': []}
        for angle in self.angles:
            # g20s = []
            CRAs = []
            CRBs = []
            I0s = []
            Isamples = []
            for seq in self.phirange(angle):
                # g20s.append(float(
                #   df_summary.loc[df_summary.seq_number == seq].g20))
                # print(df_summary.loc[df_summary.seq_number == seq].CRA)
                CRAs.append(
                    df_summary.loc[df_summary.seq_number == seq].CRA)
                CRBs.append(
                    df_summary.loc[df_summary.seq_number == seq].CRB)
                I0s.append(
                    df_summary.loc[df_summary.seq_number == seq].I0)
                Isamples.append(
                    df_summary.loc[df_summary.seq_number == seq].I)
            # dict_sls['g20'].append(np.mean(g20s))
            dict_sls['CRA'].append(np.mean(CRAs))
            dict_sls['CRB'].append(np.mean(CRBs))
            dict_sls['I0'].append(np.mean(I0s))
            dict_sls['Isample'].append(np.mean(Isamples))
            # dict_sls['err_g20'].append(np.std(g20s))
            dict_sls['err_CRA'].append(np.std(CRAs))
            dict_sls['err_CRB'].append(np.std(CRBs))
            dict_sls['err_I0'].append(np.std(I0s))
            dict_sls['err_Isample'].append(np.std(Isamples))  # beta correction
        df = pd.DataFrame(dict_sls)
        if buf is None:
            return df
        buf.angles = self.angles
        df['Ibuf'] = buf.get_sls().Isample
        df['err_Ibuf'] = buf.get_sls().err_Isample
        df['I'] = df.Isample - df.Ibuf
        return df

    def get_average_g2(self, phi):
        seq_numbers = self.phirange(phi)
        dfout = self.get_data(seq_numbers[0])
        dfout.g2 = np.zeros_like(dfout.g2)
        dfout['err_g2'] = np.zeros_like(dfout.g2)
        dfs = []
        n = 0
        for seq_number in seq_numbers:
            df = self.get_data(seq_number)
            dfout.g2 += df.g2
            n += 1
            dfs.append(df)
        dfout.g2 = dfout.g2/n
        n = 0
        for seq_number in seq_numbers:
            df = self.get_data(seq_number)
            dfout.err_g2 += (df.g2 - dfout.g2) ** 2
            n += 1
        dfout.err_g2 = np.sqrt(dfout.err_g2 / n)
        dfout['rapp'] = self.Rfromt(dfout.t, self.qq(phi))
        return dfout, dfs

    def get_phidls_filename(self, angle, end='asc'):
        filename = join(self.phidls_path, f'T{self.get_TC()}_phi{angle}.{end}')
        return filename

    def write_phidls_file(self, phi):
        filename = self.get_phidls_filename(phi)
        dfg2, dfs = self.get_average_g2(phi)
        dfg2.t = dfg2.t
        ru.write_ALV6000(self, dfg2, dfI=None, filename=filename)

    def write_phidls_files(self):
        for phi in self.angles:
            self.write_phidls_file(phi)

    def get_res(self, sorphi, source='phidls'):
        """
        sorphi: seq number or angle, depending on source
        """
        if source == 'phidls':
            filename = self.get_phidls_filename(sorphi, end='res')
        elif source == 'joALV':
            filename = self.get_joALV_filename(sorphi, end='res')
        df = pd.read_csv(filename, sep=',', header=None,
                         names=['t', 'res', 'g2', 'fit'])
        df.t = 10 ** df.t
        if source == 'joALV':
            df.t = df.t / 1000
        return df

    def get_moA(self, sorphi, A='A', source='phidls'):
        """
        sorphi: seq number or angle, depending on source
        """
        if source == 'phidls':
            filename = self.get_phidls_filename(sorphi, end=f'mo{A}')
        elif source == 'joALV':
            filename = self.get_joALV_filename(sorphi, end=f'mo{A}')
        outdict = {'peaks': [],
                   'subpeaks': []}
        with open(filename) as f:
            for line in f:
                whatpeaks = 'peaks'
                if line[0] == '|':
                    line = line[1::]
                    whatpeaks = 'subpeaks'
                try:
                    amount, value = line.split('   ')
                    peak = [float(amount), float(value)]
                    outdict[whatpeaks].append(peak)
                except:
                    # print(line[0:6])
                    # print(line.split('   '))
                    if line[0:6] == 'SqBeta':
                        outdict['SqBeta'] = float(line.split('=')[1])
                    if line[0:8] == 'Baseline':
                        outdict['Baseline'] = float(line.split('=')[1])
        return outdict

    def get_Arl(self, sorphi, A='A', rmin=0, rmax=np.inf, source='phidls'):
        """
        sorphi: seq number or angle, depending on source
        """
        if source == 'phidls':
            filename = self.get_phidls_filename(sorphi, end=f'{A}rl')
            qq = self.qq(sorphi)
        elif source == 'joALV':
            filename = self.get_joALV_filename(sorphi, end=f'{A}rl')
            qq = self.qq(self.phifromseq(sorphi))
        df = pd.read_csv(filename, sep=',', header=None,
                         names=['t', 'dist'])
        df['logt'] = df.t
        df.t = 10 ** df.t
        df['rapp'] = self.Rfromt(df.t, qq)
        df = df.loc[df.rapp > rmin]
        df = df.loc[df.rapp < rmax]
        if source == 'joALV':
            df.t = df.t / 1000
        return df

    def get_joALV_path(self):
        return join(self.path, 'joALV')

    def get_joALV_filename(self, seq_number, end='res'):
        filename = self.filename.split('#')[0]
        filename = filename + '{0:03}'.format(seq_number)
        filename = filename + '_joALV.' + end
        return join(self.get_joALV_path(), filename)


    def get_phidlstable(self, fit):
        return fit.get_phidlstable(self)

    def get_fitpar(self, fit, parameter):
        pass

    def get_gammastar(self):
        pass

    def phirange(self, phi):
        "returns all the seq_numbers of phi. Exceptions excepted."
        index = np.where(self.allangles == phi)[0]
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
        for phi in self.allangles:
            if self.seqinphi(seq_number, phi):
                return phi

    def RfromD(self, D):
        TK = self.get_TK()
        # visc = self.sample.get_viscosity(TK)
        visc = self.get_visc()
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
        visc = self.get_visc()
        out = constants.Boltzmann*TK / \
            (6*np.pi*visc*rh)
        return out

    def GfromR(self, R, qq):
        return self.DfromR(R) * qq


    def tfromR(self, rh, qq):
        D = self.DfromR(rh)
        return 1 / (D * qq)

    def q(self, phi):
        """
        calculates the scattering vector q[m^-1]
        from the scattering angle 2\\Theta.
        """
        wl = self.get_wl()
        n = self.get_n()
        out = n * 4 * np.pi * np.sin((phi*np.pi/360))/wl
        return out

    def qq(self, phi):
        """
        calculates the square of scattering vector q^2[m^-2]
        from the scattering angle 2\\Theta.
        """
        return self.q(phi)**2

    def plot_data(self, seq_number, ax=None, **kwargs):
        """
        This should most likely not be a function of the lsi class.
        """
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_xlabel('$\\tau$ [s]')
            ax.set_ylabel('$g_2-1$')
        if 'marker' not in kwargs:
            kwargs['marker'] = '.'
        if 'linestyle' not in kwargs:
            kwargs['linestyle'] = ''
        df = self.get_data(seq_number)

        # Autocorrelation Plots:
        ax.plot(df.t, df.g2, **kwargs)
        ax.set_xscale('log')
        ax.set_ylim(0, 1)
        if self.mode == '3dcross':
            ax.set_ylim(0, 0.25)
        if self.mode == 'mod3d':
            ax.set_ylim(0, 0.85)
        return ax

    def plot_fit(self, seq_number, fit, fitcolor='black', ax=None,
                 showres=True, **kwargs):
        ax = self.plot_data(seq_number, ax=ax, **kwargs)
        if fit is None:
            return ax
        fitdf = fit.get_fit(self, seq_number)
        fitdf = fitdf[fitdf.fit < 1]
        ax.plot(fitdf.t, fitdf.fit, '-', color=fitcolor)
        rescolor = fitcolor
        if 'color' in kwargs:
            rescolor = kwargs['color']
        if showres:
            axres = ax.twinx()
            axres.plot(fitdf.t, fitdf.res, alpha=0.2, color=rescolor)
            axres.set_ylabel('Relative Residuals')
            axres.set_ylim(-0.05, 0.05)
        return ax

    def plot_dist(self, seq_number, fit, ax=None, xspace='t', **kwargs):
        "plot distribution, if exists, in equal area representation."
        if ax is None:
            fig, ax = plt.subplots()
        dfd = fit.get_dist(self, seq_number)
        if xspace == 'rh':
            # rapp in nm as xaxis
            ax.plot(dfd.rh * 10**9, dfd.dist, **kwargs)
        elif xspace == 't':
            ax.plot(dfd.t, dfd.dist, **kwargs)
        ax.set_xscale('log')
        return ax

    def qplot(self, fit, parameters, function=lambda x: x,
              err_function=lambda x, err_x: err_x, ax=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_ylabel(fit.pardict[parameters[0]][0])
            ax.set_xlabel("$q^2$ $\\mathrm{[\\mu m^{-2}]}$")
        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'
        if 'linestyle' not in kwargs:
            kwargs['linestyle'] = ''

        df = fit.get_phitable(self)
        xdata = df.qq
        args = [df[par] for par in parameters]
        ydata = function(*args)
        err_args = [df[f'err_{par}'] for par in parameters]
        err_ydata = err_function(*args, *err_args)
        ax.errorbar(xdata, ydata, err_ydata, **kwargs)
        ax.set_xlim(left=0)
        return ax
