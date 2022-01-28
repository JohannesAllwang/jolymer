#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 08:10:31 2020

@author: johannes
"""


# import jscatter as js
from .cun import Cun
from . import lsi
from .. import database_operations as dbo
from .. import os_utility as osu
from . import CONTINwrapper

import numpy as np
import pandas as pd
import os
from scipy import optimize

class _peak:

    def __init__(self, dfd, fromto, m, seq_number):
        self.fromto = fromto
        self.frm, self.to = fromto
        self.dfd = dfd[dfd.t<=self.frm]
        self.dfd = dfd[dfd.t>=self.to]
        self.phi = m.phifromseq(seq_number)
        self.qq = m.qq(self.phi)
        self.Rfromto = [m.Rfromt(ft, self.qq) for ft in fromto]
        self.Ffrom, self.Rto = self.Rfromto
        self.dfd['R'] = m.Rfromt(self.dfd.t, self.qq)

    def get_m(self, n):
        if n==0:
            return np.sum(dfd.dist)
        dfd = self.dfd
        out = np.sum(dfd.dist * dfd.t**n)/self.get_m(0)
        return out

    def get_D(self, n):
        return 1 / self.get_m(n) / self.qq

    def Rin(self, R):
        if R>self.Rfrom and R<self.Rto:
            return True
        return False

class Contin(Cun):

    def __init__(self, name, parameters, labels):
        self.name = name
        self.parameters = parameters
        self.labels = labels

    def get_fitfile(self, m, seq_number):
        out = os.path.join(m.get_fitpath(self),
                f'fit_{seq_number}.contin')
        return out

    def get_distfile(self, m, seq_number):
        out = os.path.join(m.get_fitpath(self),
                f'dist_{seq_number}.contin')
        return out

    def get_rawoutputfile(self, m, seq_number):
        out = os.path.join(m.get_fitpath(self),
                f'full_output_{seq_number}.contin')
        return out

    def run_contin(self, m, seq_number):
        rawoutputfile = self.get_rawoutputfile(m, seq_number)
        CONTINwrapper.continfit(m, seq_number, rawoutputfile=rawoutputfile)

    def evaluate_continfile(self, m, seq_number):
        df0 = m.get_data(seq_number)
        lenX = len(m.get_data(seq_number))
        ngrid = 256
        with open(contin.get_rawoutputfile(m, seq_number)) as f:
            i=0
            p = False
            lines = list(f)
            for line in lines:
                i+=1
                if len(line.split('++  CHOSEN'))>1:
                    p=True
                    ichosen = i
                if len(line.split('PEAK')) > 1 and p:
                    ipeak = i
                    p=False
            n0 = 0
            ordinate_fit = []
            abscissa_fit = []
            for line in lines[ichosen - 2*lenX -2: ichosen-1]:
                n0 +=1
                line = line.replace('O', ' ')
                line = line.replace('X', ' ')
                line = line.replace('*', ' ')
                if 'ABSCISSA' in line.split():
                    print('ABSCISSa', n0, lenX)
                if not line.isspace():
                    ordinate_fit.append(float(line.split()[0].replace('D', 'E')))
                    abscissa_fit.append(float(line.split()[1].replace('D', 'E')))
            ordinate_dist = []
            error_dist = []
            abscissa_dist = []
            for line in lines[ichosen+14:ichosen+14+ ngrid*2]:
                line = line.replace('O', ' ')
                line = line.replace('..', ' ')
                line = line.replace('.X', ' ')
                line = line.replace('X', ' ')
                line = line.replace('*', ' ')
                if not line.isspace():
                    if line[0]=='.':
                        ordinate_dist.append(0.0)
                        error_dist.append(float(line.split()[1].replace('D', 'E')))
                        abscissa_dist.append(float(line.split()[2].replace('D', 'E')))
                    else:
                        ordinate_dist.append(float(line.split()[0].replace('D', 'E')))
                        error_dist.append(float(line.split()[1].replace('D', 'E')))
                        abscissa_dist.append(float(line.split()[2].replace('D', 'E')))
            fit_quality = lines[ichosen+9].split()

            peak_fromtos = []
            for line in lines[ichosen+14+ngrid*2::]:
                if line[0:5]=='0PEAK':
                    line = line.split('FROM')[1].split('J')[0]
                    fromto = [float(x) for x in line.split('TO')]
                    peak_fromtos.append(fromto)
        df_fit = pd.DataFrame(
                {'t' : df0.t,
                'abscissa': np.array(abscissa_fit),
                'g2' : df0.g2,
                'fit': np.array(ordinate_fit)**2})
        df_fit['res'] = df_fit.g2 - df_fit.fit
        df_dist = pd.DataFrame(
                {'t' : np.array(abscissa_dist),
                'dist': np.array(ordinate_dist),
                'err_dist' : np.array(error_dist)})
        df_fit.to_csv(self.get_fitfile(m, seq_number))
        df_dist.to_csv(self.get_distfile(m, seq_number))
        return df_fit, df_dist, peak_fromtos, fit_quality

    def get_fit(self, m, seq_number):
        df = pd.read_csv(self.get_fitfile(m, seq_number), index_col=(0))
        return df

    def get_dist(self, m, seq_number):
        df = pd.read_csv(self.get_distfile(m, seq_number), index_col=(0))
        return df

    def analyse_phi(self, m, df_avg, df_par, peaks_dict, phi):
        par_list=[]
        qq = m.qq(phi)
        for seq_number in m.phirange(phi):
            if seq_number in m.exceptions:
                continue
            print(m.get_filename(seq_number))
            self.run_contin(m, seq_number)
            dff, dfd, peak_fromtos, fit_quality = self.evaluate_continfile(m, seq_number)
            # Correctly weighting the distribution:
            dfd.dist = dfd.dist * dfd.t / np.sum(dfd.dist * dfd.t)
            dfd['rh'] = m.Rfromt(dfd.t, qq)
            dff.to_csv(self.get_fitfile(m, seq_number))
            dfd.to_csv(self.get_distfile(m, seq_number))
            # Get the parameters about the fit firts:
            alpha, alpha_by_s, obj_fctn, variance, std_dev, deg_freedon, probreg1, probreg2= \
                    fit_quality
            # Now we need to calculate some parameters:
            mean_tau = np.sum(dfd.dist * dfd.t)
            Dapp = 1 / (qq * mean_tau)
            # Create a parameter dict to append to the df_par dataframe:
            dict_par = {'seq_number' : seq_number}
            for value, parameter in zip([qq, Dapp] + fit_quality, self.parameters):
                dict_par[parameter] = value
            df_par = df_par.append(dict_par, ignore_index=True)
            # Add all the peaks to the peaks dictionary:
            for frmto in peak_fromtos:
                peaks_dict['seq_number'].append(seq_number)
                peaks_dict['from'].append(frmto[0])
                peaks_dict['to'].append(frmto[1])
            # Append parameters to parlist for later averaging:
            par_list.append([Dapp] + [float(fq) for fq in fit_quality])
        dict_avg ={'phi' : phi,
                   'qq' : qq}
        for index, parameter in enumerate(self.parameters[1::]):
            dict_avg[parameter] = np.mean([x[index] for x in par_list])
            dict_avg['err_' + parameter] = np.std([x[index] for x in par_list])
        df_avg = df_avg.append(dict_avg, ignore_index=True)
        return df_avg, df_par, peaks_dict

    def phis_tablename(self, m):
        out = f"phis_{self.name}_{m.name}"
        return out

    def seqs_tablename(self, m):
        out = f"seqs_{self.name}_{m.name}"
        return out

    def peaks_tablename(self, m):
        out = f"peaks_{self.name}_{m.name}"
        return out

    def analyse(self, m):
        osu.create_path(m.get_fitpath(self))
        df_avg = pd.DataFrame(columns = parameters)
        df_par = pd.DataFrame(columns = parameters)
        peaks_dict = {'seq_number' : [],
                      'from' : [],
                      'to' : []}
        for phi in m.angles:
            df_avg, df_par, peaks_dict = self.analyse_phi(m, df_avg, df_par, peaks_dict, phi)
        df_peaks = pd.DataFrame(peaks_dict)
        dict_avg ={'phi' : ['constant', 'y_intercept', 'slope'],
                   'qq' : [-1, -1, -1]}
        for index, par in enumerate(self.parameters[1::]):
            # Constant fit:
            def func(qq, Davg):
                return Davg * np.ones(len(qq))
            popt, pcov = [[None], [None]]
            try:
                popt, pcov = optimize.curve_fit(func, df_avg.qq, df_avg[par],
                                    sigma=df_avg[f'err_{par}'], bounds=(-np.inf ,np.inf))
                pstd = np.sqrt(np.diag(pcov))
            except Exception as e:
                print(e)
                pstd = [None]
            constant = popt[0]
            err_constant = pstd[0]

            # Linear fit:
            def func(qq, D0, D0chrg2):
                out = D0*np.ones(len(qq)) + D0chrg2*qq
                return out
            try:
                popt, pcov = optimize.curve_fit(func, df_avg.qq, df_avg[par],
                                    sigma=df_avg[f'err_{par}'],
                                    bounds=((-np.inf,-np.inf), (np.inf, np.inf)))
                pstd = np.sqrt(np.diag(pcov))
            except Exception as e:
                print(e)
                popt, pcov = [[None, None], [None, None]]
                pstd = [None, None]
            yintercept, slope = popt
            err_yintercept, err_slope = pstd

            dict_avg[par] = [constant, yintercept, slope]
            dict_avg['err_'+par] = [err_constant, err_yintercept, err_slope]
        ap = pd.DataFrame(dict_avg)
        df_avg = df_avg.append(ap, ignore_index=True)
        with dbo.dbopen() as c:
            conn = c.connection
            df_avg.to_sql(self.phis_tablename(m), conn, if_exists='replace')
            df_par.to_sql(self.seqs_tablename(m), conn, if_exists='replace')
            df_peaks.to_sql(self.peaks_tablename(m), conn, if_exists='replace')

    def analyse_phidls(self, m):
        osu.create_path(m.get_phidls_path(self))
        df_avg = pd.DataFrame(columns=parameters)
        df_par = pd.DataFrame(columns=parameters)
        peaks_dict = {'seq_number': [],
                      'from': [],
                      'to': []}
        # for phi in m.angles:
        #     df_avg, df_par, peaks_dict = self.analyse_phi(m, df_avg, df_par, peaks_dict, phi)
        # df_peaks = pd.DataFrame(peaks_dict)
        # dict_avg ={'phi' : ['constant', 'y_intercept', 'slope'],
        #            'qq' : [-1, -1, -1]}
        # for index, par in enumerate(self.parameters[1::]):
        #     # Constant fit:
        #     def func(qq, Davg):
        #         return Davg * np.ones(len(qq))
        #     popt, pcov = [[None], [None]]
        #     try:
        #         popt, pcov = optimize.curve_fit(func, df_avg.qq, df_avg[par],
        #                             sigma=df_avg[f'err_{par}'], bounds=(-np.inf ,np.inf))
        #         pstd = np.sqrt(np.diag(pcov))
        #     except Exception as e:
        #         print(e)
        #         pstd = [None]
        #     constant = popt[0]
        #     err_constant = pstd[0]

        #     # Linear fit:
        #     def func(qq, D0, D0chrg2):
        #         out = D0*np.ones(len(qq)) + D0chrg2*qq
        #         return out
        #     try:
        #         popt, pcov = optimize.curve_fit(func, df_avg.qq, df_avg[par],

    def get_peaks(self, m, seq_number):
        dfd = self.get_dist(m, seq_number)
        query = f"""
        SELECT * FROM {self.peaks_tablename(m)}
        WHERE seq_number = {seq_number}
        """
        print(query)
        peakargs = dbo.execute(query)
        for peakarg in peakargs:
            seq_number, fromto = peak
            peak = _peak(dfd, fromto, m, seq_number)

    def get_average(self, m, seq_number):
        df = self.get_fit(m, seq_number)
        out = np.avg(df.t * df.g2)
        return out

    def get_par(self, m, seq_number):
        df = self.seq_number

parameters = []
parameters = ["qq", "Dapp", 'alpha', 'alpha_by_s', 'obj_fctn', 'variance', 'std_dev', 'deg_freedom', 'probrej1', 'probrej2']
contin = Contin('contin', parameters, None)
