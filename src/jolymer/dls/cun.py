#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 12:35:09 2020

@author: johannes
"""

import pandas as pd
import numpy as np
import datetime as dt
from scipy import optimize, constants

from .. import database_operations as dbo


class Cun:

    def __init__(self, name, pardict, create_fitfunc, create_bounds):
        self.name = name
        self.parameters = pardict.keys()
        self.pardict = pardict
        self.create_fitfunc = create_fitfunc
        self.create_bounds = create_bounds
        self.seq_columns = ['seq_number']
        self.phi_columns = ['phi', 'qq']
        for par in self.parameters:
            self.seq_columns.append(par)
            self.seq_columns.append(f'std_{par}')
            self.phi_columns.append(par)
            self.phi_columns.append(f'err_{par}')

    @staticmethod
    def fit_data(df, fitfunc, bounds):
        popt, pcov = optimize.curve_fit(fitfunc, df.t, df.g2, bounds=bounds)
        return popt, pcov

    def fit_m(self, m, seq_number, fitfunc, bounds):
        df = m.get_data(seq_number)
        popt, pcov = self.fit_data(df, fitfunc, bounds)
        return popt, pcov

    def fit_phidls(self, m, phi, fitfunc, bounds):
        df, dfs = m.get_average_g2(phi)
        popt, pcov = optimize.curve_fit(fitfunc, df.t, df.g2, sigma=df.err_g2,
                                        bounds=bounds)
        return popt, pcov

    def get_seqtable(self, m):
        fitpars = dbo.get_table(self.seqs_tablename(m))
        fitpars = fitpars.set_index('seq_number', drop=True)
        return fitpars

    def get_phitable(self, m):
        fitpars = dbo.get_table(self.phis_tablename(m))
        return fitpars

    def get_phidlstable(self, m):
        fitpars = dbo.get_table(self.phidls_tablename(m))
        return fitpars

    def get_fit(self, m, seq_number):
        data = m.get_data(seq_number)
        with dbo.dbopen() as c:
            query = f"""SELECT * FROM {self.seqs_tablename(m)}
                    WHERE seq_number={seq_number}"""
            conn = c.connection
            fitpars = pd.read_sql_query(query, conn)

        fitpars = fitpars.set_index('seq_number', drop=True)
        popt = []
        for parameter in self.parameters:
            popt.append(fitpars.loc[seq_number, parameter])
        qq = m.qq(m.phifromseq(seq_number))
        fitfunc = self.create_fitfunc(m, qq)
        data['fit'] = fitfunc(data.t, *popt)
        data['res'] = data.g2 - data.fit
        return data

    def get_phidlsfit(self, m, phi):
        data, dfs = m.get_average_g2(phi)
        with dbo.dbopen() as c:
            query = f"""SELECT * FROM {self.phidls_tablename(m)}
                    WHERE phi={phi}"""
            conn = c.connection
            fitpars = pd.read_sql_query(query, conn)
        fitpars = fitpars.set_index('phi', drop=True)
        popt = []
        for parameter in self.parameters:
            popt.append(fitpars.loc[phi, parameter])
        qq = m.qq(phi)
        fitfunc = self.create_fitfunc(m, qq)
        data['fit'] = fitfunc(data.t, *popt)
        data['res'] = data.g2 - data.fit
        # print(data)
        return data

    def analyse_phi(self, m, df_avg, df_par, phi):
        par_list = []
        qq = m.qq(phi)
        fitfunc = self.create_fitfunc(m, qq)
        bounds = self.create_bounds(m, phi)
        for seq_number in m.phirange(phi):
            if seq_number in m.exceptions:
                continue
            popt, pcov = self.fit_phidls(m,
                                         seq_number, fitfunc, bounds)
            std_popt = np.sqrt(np.diag(pcov))
            par_list.append(popt)
            dict_par = {'seq_number': seq_number}
            for index, parameter in enumerate(self.parameters):
                dict_par[parameter] = popt[index]
                dict_par['std_' + parameter] = std_popt[index]
            df_par = df_par.append(dict_par, ignore_index=True)
        dict_avg = {'phi': phi,
                    'qq': m.qq(phi)}
        for index, parameter in enumerate(self.parameters):
            dict_avg[parameter] = np.mean([x[index] for x in par_list])
            dict_avg['err_' + parameter] = np.std([x[index] for x in par_list])
        df_avg = df_avg.append(dict_avg, ignore_index=True)
        return df_avg, df_par

    def phis_tablename(self, m):
        out = f"phis_{self.name}_{m.id}"
        return out

    def seqs_tablename(self, m):
        out = f"seqs_{self.name}_{m.id}"
        return out

    def phidls_tablename(self, m):
        out = f'phidls_{self.name}_{m.id}'
        return out

    def get_R(self, m):
        df = self.get_phitable(m)
        D = df.loc[df.phi == 'y_intercept'].Dapp
        err_D = df.loc[df.phi == 'y_intercept'].err_Dapp
        R = float(m.RfromD(D))
        err_R = float(R * err_D / D)
        return R, err_R

    def str_R(self, m):
        R, err_R = self.get_R(m)
        # return '{0:.2f} ± {1:.2f} nm'.format(R * 1e9, err_R * 1e9)
        return '{0:.2f} $\\pm$ {1:.2f} nm'.format(R * 1e9, err_R * 1e9)

    def analyse(self, m):
        df_avg = pd.DataFrame(columns=self.phi_columns)
        df_par = pd.DataFrame(columns=self.seq_columns)
        for phi in m.angles:
            print(f'Analyzing angle {phi}')
            df_avg, df_par = self.analyse_phi(m, df_avg, df_par, phi)
        # Now linar and constant fits on the qdata:

        dict_avg = {'phi': ['constant', 'y_intercept', 'slope'],
                    'qq': [-1, -1, -1]}
        for index, par in enumerate(self.parameters):
            # Constant fit:
            def func(qq, Davg):
                return Davg * np.ones(len(qq))
            print(par)
            popt, pcov = [[None], [None]]
            try:
                popt, pcov = optimize.curve_fit(func, df_avg.qq, df_avg[par],
                                                sigma=df_avg[f'err_{par}'],
                                                bounds=(-np.inf, np.inf))
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
                                                bounds=((-np.inf, -np.inf),
                                                        (np.inf, np.inf)))
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
            df_avg.to_sql(self.phis_tablename(m), conn,
                          if_exists='replace', index=False)
            df_par.to_sql(self.seqs_tablename(m), conn,
                          if_exists='replace', index=False)

    def analyse_phidls(self, m):
        par_list = []
        dflist = []
        df_par = pd.DataFrame(columns=self.phi_columns)
        for phi in m.angles:
            qq = m.qq(phi)
            fitfunc = self.create_fitfunc(m, qq)
            bounds = self.create_bounds(m, phi)
            popt, pcov = self.fit_phidls(m, phi, fitfunc, bounds)
            std_popt = np.sqrt(np.diag(pcov))
            par_list.append(popt)
            dict_par = {'phi': phi,
                        'qq': qq}
            for index, parameter in enumerate(self.parameters):
                dict_par[parameter] = popt[index]
                dict_par['std_' + parameter] = std_popt[index]
            df_par = df_par.append(dict_par, ignore_index=True)
        with dbo.dbopen() as c:
            conn = c.connection
            df_par.to_sql(self.phidls_tablename(m), conn,
                          if_exists='replace', index=False)
        return df_par, dflist


def c2_create_fitfunc(m, qq):
    def inner(x, Dapp, Dapp2, beta):
        g1 = np.exp(-Dapp*x * qq) * \
            (1 + Dapp2 * qq**2 * x**2 / 2)
        g2 = beta * g1**2
        return g2
    return inner


def cu2_create_bounds(m, phi):
    Dmin = m.DfromR(m.rmax)
    Dmax = m.DfromR(m.rmin)
    # Dmin = 0
    # Dmax = 10000
    D2min = 0
    D2max = Dmax**2 / 10000000
    # D2max = 1000_000_000_000
    betamin = 0
    betamax = 1.2
    if m.mode == '3dcross':
        betamax = 0.3
    if m.mode == 'mod3d':
        betamax = 0.9
    return ((Dmin, D2min, betamin), (Dmax, D2max, betamax))


cu2_pardict = {
    'Dapp': ["$D_{app}$ $\\mathrm{[m^2/s}]$"],
    'varDapp': ['var$(D_{app})$'],
    'beta': ['$\\beta$']
    }

cu2 = Cun('cu2', cu2_pardict, c2_create_fitfunc, cu2_create_bounds)
