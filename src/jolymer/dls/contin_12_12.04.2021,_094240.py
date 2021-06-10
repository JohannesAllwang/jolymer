#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 08:10:31 2020

@author: johannes
"""


# import jscatter as js
from . import DLS_Measurement
import numpy as np
import pandas as pd
import os
from .cun import Cun
from .. import database_operations as dbo

from . import CONTINwrapper

def create_path(path):
    "Tryes to create path and prints success"
    try:
        os.mkdir(path)
    except OSError:
        print ("Creation of the directory %s failed" % path)
    else:
        print ("Successfully created the directory %s " % path)

class Contin(Cun):
    
    def __init__(self, name, parameters, labels):
        self.name = name
        self.parameters = parameters
        self.labels = labels
        
    def get_fitfile(self, measurement, seq_number):
        path = self.get_continpath(measurement)
        return f'{path}contin_fit{seq_number}.csv'
    def get_distfile(self, measurement, seq_number):
        path = self.get_continpath(measurement)
        return f'{path}contin_dist{seq_number}.csv'
    def fit_measurement(self, measurement, seq_number):
        pars = []
        fit, distribution, alpha = CONTINwrapper.continfit(measurement, seq_number)
        distt, distdist = [np.array(x) for x in distribution[0:2]]
        distdist = distdist/np.sum(distdist)
        dist = pd.DataFrame({'t':distt, 'dist':distdist})
        mean_of_logdist = np.sum(np.log(dist.t) * dist.dist)
        # plt.plot(dist.t, dist.dist)
        # plt.xscale('log')
        # plt.show()
        print(seq_number, np.sum(dist.dist))
        geomean = np.exp(mean_of_logdist)
        armean = np.sum(dist.t * dist.dist)
        harmean = 1 / np.sum(dist.dist/dist.t)
        pars.append(alpha)
        pars.append(geomean)
        pars.append(armean) # 1/ armean(tau) = harmeanD
        pars.append(harmean) # 1 / harmean(tau) = armeanD
        fit.to_csv(self.get_fitfile(measurement, seq_number))
        dist.to_csv(self.get_distfile(measurement, seq_number))
        return pars
    
    def df_fit(self, measurement, seq_number):
        df = pd.read_csv(self.get_fitfile(measurement, seq_number), index_col=(0))
        return df
    def df_dist(self, measurement, seq_number):
        df = pd.read_csv(self.get_distfile(measurement, seq_number), index_col=(0))
        return df
    
    def analyse_phi(self, measurement, df_avg, df_par, phi):
        par_list=[]
        qq = measurement.qq(phi)
        for seq_number in measurement.phirange(phi):
            if seq_number in measurement.exceptions:
                continue
            pars = self.fit_measurement(measurement, seq_number)
            par_list.append([1/(pars[0] * qq), 1/(pars[1] * qq), 1/( pars[2]*qq), pars[3]])
            dict_par = {'seq_number' : seq_number}
            for index, parameter in enumerate(self.parameters):
                dict_par[parameter] = pars[index]
            df_par = df_par.append(dict_par, ignore_index=True)
        dict_avg ={'qq' : measurement.qq(phi)}
        for index, parameter in enumerate(self.parameters):
            dict_avg[parameter] = np.mean([x[index] for x in par_list])
            dict_avg['err_' + parameter] = np.std([x[index] for x in par_list])
        df_avg = df_avg.append(dict_avg, ignore_index=True)
        return df_avg, df_par
    
    def phis_tablename(self, measurement):
        out = f"phis_{self.name}_{measurement.id}"
        return out
    def seqs_tablename(self, measurement):
        out = f"seqs_{self.name}_{measurement.id}"
        return out
    def analyse(self, measurement):
        create_path(self.get_continpath(measurement))
        df_avg = pd.DataFrame(columns = ["qq", "Dgeo", "Dar", 'Dhar'])
        df_par = pd.DataFrame(columns = ["seq_number", "D", "std_D"])
        for phi in measurement.angles:
            df_avg, df_par = self.analyse_phi(measurement, df_avg, df_par, phi)
        df_avg.to_sql(self.phis_tablename(measurement), dbo.conn, if_exists='replace')
        df_par.to_sql(self.seqs_tablename(measurement), dbo.conn, if_exists='replace')
        print(df_avg)
    
    
contin = Contin('contin', None, None)
