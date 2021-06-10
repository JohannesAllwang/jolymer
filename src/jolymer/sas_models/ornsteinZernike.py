# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 17:00:15 2020

@author: xcill
"""

from ..SAXS_Model import SAXS_Model
import numpy as np

def ornsteinZernike(q, xi, lorentz_exp, porod_exp, lorentz_scale, 
                    porod_scale, bg):
    q = np.array(q)
    I = lorentz_scale / (1+ (q*xi)**lorentz_exp) 
    I += porod_scale/q**porod_exp
    I += bg
    return I

class oz(SAXS_Model):
    
    def __init__(self):
        self.name = 'oz'
        self.longname = 'Ornstein Zernike'
        self.parameters = ['xi', 'lorentz_exp', 'porod_exp', 
                           'lorentz_scale', 'porod_scale', 'bg']
        self.pdict = {'xi' : ['$\\xi$', 'nm']}
        self.fitfunc = ornsteinZernike
        
    def get_clustering_strength(self, fitdict, q):
        A = fitdict['porod_scale']
        n = fitdict['porod_exp']
        out = A/q**n
        return out
        
        