# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 14:12:42 2020

@author: xcill
"""

from scipy import optimize
import numpy as np
from .. import database_operations as dbo
from ..Measurement import Measurement

# import database_operations as dbo

class SAXS_Measurement(Measurement):
    
    def __init__(self, id):
        
        self.id = id
        
    
    def get_data(self):
        pass
    
    @staticmethod
    def get_distribution(df):
        
        pass
    
    
    
def gen_guinier_fitfunc(alpha):
    def inner(q, Rg, A):
        if alpha == 0:
            pre = 1
        elif alpha == 1 or alpha == 2:
            pre = alpha * np.pi * q ** -alpha
        else:
            raise TypeError('alpha needs to be in 0,1,2')
            
        I = pre * A * np.exp(-Rg ** 2 * q ** 2 / (3 - alpha))
        return I
    return inner

def guinier_porod_3D(q, Rg1, s1, Rg2, s2, G2, dd):
    q = np.atleast_1d(q)

    # define parameters for smooth transitions
    Q1 = (1 / Rg1) * ((dd - s1) * (3 - s1) / 2) ** 0.5
    Q2 = ((s1 - s2) / (2 / (3 - s2) * Rg2 ** 2 - 2 / (3 - s1) * Rg1 ** 2)) ** 0.5
    G1 = G2 / (np.exp(-Q2 ** 2 * (Rg1 ** 2 / (3 - s1) - 
                                  Rg2 ** 2 / (3 - s2))) * Q2 ** (s2 - s1))
    D = G1 * np.exp(-Q1 ** 2 * Rg1 ** 2 / (3 - s1)) * Q1 ** (dd - s1)

    # define functions in different regions
    def _I1_3regions(q):
        res = G2 / q ** s2 * np.exp(-q ** 2 * Rg2 ** 2 / (3 - s2))
        return res

    def _I2_3regions(q):
        res = G1 / q ** s1 * np.exp(-q ** 2 * Rg1 ** 2 / (3 - s1))
        return res

    def _I3_3regions(q):
        res = D / q ** dd
        return res

    I = np.piecewise(q, [q < Q2, (Q2 <= q) & (q < Q1), q >= Q1],
                     [_I1_3regions, _I2_3regions, _I3_3regions])
    return I

