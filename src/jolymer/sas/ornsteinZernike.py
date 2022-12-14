# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 17:00:15 2020

@author: xcill
"""

from . import SAXS_Model as sasmodel
import numpy as np

def ornsteinZernike(q, xi, oz_exp, oz_scale):
    q = np.array(q)
    I = oz_scale / (1+ (q*xi)**oz_exp) 
    return I

class OZ(sasmodel.SAXS_Model):
    
    def __init__(self):
        self.name = 'oz'
        self.longname = 'Ornstein Zernike'
        self.parameters = ['xi', 'oz_exp', 
                           'oz_scale']
        self.pdict = {'xi' : ['$\\xi$', 'nm']}
        self.fitfunc = ornsteinZernike

    def get_text(self, fit_dict):
        text = """
        $\\xi = {:.3f} \\pm {:.3f}$ nm
        $n = {:.3f} \\pm {:.3f}$
        $A = {:.3f}
        """.format(fit_dict['xi'], fit_dict['std_xi'], 
                   fit_dict['oz_exp'], fit_dict['std_oz_exp'], 
                   fit_dict['oz_scale'])
        return text
        
    def get_clustering_strength(self, fitdict, q):
        A = fitdict['porod_scale']
        n = fitdict['porod_exp']
        out = A/q**n
        return out
        
 
background = sasmodel.Background()
oz = OZ()
oz_bg = oz.plusmodel(background)
def get_text_oz(fit_dict):
    text="""
    $n =$ {0:.2f} $\\pm$ {1:.2f}
    $\\xi =$ {2:.2f} $\\pm$ {3:.2f} nm
    $C = $ {4:.2E} 
    $\chi^2=$ {5:.2f}
    """.format(fit_dict['oz_exp'], fit_dict['std_oz_exp'], 
                  fit_dict['xi'], fit_dict['std_xi'], 
                  fit_dict['oz_scale'], 
                  fit_dict['chi2'])   
    return text
oz_bg.get_text = get_text_oz

forward = sasmodel.Porod()

oz_porod = oz.plusmodel(forward, name='oz_forward', longname = 'Ornstein Zernike')       

def get_text_ozp(fit_dict):
    text="""
    $n =$ {0:.2f} $\\pm$ {1:.2f}
    $m_F =$ {2:.2f} $\\pm$ {3:.2f}
    $\\xi =$ {4:.2f} $\\pm$ {5:.2f} nm
    $C = $ {6:.2E} 
    $A_F = $ {7:.2E}
    $\chi^2=$ {8:.2f}
    """.format(fit_dict['oz_exp'], fit_dict['std_oz_exp'], 
                  fit_dict['porod_exp'],fit_dict['std_porod_exp'],
                  fit_dict['xi'], fit_dict['std_xi'], 
                  fit_dict['oz_scale'], fit_dict['porod_scale'],
                  fit_dict['chi2'])
    return text
oz_porod.get_text = get_text_ozp
