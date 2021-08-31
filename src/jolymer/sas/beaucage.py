# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 10:16:44 2021

@author: xcill
"""

from . import SAXS_Model as sasmodel
import numpy as np
from scipy import special

def beaucage(q, Rg, G, d):
    r"""
    Beaucage introduced a model based on the polymer fractal model.

    Beaucage used the numerical integration form (Benoit, 1957) although the analytical
    integral form was available [1]_. This is an artificial connection of Guinier and Porod Regime .
    Better use the polymer fractal model [1]_ used in gaussianChain.
    For absolute scattering see introduction :ref:`formfactor (ff)`.

    Parameters
    ----------
    q : array
        Wavevector
    Rg : float
        Radius of gyration in 1/q units
    G : float
        Guinier scaling factor, transition between Guinier and Porod
    d : float
        Porod exponent for large wavevectors

    Returns
    -------
    dataArray
        Columns [q,Fq]

    Notes
    -----

    .. math:: I(q) &= G e^{-q^2 R_g^2 / 3.} + C q^{-d} \left[erf(qR_g / 6^{0.5})\right]^{3d}

                C &= \frac{G d}{R_g^d} \left[\frac{6d^2}{(2+d)(2+2d)}\right]^{d / 2.} \Gamma(d/2)

    with the Gamma function :math:`\Gamma(x)` .

    Polymer fractals:

    | d = 5/3    fully swollen chains,
    | d = 2      ideal Gaussian chains and
    | d = 3      globular e.g. collapsed chains. (volume scattering)
    | d = 4      surface scattering at a sharp interface/surface
    | d = 6-dim  rough surface area with a dimensionality dim between 2-3 (rough surface)
    | d < r      mass fractals (eg gaussian chain)

    The Beaucage model is used to analyze small-angle scattering (SAS) data from
    fractal and particulate systems. It models the Guinier and Porod regions with a
    smooth transition between them and yields a radius of gyration and a Porod
    exponent. This model is an approximate form of an earlier polymer fractal
    model that has been generalized to cover a wider scope. The practice of allowing
    both the Guinier and the Porod scale factors to vary independently during
    nonlinear least-squares fits introduces undesired artefact's in the fitting of SAS
    data to this model.

    .. [1] Analysis of the Beaucage model
            Boualem Hammouda  J. Appl. Cryst. (2010). 43, 1474–1478
            http://dx.doi.org/10.1107/S0021889810033856

    """
    q = np.atleast_1d(q)
    Rg = float(Rg)
    C = G * d / Rg ** d * (6 * d ** 2 / ((2. + d) * (2. + 2. * d))) ** (d / 2.) * special.gamma(d / 2.)
    I = G * np.exp(-q ** 2 * Rg ** 2 / 3.) + C / q ** d * (special.erf(q * Rg / 6 ** 0.5)) ** (3 * d)
    I[q == 0] = 1
    return I

beaucage_dict = {
    'Rg': {
        'unit': 'nm',
        'tex': '$R_g$',
        'p0' : 1,
        'bounds': [0,10]
    },
    'scale' : {
        'unit' : '',
        'tex' : '$A$',
        'p0' : 0.0002,
        'bounds': [0, 1]
    },
    'porod_exp' : {
        'unit':'',
        'tex': 'n',
        'p0':3,
        'bounds':[1, 6]
    }
}

class Beaucage(sasmodel.SAXS_Model):
    
    def __init__(self):
        self.name = 'beaucage'
        self.longname = 'Beaucage'
        self.parameters = ['Rg', 'scale', 'porod_exp']
        self.pdict = {'xi' : ['$\\xi$', 'nm']}
        self.fitfunc = beaucage

    def get_text(self, fit_dict):
        text = """
        $R_g =$ {0:.2f} $\\pm$ {1:.2f} nm
        $m =$ {2:.2f} $\\pm$ {3:.2f}
        $A =$ {4:.2E}
        $\\chi^2 = $ {6:.4}
        """.format(fit_dict['Rg'], fit_dict['std_Rg'],
                   fit_dict['porod_exp'], fit_dict['std_porod_exp'], 
                   fit_dict['scale'], fit_dict['std_scale'],
                  fit_dict['chi2'])
        return text
        

background = sasmodel.Background()
beaucage = Beaucage().plusmodel(background)
def get_text_be(fit_dict):
        text = """
        $R_g =$ {0:.2f} $\\pm$ {1:.2f} nm
        $n =$ {2:.2f} $\\pm$ {3:.2f}
        $G =$ {4:.2E}
        $\\chi^2 = $ {5:.4}
        """.format(fit_dict['Rg'], fit_dict['std_Rg'],
                   fit_dict['porod_exp'], fit_dict['std_porod_exp'], 
                   fit_dict['scale'],
                  fit_dict['chi2'])
        return text
beaucage.get_text = get_text_be

forward = sasmodel.Porod()
forward.parameters = ['fw_scale', 'fw_exp']

beaucage_forward = beaucage.plusmodel(forward, name='beaucage_fw', longname = 'Beaucage')
def get_text_bef(fit_dict):
        text = """
        $R_g =$ {0:.2f} $\\pm$ {1:.2f} nm
        $n =$ {2:.2f} $\\pm$ {3:.2f}
        $m_F = $ {4:.2f} $\\pm$ {5:.2f}
        $G =$ {6:.2E}
        $A_F = $ {7:.2E}
        $bg = $ {9:.2E}
        $\\chi^2 = $ {8:.4}
        """.format(fit_dict['Rg'], fit_dict['std_Rg'],
                   fit_dict['porod_exp'], fit_dict['std_porod_exp'], 
                   fit_dict['fw_exp'], fit_dict['std_fw_exp'],
                   fit_dict['scale'], fit_dict['fw_scale'],
                  fit_dict['chi2'],
                  fit_dict['bg'])
        return text
beaucage_forward.get_text = get_text_bef
