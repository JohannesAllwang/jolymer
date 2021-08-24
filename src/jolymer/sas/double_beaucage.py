# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 10:16:44 2021

@author: xcill
"""

from . import SAXS_Model as sasmodel
import numpy as np
from scipy import special

def double_beaucage(q, Rg2, Rg3, G2, G3, d1, d2, d3):
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
    Rg2 = float(Rg2)
    Rg3 = float(Rg3)
    Rcut1 = Rg2
    Rcut2 = Rg3
    G1 = 1e-9
    B2 = G2 * d2 * Rg2**(-d2) * (6 * d2**2 * float(2+d2)**-1 * float(2 + 2*d2)**-1)**(d2/2)
    B3 = G3 * d3 * Rg3**(-d3) * (6 * d3**2 * float(2+d3)**-1 * float(2 + 2*d3)**-1)**(d3/2)
    I1 = G1 * np.exp(-q**2 * Rcut1**2 / 3) * q**-d1
    I2 = G2 * np.exp(-q**2 * Rg2**2 / 3) + B2 * np.exp(-q**2 * Rcut2**2 / 3) * special.erf(q * Rg2 / np.sqrt(6))**(3*d2) * q**(-d2)
    I3 = G3 * np.exp(-q**2 * Rg3**2 / 3) + B3 *  special.erf(q * Rg3 / np.sqrt(6))**(3*d3) * q**(-d3)
    I = I1 + I2 + I3
    I[q == 0] = 1
    return I

beaucage_dict = {
    'Rg2': {
        'unit': 'nm',
        'tex': '$R_g$',
        'p0' : 1,
        'bounds': [0,1000]
    },
    'Rg3': {
        'unit': 'nm',
        'tex': '$R_g$',
        'p0' : 1,
        'bounds': [0,10]
    },
    'G2' : {
        'unit' : '',
        'tex' : '$G_2$',
        'p0' : 0.0002,
        'bounds': [0, 1]
    },
    'd1' : {
        'unit':'',
        'tex': 'd1',
        'p0':2,
        'bounds':[1, 6]
    },
    'd2' : {
        'unit':'',
        'tex': 'd2',
        'p0':3,
        'bounds':[1, 6]
    },
    'd3' : {
        'unit':'',
        'tex': 'd3',
        'p0':4,
        'bounds':[1, 6]
    }
}

class DoubleBeaucage(sasmodel.SAXS_Model):
    
    def __init__(self):
        self.name = 'beaucage'
        self.longname = 'Beaucage'
        self.parameters = ['Rg2',  'Rg3', 'G2', 'G3', 'd1', 'd2', 'd3']
        self.pdict = {'xi' : ['$\\xi$', 'nm']}
        self.fitfunc = double_beaucage

    def get_text(self, fit_dict):
        text = """
        $R_g2 =$ {0:.2f} $\\pm$ {1:.2f} nm
        $R_g3 =$ {2:.2f} $\\pm$ {3:.2f} nm
        $d_1 =$ {4:.2f} $\\pm$ {5:.2f}
        $d_2 =$ {6:.2f} $\\pm$ {7:.2f}
        $d_3 =$ {8:.2f} $\\pm$ {9:.2f}
        $G_2 =$ {10:.2E}
        $G_3 =$ {11:.2E}
        $\\chi^2 = $ {12:.4}
        """.format(fit_dict['Rg2'], fit_dict['std_Rg2'],
               fit_dict['Rg3'], fit_dict['std_Rg3'],
               fit_dict['d1'], fit_dict['std_d1'], 
               fit_dict['d2'], fit_dict['std_d2'], 
               fit_dict['d3'], fit_dict['std_d3'], 
               fit_dict['G2'], 
               fit_dict['G3'],
               fit_dict['chi2'])
        return text
        

background = sasmodel.Background()
dbeaucage = DoubleBeaucage().plusmodel(background)

# background = sasmodel.Background()
# beaucage = Beaucage().plusmodel(background)
# def get_text_be(fit_dict):
#         text = """
#         $R_g =$ {0:.2f} $\\pm$ {1:.2f} nm
#         $n =$ {2:.2f} $\\pm$ {3:.2f}
#         $C =$ {4:.2E}
#         $\\chi^2 = $ {5:.4}
#         """.format(fit_dict['Rg'], fit_dict['std_Rg'],
#                    fit_dict['porod_exp'], fit_dict['std_porod_exp'], 
#                    fit_dict['scale'],
#                   fit_dict['chi2'])
#         return text
# beaucage.get_text = get_text_be

# forward = sasmodel.Porod()
# forward.parameters = ['fw_scale', 'fw_exp']

# beaucage_forward = beaucage.plusmodel(forward, name='beaucage_fw', longname = 'Beaucage')
# def get_text_bef(fit_dict):
#         text = """
#         $R_g =$ {0:.2f} $\\pm$ {1:.2f} nm
#         $n =$ {2:.2f} $\\pm$ {3:.2f}
#         $m_F = $ {4:.2f} $\\pm$ {5:.2f}
#         $C =$ {6:.2E}
#         $A_F = $ {7:.2E}
#         $\\chi^2 = $ {8:.4}
#         """.format(fit_dict['Rg'], fit_dict['std_Rg'],
#                    fit_dict['porod_exp'], fit_dict['std_porod_exp'], 
#                    fit_dict['fw_exp'], fit_dict['std_fw_exp'],
#                    fit_dict['scale'], fit_dict['fw_scale'],
#                   fit_dict['chi2'])
#         return text
# beaucage_forward.get_text = get_text_bef
