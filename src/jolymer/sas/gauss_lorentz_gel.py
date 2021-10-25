r"""
This model calculates the scattering from a gel structure,
but typically a physical rather than chemical network.
It is modeled as a sum of a low-q exponential decay (which happens to
give a functional form similar to Guinier scattering, so interpret with
care) plus a Lorentzian at higher-q values. See also the gel_fit model.

Definition
----------

The scattering intensity $I(q)$ is calculated as (Eqn. 5 from the reference)

.. math:: I(q) = I_G(0) \exp(-q^2\Xi ^2/2) + I_L(0)/(1+q^2\xi^2)

$\Xi$ is the length scale of the static correlations in the gel, which can
be attributed to the "frozen-in" crosslinks. $\xi$ is the dynamic correlation
length, which can be attributed to the fluctuating polymer chains between
crosslinks. $I_G(0)$ and $I_L(0)$ are the scaling factors for each of these
structures. Think carefully about how these map to your particular system!

.. note::
    The peaked structure at higher $q$ values (Figure 2 from the reference)
    is not reproduced by the model. Peaks can be introduced into the model
    by summing this model with the :ref:`gaussian-peak` model.

For 2D data the scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math:: q = \sqrt{q_x^2 + q_y^2}

References
----------

#. G Evmenenko, E Theunissen, K Mortensen, H Reynaers,
   *Polymer*, 42 (2001) 2907-2913

Authorship and Verification
----------------------------

* **Author:**
* **Last Modified by:**
* **Last Reviewed by:**
"""

import numpy as np
from . import SAXS_Model as sasmodel
import numpy as np
from scipy import special

name = "gauss_lorentz_gel"
title = "Gauss Lorentz Gel model of scattering from a gel structure"
description = """
            Class that evaluates a GaussLorentzGel model.

            I(q) = scale_g*exp(- q^2*Z^2 / 2)+scale_l/(1+q^2*z^2)
                    + background
            List of default parameters:
                scale_g = Gauss scale factor
                Z = Static correlation length
                scale_l = Lorentzian scale factor
                z = Dynamic correlation length
                background = Incoherent background
            """
category = "shape-independent"
# pylint: disable=bad-whitespace, line-too-long
#            ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["gauss_scale",   "",    100.0,  [-np.inf, np.inf], "", "Gauss scale factor"],
              ["cor_length_static",    "Ang", 100.0,  [0, np.inf],    "", "Static correlation length"],
              ["lorentz_scale", "",     50.0,  [-np.inf, np.inf], "", "Lorentzian scale factor"],
              ["cor_length_dynamic",   "Ang",  20.0,  [0, np.inf],    "", "Dynamic correlation length"],
             ]
# pylint: enable=bad-whitespace, line-too-long

def gauss_function(q, gauss_scale, cor_length_static):
    """

    :param q:                    Input q-value
    :param gauss_scale:   Gauss scale factor
    :param cor_length_static:    Static correlation length
    :param lorentz_scale: Lorentzian scale factor
    :param cor_length_dynamic:   Dynamic correlation length
    :return:                     1-D intensity
    """

    term1 = gauss_scale *\
            np.exp(-1.0*q*q*cor_length_static*cor_length_static/2.0)

    return term1

class Gauss(sasmodel.SAXS_Model):

    def __init__(self):
        self.name = 'gauss'
        self.longname = 'Gauss'
        self.parameters = ['gauss_scale', 'cor_length_static']
        self.fitfunc = gauss_function

    def get_text(self, fit_dict):
        text = """
        $\\Xi$ = {.3f} $\\pm$ {.3f}
        $A_G$ = {.3f}
        """.format(fit_dict['cor_length_static'], fit_dict['std_cor_length_static'],
                fit_dict['gauss_scale'])
        return text


def lorentz_function(q, lorentz_scale, cor_length_dynamic):
    term2 = lorentz_scale /\
            (1.0+(q*cor_length_dynamic)*(q*cor_length_dynamic))
    return term2


class Lorentz(sasmodel.SAXS_Model):

    def __init__(self):
        self.name = 'lorentz'
        self.longname = 'Lorentz'
        self.parameters = ['lorentz_scale', 'cor_length_dynamic']
        self.fitfunc = lorentz_function

    def get_text(self, fit_dict):
        text = """
        $\\xi = {.3f} \\pm {.3f}$ nm
        $A_L = {.3f}$
        """.format(fit_dict['cor_length_dynamic'], fit_dict['std_cor_length_dynamic'], 
                fit_dict['lorentz_scale'])
        return text

def Iq(q, gauss_scale, cor_length_static, lorentz_scale, cor_length_dynamic):
    return gauss_function(q, gauss_scale, cor_length_static) + lorentz_function(q, lorentz_scale, cor_length_dynamic)


class GaussLorentzGel(sasmodel.SAXS_Model):
    
    def __init__(self):
        self.name = 'gauss_lorentz_gel'
        self.longname = 'Gauss Lorentz Gel'
        self.parameters = ['gauss_scale', 'cor_length_static' ,'lorentz_scale', 'cor_length_dynamic']
        self.pdict = {'xi' : ['$\\xi$', 'nm']}
        self.fitfunc = Iq

    def get_text(self, fit_dict):
        text = """
        $\\zeta_s =$ {0:.2f} $\\pm$ {1:.2f} nm
        $\\zeta_d =$ {2:.2f} $\\pm$ {3:.2f} nm
        $A_{{Gauss}} =$ {4:.2E}
        $A_{{Lorentz}} =$ {6:.2E}
        $\\chi^2 = $ {8:.4}
        """.format(fit_dict['cor_length_static'], fit_dict['std_cor_length_static'],
                   fit_dict['cor_length_dynamic'], fit_dict['std_cor_length_dynamic'], 
                   fit_dict['gauss_scale'], fit_dict['std_gauss_scale'],
                   fit_dict['lorentz_scale'], fit_dict['std_lorentz_scale'],
                  fit_dict['chi2'])
        return text

lorentz = Lorentz()
gauss = Gauss()

gauss_lorentz = GaussLorentzGel()
        
bg = sasmodel.Background()
gauss_lorentz_bg = GaussLorentzGel().plusmodel(bg)

forward = sasmodel.Porod()
forward.parameters = ['fw_scale', 'fw_exp']
fw_gauss_lorentz = gauss_lorentz.plusmodel(forward, name='fw_gauss_lorentz', longname = 'Gauss Lorentz Gel Forward')
fw_gauss_lorentz_bg = gauss_lorentz_bg.plusmodel(forward, name='fw_gauss_lorentz_bg', longname = 'Gauss Lorentz Gel Forward')


