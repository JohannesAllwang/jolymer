r"""
my own plugin model
"""

# from jolymer.sas.plugin_models import jbeaucage

import numpy as np
from numpy import inf, errstate
from scipy import special

name = "beaucage"
title = "Beaucage model with correction"
description = """\
      I(q) = scale_p/pow(q,exponent)+scale_l/
      (1.0 + pow((fabs(q-q_peak)*length_l),exponent_l) )+ background
      List of default parameters:
      G = scale
      exp = Porod exponent
      lorentz_scale = Lorentzian term scaling
      lorentz_length = Lorentzian screening length [A]
      peak_pos = peak location [1/A]
      lorentz_exp = Lorentzian exponent
      background = Incoherent background"""
category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["rg",    "",  1.0e-05, [0, inf], "", "Radius of Gyratoin"],
              ["exp",      "",      3.0, [1, 10], "", "Exponent of power law"],
              ["G",  "",     10.0, [0, inf], "", "Scale factor"],
             ]
# pylint: enable=bad-whitespace, line-too-long

def Iq(q,
       rg=1.0e-5,
       exp=3.0,
       G=10.0):
    """
    :param q:              Input q-value
    :param rg:  Radius of Gyration
    :param exp:      Exponent of power law
    :param G:    Scale factor
    :return:               Calculated intensity
    """
    q = np.atleast_1d(q)
    rg = float(rg)
    C = G * exp / rg ** exp * (6 * exp ** 2 / ((2. + exp) * (2. + 2. * exp))) ** (exp / 2.) * special.gamma(exp / 2.)
    I = G * np.exp(-q ** 2 * rg ** 2 / 3.) + C / q ** exp * (special.erf(q * rg / 6 ** 0.5)) ** (3 * exp)
    I[q == 0] = 1
    return I
Iq.vectorized = True  # Iq accepts an array of q values

def random():
    """Return a random parameter set for the model."""
    pars = dict(
        scale=1,
        G=10**np.random.uniform(-8, -5),
        exp=np.random.uniform(1, 6),
        rg=10**np.random.uniform(0.3, 6),
    )
    return pars
