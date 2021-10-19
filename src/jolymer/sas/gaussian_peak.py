r"""

Definition
----------

This model describes a Gaussian shaped peak on a flat background

.. math::

    I(q) = (\text{scale}) \exp\left[ -\tfrac12 (q-q_0)^2 / \sigma^2 \right]
        + \text{background}

with the peak having height of *scale* centered at $q_0$ and having a standard
deviation of $\sigma$. The FWHM (full-width half-maximum) is $2.354 \sigma$.

For 2D data, scattering intensity is calculated in the same way as 1D,
where the $q$ vector is defined as

.. math::

    q = \sqrt{q_x^2 + q_y^2}


References
----------

None.

Authorship and Verification
----------------------------

* **Author:**
* **Last Modified by:**
* **Last Reviewed by:**
"""

from . import SAXS_Model as sasmodel
import numpy as np
from scipy import special
from .ornsteinZernike import oz

name = "gaussian_peak"
title = "Gaussian shaped peak"
description = """
    Model describes a Gaussian shaped peak including a flat background
    Provide F(q) = scale*exp( -1/2 *[(q-peak_pos)/sigma]^2 )+ background
"""
category = "shape-independent"

#             ["name", "units", default, [lower, upper], "type","description"],
parameters = [["peak_pos", "1/Ang", 0.05, [-np.inf, np.inf], "", "Peak position"],
              ["sigma", "1/Ang", 0.005, [0, np.inf], "",
               "Peak width (standard deviation)"],
             ]

def Iq(q, peak_pos, sigma, scale):
    scaled_dq = (q - peak_pos*np.ones(len(q)))/sigma
    return scale * np.exp(-0.5*scaled_dq*scaled_dq)

class GaussianPeak(sasmodel.SAXS_Model):
    
    def __init__(self):
        self.name = 'gaussian_peak'
        self.longname = 'Gaussian Peak'
        self.parameters = ['peak_pos', 'sigma', 'scale']
        self.pdict = {'xi' : ['$\\xi$', 'nm']}
        self.fitfunc = Iq

    def get_text(self, fit_dict):
        text = """
        $R =$ {0:.4f} $\\pm$ {1:.2f} nm
        $\\sigma =$ {2:.2f} $\\pm$ {3:.2f}
        $A =$ {4:.2E}
        $\\chi^2 = $ {6:.4}
        """.format(fit_dict['peak_pos'], fit_dict['std_peak_pos'],
                   fit_dict['sigma'], fit_dict['std_sigma'], 
                   fit_dict['scale'], fit_dict['std_scale'],
                  fit_dict['chi2'])
        return text
        


background = sasmodel.Background()
gaussian_peak = GaussianPeak()
# gaussian_peak.get_text = get_text_be

forward = sasmodel.Porod()
forward.parameters = ['fw_scale', 'fw_exp']

gaussian_peak_forward = gaussian_peak.plusmodel(forward, name='gaussian_peak_fw', longname = 'Beaucage')

def get_text_gel(fit_dict):
    text = """
    $R =$ {0:.5f} $\\pm$ {1:.2f} nm
    $\\sigma =$ {2:.2f} $\\pm$ {3:.2f}
    $\\xi =$ {4:.2f} $\\pm$ {5:.2f}
    $\\n_oz =$ {6:.2f} $\\pm$ {7:.2f}
    $A =$ {8:.2E}
    $A =$ {10:.2E}
    $\\chi^2 = $ {12:.4}
    """.format(fit_dict['peak_pos'], fit_dict['std_peak_pos'],
               fit_dict['sigma'], fit_dict['std_sigma'], 
               fit_dict['xi'], fit_dict['std_xi'], 
               fit_dict['lorentz_exp'], fit_dict['std_lorentz_exp'], 
               fit_dict['scale'], fit_dict['std_scale'],
               fit_dict['lorentz_scale'], fit_dict['std_lorentz_scale'],
              fit_dict['chi2'])
    return text
gel_model = gaussian_peak_forward.plusmodel(oz, name = 'get_model', longname = 'forward gaussian peak and oz')
gel_model.get_text = get_text_gel

