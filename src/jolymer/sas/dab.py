from . import SAXS_Model as sasmodel
import numpy as np
from scipy import special


def dab_function(q, scale, cor_length):
    out = scale * cor_length**3 / (1 + q*q*cor_length*cor_length)**2
    return out

class DAB(sasmodel.SAXS_Model):

    def __init__(self):
        self.name = 'dab'
        self.longname = 'Debye Anderson Brumberger'
        self.parameters = ['dab_scale', 'dab_cor_length']
        self.fitfunc = dab_function

    def get_text(self, fit_dict):
        text = """
        $\\xi = {.3f} \\pm {.3f}$ nm
        $A_{{dab}} = {.3f}
        """.format(fit_dict['dab_cor_length'], fit_dict['std_dab_cor_length'],
                fit_dict['dab_scale'])

dab = DAB()
