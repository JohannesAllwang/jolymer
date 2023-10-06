import numpy as np
import subprocess
import getpass
import io
import os
import re
from os.path import join
import pandas as pd

from sasmodels import data, core
# import sasmodels.compare as sascomp
from sasmodels.bumps_model import Experiment, Model

from .sas import SAXS_Model as sasmodel

sasview_bin = f"C:\\Users\\{getpass.getuser()}\\sasview"


def run_program(name, *args):
    cmd = ['powershell', f'./{name}.exe', *args]
    return subprocess.run(cmd, cwd=sasview_bin)


def run_sasview(*args):
    out = run_program('sasview', *args)
    return out


def get_fit_df(M):
    fit_df = pd.DataFrame({
        'q': np.array(M._data.x),
        'I': np.array(M.Iq),
        'err_I': np.array(M.dIq),
        'fit': M.theory(),
        'res': M.theory() - np.array(M.Iq)})
    return fit_df

class SasModel(sasmodel.SAXS_Model):

    def __init__(self, name, usename=None):
        self.sasmodel = Model(core.load_model(name))
        self.name = name
        if usename is None:
            pass
        else:
            self.sasmodel.name = usename
            self.name = usename
        self.longname = name
        self.parameters = list(self.sasmodel.parameters().keys())
        self.pdict = self.sasmodel.parameters()

    def get_Experiment(self, m):
        path = m.get_filename()
        mdata = data.load_data(path)
        return Experiment(mdata, self.sasmodel)

    def load_fit(self, m, name=None):
        fit_dict = {}
        if name is None:
            name = self.name
        M = self.get_Experiment(m)
        with open(join(m.path, f'{name}.err')) as f:
            for line in f:
                if len(line) < 2:
                    continue
                elif line[0] == '.':
                    split1 = line.split('=')
                    par = split1[0].replace('.', '').replace(' ', '')
                    value = split1[1].split('in')[0]
                    value = float(value)
                    try:
                        M.parameters()[par].value = value
                    except: pass
                    fit_dict[par] = value
                    fit_dict[f'std_{par}'] = 'fixed'
                elif line[0] == '[':
                    fit_dict['chi2'] = float(line.split('=')[1].split('(')[0])
                elif line[1].isdigit():
                    split1 = re.split('\s+', line.replace('[', ' ').replace(']', ' '))
                    _, _, par, mean, median, best, i68lower, i68upper,_,_,_= split1
                    std_par = (float(i68upper) - float(i68lower))/2
                    # fit_dict[par] = list()
                    # fit_dict[f'std_{par}'] = list()
                    # fit_dict[par].append(float(best))
                    # fit_dict[f'std_{par}'].append(std_par)
                    fit_dict[par] = float(best)
                    fit_dict[f'std_{par}'] = std_par
        return M, fit_dict

    def fit(self, m, iqmin=0, iqmax=np.inf, **kwargs):
        """
        That this kind of model is compatible with the following kind of code:
        m.fit_dict, m.fit_df = m.model.fit(m, bounds=m.bounds,
                                           iqmax=m.iqmax, p0=m.p0,
                                           iqmin=m.iqmin,
                                           fixed_parameters=m.fixed_pars)
        """
        M, fit_dict = self.load_fit(m, name=None)
        fit_df = get_fit_df(M)
        # fit_df = fit_df.loc[fit_df.q < iqmax]
        # fit_df = fit_df.loc[fit_df.q > iqmin]
        fit_df = fit_df[iqmin:iqmax]
        return fit_dict, fit_df

    def load_pars(self):
        pass

    def plot_fit():
        pass

    def get_text(self, fit_dict):
        pass
