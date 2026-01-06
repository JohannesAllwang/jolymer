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

    def __init__(self, name, usename=None, samplename=None):
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

    def get_Experiment(self, m, path=None):
        if path is None:
            path = m.get_filename()
        mdata = data.load_data(path)
        return Experiment(mdata, self.sasmodel)

    def load_fit(self, m, datapath=None, errpath=None, name=None):
        fit_dict = {}
        if errpath is None:
            path = join(m.path, f'{name}.err')
        if datapath is None:
            path = join(m.path, f'{name}.dat')
        if name is None:
            name = self.name
        M = self.get_Experiment(m, path=datapath)
        print('load_fit', errpath)
        with open(errpath) as f:
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

    def load_json_fit(self, m, datapath=None, jsonpath=None, name=None):
        """
        Load DREAM JSON results and apply .best values and std ~ p68 half-width.
        Returns (experiment, fit_dict) similar to load_fit().
        """
        import json
        if name is None:
            name = self.name

        # Determine paths
        if datapath is None:
            datapath = join(m.path, f"{name}.dat")
        if jsonpath is None:
            jsonpath = join(m.path, f"{name}.json")

        # Create experiment/model
        M = self.get_Experiment(m, path=datapath)

        # Container for returned values
        fit_dict = {}

        # Load json file
        with open(jsonpath, "r") as f:
            results = json.load(f)

        # Loop through JSON parameters
        for par, stats in results.items():
            best = stats.get("best", None)
            p68 = stats.get("p68", None)
            mean = stats.get("mean", None)

            if best is not None:
                # update model value if exists
                try:
                    M.parameters()[par].value = float(best)
                except Exception:
                    pass

                fit_dict[par] = float(best)

            # standard deviation approximation: half width of 68% credible interval
            if p68 is not None and len(p68) == 2:
                std = (float(p68[1]) - float(p68[0])) / 2
            else:
                std = stats.get("std", "NA")

            fit_dict[f"std_{par}"] = std

        # chi2 if present
        if "chi2" in results:
            fit_dict["chi2"] = float(results["chi2"])

        return M, fit_dict


    def fit(self, m, iqmin=0, iqmax=np.inf, **kwargs):
        """
        That this kind of model is compatible with the following kind of code:
        m.fit_dict, m.fit_df = m.model.fit(m, bounds=m.bounds,
                                           iqmax=m.iqmax, p0=m.p0,
                                           iqmin=m.iqmin,
                                           fixed_parameters=m.fixed_pars)
        """
        if 'errpath' in kwargs:
            errpath = kwargs.pop('errpath')
        else:
            errpath = None
        if 'datapath' in kwargs:
            datapath = kwargs.pop('datapath')
        else:
            datapath = None
        print('fit', errpath)
        M, fit_dict = self.load_fit(m, errpath=errpath,
                                    datapath=datapath,
                                    name=None)
        fit_df = get_fit_df(M)
        # fit_df = fit_df.loc[fit_df.q < iqmax]
        # fit_df = fit_df.loc[fit_df.q > iqmin]
        # fit_df = fit_df[iqmin:iqmax]
        return fit_dict, fit_df

    def write_bumpsfile(self, datapath, datafile, model='sphere', bumpsinpath=None):
        # datapath = join(datapath, datafile) # defined twice
        if bumpsinpath is None:
            bumpsinpath = join(datapath, f'{self.name}_{datafile}.py')
        bumpsoutpath = join(f'{self.name}_{datafile}')

        with open(bumpsinpath, 'w') as f:
            f.write('import matplotlib.pyplot as plt\n')
            f.write('from bumps import names\n')
            f.write('import subprocess\n')
            f.write('import numpy as np\n')
            f.write('from os.path import join\n')
            f.write('from sasmodels import compare, data, core\n')
            f.write('import sasmodels.compare as sascomp\n')
            f.write('import sasmodels.data as data\n')
            f.write('from sasmodels.bumps_model import Experiment, Model\n')

            f.write(f"model = core.load_model('{model}')\n")
            f.write("model.name = 'modelname'\n")
            f.write("jspheregel = Model(model)\n")

            f.write(f"path = '{datapath}/{datafile}'\n\n")
            f.write("testdata = data.load_data(path)\n")
            f.write("M = Experiment(testdata, jspheregel)\n")
            for par in self.parameters:
                f.write(f'# M.parameters()["{par}"].range(0, np.inf)\n')
                f.write(f'# M.parameters()["{par}"].value = {self.sasmodel.parameters()[par].value}\n')

            f.write('M.plot()\n')
            f.write('plt.show()\n')
            f.write('problem = names.FitProblem(M)\n')
            f.write(f'problem.store = "{self.name}_{datafile.split(".")[0]}"\n')
            f.write('problem.plot()\n')
            f.write('plt.show()\n')

    def load_pars(self):
        pass

    def plot_fit(self, **kwargs):
        df = self.fit(**kwargs)

    def get_text(self, fit_dict):
        pass
