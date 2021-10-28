#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 12:17:59 2020

@author: johannes
"""


from scipy import optimize
import numpy as np
from .. import database_operations as dbo



class SAXS_Model:
    
    def __init__(self, name, longname, parameters, pardict, fitfunc):
        self.name = name
        self.longname = longname
        self.parameters = parameters
        self.fitfunc = fitfunc
        self.pardict = pardict
    
    def fit(self, measurement, X = 'q', Y = 'I', **kwargs):
        """
        fits the model and returns a dictionary

        Parameters
        ----------
        measurement : SAXS_Measurement object
        **kwargs : TYPE
            iqmin/iqmax : index of lowest/highest q to use
            qmin/qmax : fit only regards data where qmin < q < qmax
            bounds : get passed to curve_fit. But give as a dictionary with the parnames as keys
            p0 : as bounds
            fixed_parameters : dict {'parameter' : set_value}

        Returns
        -------
        fit_dict : dictionary
            Contains all the relevant information about the fit.
        df : pandas dataframe:
            q, I, fit, res
        """
        df = measurement.get_data(cout=False)
        iqmin = 0
        iqmax = len(df)
        if 'iqmin' in kwargs:
            iqmin = kwargs['iqmin']
            qmin = df.q[iqmin]
            if 'qmin' in kwargs:
                print('qmin will overwrite iqmin...')
        if 'iqmax' in kwargs:
            iqmax = kwargs['iqmax']
            qmax = df.q[iqmax]
            if 'qmax' in kwargs:
                print('qmax will overwrite iqmax...')
        if 'qmin' in kwargs:
            qmin = kwargs['qmin']
            iqmin = np.where(np.diff(np.sign(df.q - qmin)))[0][0] + 1
        if 'qmax' in kwargs:
            qmax = kwargs['qmax']
            iqmax = np.where(np.diff(np.sign(df.q - qmax)))[0][0]
        
        _pfit = []
        _pfix = []
        if 'fixed_parameters' in kwargs:
            _pfix = kwargs['fixed_parameters']
            fitfunc = self.fix_parameters(_pfix)
            for par in self.parameters:
                if not par in _pfix:
                    _pfit.append(par)
        else:
            _pfit = self.parameters
        
        bounds = (-np.inf, np.inf)
        if 'bounds' in kwargs:
            bounds= [[], []]
            bounds_dict = kwargs['bounds']
            for p in _pfit:
                if p in bounds_dict:
                    bounds[0].append(bounds_dict[p][0])
                    bounds[1].append(bounds_dict[p][1])
                else:
                    bounds[0].append(-np.inf)
                    bounds[1].append(np.inf)
        _p0 = None
        if 'p0' in kwargs:
            _p0=[]
            for par in _pfit:
                if par in kwargs['p0']:
                    _p0.append(kwargs['p0'][par])
                else:
                    _p0.append(0)
            
        df = df[iqmin : iqmax + 1]
        popt, pcov = optimize.curve_fit(fitfunc, df.q, df.I, sigma = df.err_I, 
                                        p0 = _p0, bounds=bounds)
        pstd = np.sqrt(np.diag(pcov))
        normalized_pcov = pcov / np.outer(pstd, pstd)
      #  print('normalized covariance matrix:', normalized_pcov)
        df['fit'] = fitfunc(df.q, *popt)
        df['res'] = df.fit - df.I
        chi2 = np.sum(((df.fit - df.I) / df.err_I)**2 / 
                                 (len(df) - len(_pfit))) 
        fit_dict = {'iqmin':iqmin, 'iqmax':iqmax}
        for i, p in enumerate(self.parameters):
            if p in _pfix:
                fit_dict[p] = _pfix[p]
                fit_dict[f'std_{p}'] = 'fixed'
            else:
                index = _pfit.index(p)
                fit_dict[p] = popt[index]
                fit_dict[f'std_{p}'] = pstd[index]
        fit_dict['chi2'] = chi2
        # fit_dict['normalized_pcov'] = normalized_pcov
        fit_dict['measurement'] = measurement
        return fit_dict, df
    
    
    def fix_parameters(self, parameters):
        func=self.fitfunc
        def inner(*args):
            inner_args = [args[0]]
            i = 1
            for parameter in self.parameters:
                if parameter in parameters:
                    inner_args.append(parameters[parameter])
                else:
                    inner_args.append(args[i])
                    i+=1
            return func(*inner_args)
        return inner
    
    @staticmethod
    def plot_fit(fit_df, figure, scale = 1, **kwargs):
        
        fig, ax = figure
        
        ax.errorbar(fit_df.q, fit_df.fit * scale, **kwargs)

    def plusmodel(self, model, **kwargs):
        "TODO"
        if 'name' in kwargs:
            name = kwargs['name']
        else:
            name = f'{self.name}'
        
        if 'longname' in kwargs:
            longname = kwargs['longname']
        else:
            longname = f'{self.longname}'
        
        for par in self.parameters:
            if par in model.parameters:
                raise Exception("Model parameters need to have different names.")
        new_parameters = self.parameters.copy()
        len_self_pars = len(self.parameters)
        len_model_pars = len(model.parameters)
        new_pardict = {} #TODO
        for par in model.parameters:
            new_parameters.append(par)
        self_fitfunc = self.fitfunc
        model_fitfunc = model.fitfunc
        def new_fitfunc(*args):
            q = args[0]
            args_oldfunc = args[1:len_self_pars + 1]
            args_newfunc = args[1 + len_self_pars::]
            return self_fitfunc(q, *args_oldfunc) + model_fitfunc(q, *args_newfunc)
        self_get_text = self.get_text
        model_get_text = model.get_text
        def new_get_text(fit_dict):
            try:
                return self_get_text(fit_dict) + model_get_text(fit_dict)
            except:
                print('Something is wrong with new_get_text')
                return model_get_text(fit_dict)
        newmodel = SAXS_Model(name, longname, new_parameters, new_pardict, new_fitfunc)
        newmodel.get_text = new_get_text
        return newmodel

    def get_text(self, fit_dict={}):
        return ''

class Background(SAXS_Model):

    def __init__(self):
        self.name  = 'background'
        self.longname = 'Background'
        self.parameters = ['bg']
        self.fitfunc = lambda q, const: const*np.ones(len(q))

    def get_text(self, fit_dict):
        text = """
        bg = {:.3f}
        """.format(fit_dict['bg'])
        return text

class Porod(SAXS_Model):

    def __init__(self):
        self.name = 'porod'
        self.longname = 'Porod'
        self.parameters = ['porod_scale', 'porod_exp']
        self.fitfunc = lambda q, A, m: A*np.array(q)**-m

    def get_text(self, fit_dict):
        text = """
        $n = {:.3f} \\pm {:.3f}$
        $A = {:.3f}$
        """.format(fit_dict['porod_exp'], fit_dict['std_porod_exp'],
                fit_dict['porod_scale'])
        return text

porod = Porod()
fw = Porod()
fw.parameters = ['fw_scale', 'fw_exp']
bg = Background()

def par_from_measurement():
    pass
    
def combine_fitresults(fitdicts):
    out = {}
    for fitdict in fitdicts:
        for key in fitdict:
            if key in out:
                out[key].append(fitdict[key])
            else:
                out[key] = [fitdict[key]]
    return out

def treated_untreated(fitdicts):
    tlist = [dic for dic in fitdicts if dic['measurement'].sample.istreated]
    ulist = [dic for dic in fitdicts if not dic['measurement'].sample.istreated]
    treated = combine_fitresults(tlist)
    untreated = combine_fitresults(ulist)
    return treated, untreated

class AddSaxsModels(SAXS_Model):

    def __init__(self, models):
        self.models = models
        name = ''
        parameters = []
        for model in self.models:
            model.parloc = list( range(
                len(parameters) + 1, len(parameters) + len(model.parameters) + 1
                ) )
            parameters += model.parameters
            toname = model.name if name=='' else f'_{model.name}'
            name += toname
        self.parameters = parameters
        self.name = name
        self.longname = name

    def fitfunc(self, *args):
        q = args[0]
        out = np.zeros(len(q))
        for model in self.models:
            arguments = [args[i] for i in model.parloc]
            out += model.fitfunc(q, *arguments)
        return out

    def get_text(self, fit_dict):
        text = ''
        for model in self.models:
            text += model.get_text(fit_dict)
        return text



def add_models(models):
    name = ''
    parameters = []
    for model in models:
        parameters += model.parameters
        toname = model.name if name=='' else f'_{model.name}'
        name += toname

        new_parameters = self.parameters.copy()
        len_self_pars = len(self.parameters)
        len_model_pars = len(model.parameters)
        for par in model.parameters:
            new_parameters.append(par)
        self_fitfunc = self.fitfunc
        model_fitfunc = model.fitfunc
        def new_fitfunc(*args):
            q = args[0]
            args_oldfunc = args[1:len_self_pars + 1]
            args_newfunc = args[1 + len_self_pars::]
            return self_fitfunc(q, *args_oldfunc) + model_fitfunc(q, *args_newfunc)
        def new_get_text(fit_dict):
            out = ''
            for model in models:
                out += model.get_text(fit_dict)
            return out
        newmodel = SAXS_Model(name, longname, new_parameters, new_pardict, new_fitfunc)
        newmodel.get_text = new_get_text
    return outmodel
