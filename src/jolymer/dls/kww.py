import subprocess
import getpass
import io
import os
from os.path import join
import pandas as pd
import numpy as np
from scipy import optimize, constants
from scipy import special

from .dlsmodel import DLSmodel
from .lsi import lsi
from .. import database_operations as dbo

_name = 'kww'
_pardict = {
    # 'GammaK': ["$\\Gamma_{K}$ $\\mathrm{[1/s}]$"],
    'DK' : ["$D_K$ $\\mathrm{[m^2/s}$"],
    'nu': ['$\\nu$'],
    'beta': ['$\\beta$']
    }


def _create_bounds(m: lsi, qq):
    DKmin = m.DfromR(m.rmax)
    print('lower DK:', DKmin)
    DKmax = m.DfromR(m.rmin)
    print('upper DK:', DKmax)
    numin = 0.2
    numax = 1
    betamin = 0
    betamax = 1.2
    if m.mode == '3dcross':
        betamax = 0.3
    if m.mode == 'mod3d':
        betamax = 0.3
    return ((DKmin, numin, betamin), (DKmax, numax, betamax))

def _create_fitfunc(m, qq):
    def inner(x, DK, nu, beta):
        GammaK = DK * qq
        exponent = (GammaK * x) ** nu
        g1 = np.exp(-exponent)
        g2 = beta * g1**2
        return g2
    return inner


class KWW(DLSmodel):

    def __init__(self, name=_name, create_fitfunc=_create_fitfunc,
                 create_bounds=_create_bounds, pardict=_pardict,
                 use_db=False):
        self.name = name
        self.parameters = pardict.keys()
        self.pardict = pardict
        self.create_fitfunc = create_fitfunc
        self.create_bounds = create_bounds
        self.seq_columns = ['seq_number']
        self.phi_columns = ['phi', 'qq']
        self.use_db = use_db
        for par in self.parameters:
            self.seq_columns.append(par)
            self.seq_columns.append(f'std_{par}')
            self.phi_columns.append(par)
            self.phi_columns.append(f'err_{par}')

    def get_phidlstable(self, m: lsi, phidls=True, **kwargs):
        if self.use_db:
            fitpars = dbo.get_table(self.phidls_tablename(m))
        elif phidls:
            fitpars = pd.read_csv(self.phidls_tablepath(m))
        elif phidls is False:
            fitpars = self.get_seqtable(m)
            fitpars['phi'] = [m.phifromseq(seq) for seq in fitpars.seq_number]
            fitpars['qq'] = [m.qq(m.phifromseq(seq)) for seq in fitpars['seq_number']]
        fitpars['GammaK'] = fitpars.DK * fitpars.qq
        fitpars['tauK'] = 1 / fitpars.GammaK
        fitpars['tau'] = (fitpars.tauK / fitpars.nu) * special.gamma(1/fitpars.nu)
        fitpars['Gamma'] = 1 / fitpars.tau
        fitpars['Dapp'] = fitpars.Gamma / fitpars.qq
        fitpars['Rapp'] = m.RfromD(fitpars.Dapp)
        fitpars['RK'] = m.RfromD(fitpars.DK)
        return fitpars



def c2_create_fitfunc(m, qq):
    def inner(x, Dapp, Dapp2, beta):
        g1 = np.exp(-Dapp*x * qq) * \
            (1 + Dapp2 * qq**2 * x**2 / 2)
        g2 = beta * g1**2
        return g2
    return inner


def cu2_create_bounds(m, phi):
    Dmin = m.DfromR(m.rmax)
    Dmax = m.DfromR(m.rmin)
    # Dmin = 0
    # Dmax = 10000
    D2min = 0
    D2max = Dmax**2 / 10000000
    # D2max = 1000_000_000_000
    betamin = 0
    betamax = 1.2
    if m.mode == '3dcross':
        betamax = 0.5
    if m.mode == 'mod3d':
        betamax = 0.7
    return ((Dmin, D2min, betamin), (Dmax, D2max, betamax))


cu2_pardict = {
    'Dapp': ["$D_{app}$ $\\mathrm{[m^2/s}]$"],
    'varDapp': ['var$(D_{app})$'],
    'beta': ['$\\beta$']
    }

kww = KWW()
