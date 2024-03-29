import getpass

from .. import database_operations as dbo
import pandas as pd
import numpy as np
import datetime as dt
import os
from os.path  import join
import matplotlib.pyplot as plt
from scipy import optimize, constants
from ..Sample import Sample
from .. import plot_utility as plu
# from .DLS_Measurement import DL
from .. import os_utility as osu
from ..Measurement import Measurement

from matplotlib.backends.backend_pdf import PdfPages
from .repes_utility import repes


def qlabel(ax, par):
    parstring = par
    ax.set_xlabel('$q\\,\\mathrm{[nm^{-1}]}$')
    ax.set_ylabel(f'${parstring}\\,\\mathrm{{[A.U.]}}$')


def qqlabel(ax, par, unit='A.U.'):
    parstring = par
    ax.set_xlabel('$q^2\\,\\mathrm{[\mu m^{-2}]}$')
    ax.set_ylabel(f'${parstring}\\,\\mathrm{{[{unit}]}}$')
    if unit is None:
        ax.set_ylabel(f'${parstring}$')


def make_plot(kwargs):
    if 'ax' not in kwargs:
        fig, ax = plt.subplots()
    else:
        ax = kwargs.pop('ax')
    return ax, kwargs


def constant_fit(x, A):
    return A*x/x

def linear(x, A, B):
    return A*x + B

def quadratic(x, A, B, C):
    return A*x**2 + B*x + C

def quadratic_intercept0(x, A, B):
    return A*x**2 + B*x

def get_qqfit(xdata, ydata, fitfunc=linear, bounds=[0, np.inf], **kwargs):
    xdata_temp = xdata  # * 1e-14
    popt, pcov = optimize.curve_fit(fitfunc, xdata_temp, ydata,
                                    #  sigma=df.err_g2,
                                    bounds=bounds)
    # popt[0] = popt * 1e14
    xdata = np.append(np.array([0]), xdata)
    return popt, pcov

def plot_qqfit(xdata, ydata, fitfunc=linear, bounds=[0, np.inf], **kwargs):
    ax, kwargs = make_plot(kwargs)
    xdata_temp = xdata  # * 1e-14
    popt, pcov = optimize.curve_fit(fitfunc, xdata_temp, ydata,
                                    #  sigma=df.err_g2,
                                    bounds=bounds)
    # popt[0] = popt * 1e14
    xdata = np.linspace(0, xdata[len(xdata)-1], num=1000)
    ax.plot(xdata, fitfunc(xdata, *popt), **kwargs)
    return popt, pcov


def plot_par(m, par, fit=None, ax=None, **kwargs):
    ax, kwargs = make_plot(kwargs)
    df_sls = m.get_sls()
    # df, dfs = m.get_average_g2(phi)
    ax.plot(df_sls.q, df_sls.Isample, **kwargs)


def rayley(m1, m2, m3, times=1, fit=repes, useq2=False,
           qsquared=False, logI=False, **kwargs):
    ax, kwargs = make_plot(kwargs)
    # df = m1.get_rayley_ratio(m_buffer, m_toluene)
    df1 = m1.get_sls()
    df2 = m2.get_sls()
    df3 = m3.get_sls()
    dfbeta = fit.get_phidlstable(m1)
    # print(dfbeta)
    df1mod = df1[df1.q.isin(df2.q)].set_index('angle')
    df2mod = df2[df2.q.isin(df1.q)].set_index('angle')
    # print(df2mod.index)
    df3mod = df3[df3.q.isin(df2mod.q)].set_index('angle')
    dfbetamod = dfbeta[dfbeta.phi.isin(df2mod.index)].set_index('phi')
    beta = dfbetamod.beta
    if 'set_beta' in kwargs:
        set_beta = kwargs.pop( 'set_beta' )
        beta = np.ones(beta.shape)*set_beta
    if 'qunit' in kwargs:
        unit = kwargs.pop( 'qunit' )
        if unit == 'nm':
            xdat = df3mod.q / 1e9 # nm
        elif unit == 'mum':
            xdat = df3mod.q / 1e6 # nm
    else:
        xdat = df3mod.q
    if qsquared:
        xdat = xdat**2
    ydat = beta*times*0.01 *(df1mod.Isample - df2mod.Isample) / df3mod.Isample
    err_ydat = ydat * (df1mod.err_Isample/df1mod.Isample + df2mod.err_Isample/df2mod.Isample + df3mod.err_Isample/df3mod.Isample)
    if logI:
        percent_err = err_ydat / ydat
        ydat = np.log(ydat)
        err_ydat = ydat * percent_err
    if useq2:
        xdat=xdat**2 *1e6 #mum
    ax.errorbar(xdat, ydat,yerr=err_ydat, **kwargs)
    return ax, (xdat, ydat, err_ydat)


def guinier(m1, m2, m3, **kwargs):
    ax, kwargs = make_plot(kwargs)
    df1 = m1.get_sls()
    df2 = m2.get_sls()
    df3 = m3.get_sls()
    df1mod = df1[df1.q.isin(df2.q)].set_index('angle')
    df2mod = df2[df2.q.isin(df1.q)].set_index('angle')
    df3mod = df3[df3.q.isin(df2mod.q)].set_index('angle')
    rayley_ratio = (df1mod.Isample - df2mod.Isample) / df3mod.Isample
    xdat = df3mod.q ** 2 / 1e18
    ydat = rayley_ratio
    ax.plot(xdat, ydat, **kwargs)
    ax.set_yscale('log')
    return ax


def Gamma(m, fit='repes', rmin=0, rmax=np.inf, **kwargs):
    ax, kwargs = make_plot(kwargs)
    df = fit.get_phidlstable(m, rmin=rmin, rmax=rmax)
    xdata, ydata = df.qq * 1e-12, df.Gamma
    ax.errorbar(xdata, df.Gamma, **kwargs)
    qqlabel(ax, '\\Gamma', unit='1/s')
    return ax, (xdata, ydata)


def Dapp(m, fit=None, rmin=0, rmax=np.inf, **kwargs):
    ax, kwargs = make_plot(kwargs)
    df = fit.get_phidlstable(m, rmin=rmin, rmax=rmax)
    xdata, ydata = df.qq * 1e-12, df.Dapp * 1e12
    ax.errorbar(xdata, ydata, **kwargs)
    qqlabel(ax, 'D_{{app}}', unit='\\mu m^2/s')
    return ax, (xdata, ydata)


def Rapp(m, fit='repes', rmin=0, rmax=np.inf, **kwargs):
    ax, kwargs = make_plot(kwargs)
    df = fit.get_phidlstable(m, rmin=rmin, rmax=rmax)
    xdata, ydata = df.qq * 1e-12, df.Rapp * 10**9
    ax.errorbar(xdata, ydata, **kwargs)
    qqlabel(ax, 'R_{{app}}', unit='nm')
    return ax, (xdata, ydata)


def beta(m, fit=repes, **kwargs):
    ax, kwargs = make_plot(kwargs)
    df = fit.get_phidlstable(m)
    xdata, ydata = df.qq * 1e-12, df.beta
    ax.errorbar(xdata, ydata, **kwargs)
    qqlabel(ax, '\\beta', unit=None)
    return ax, (xdata, ydata)


