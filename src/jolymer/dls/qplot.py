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


def qlabel(ax, par):
    parstring = par
    ax.set_xlabel('$q\\,\mathrm{[nm^{-1}]}$')
    ax.set_ylabel(f'${parstring}\\,\mathrm{{[A.U.]}}$')


def make_plot(ax, par=None):
    if ax is None:
        fig, ax = plt.subplots()
        qlabel(ax, par)


def plot_par(m, par, fit=None, ax=None, **kwargs):
    make_plot(ax, par=par)
    df_sls = m.get_sls()
    # df, dfs = m.get_average_g2(phi)
    ax.plot(df_sls.q, df_sls.Isample, **kwargs)


def plot_IminI(m1, m2, m3, ax=None, **kwargs):
    make_plot(ax)
    df1 = m1.get_sls()
    df2 = m2.get_sls()
    df3 = m3.get_sls()
    df1mod = df1[df1.q.isin(df2.q)].set_index('angle')
    df2mod = df2[df2.q.isin(df1.q)].set_index('angle')
    df3mod = df3[df3.q.isin(df2mod.q)].set_index('angle')
    ydat = (df1mod.Isample - df2mod.Isample) / df3mod.Isample
    ax.plot(df1mod.q, ydat, **kwargs)
    return df1mod, df2mod
