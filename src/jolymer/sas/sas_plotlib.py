#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 11:37:48 2021

@author: johannes
"""


from .. import database_operations as dbo
import pandas as pd
import numpy as np
import datetime as dt
import os
import matplotlib.pyplot as plt
from scipy import optimize, constants
from ..Sample import Sample
from .. import plot_utility as plu
# from .DLS_Measurement import DL

from matplotlib.backends.backend_pdf import PdfPages

# plt.style.use('classic')

def _fit_compilation(m, fit, seq_numbers, cm, labelfunc, title, ax=None, **kwargs):
    citer = plu.cm_for_l(seq_numbers, cm)
    for s, color in zip(seq_numbers, citer):
        label=labelfunc(m, fit, s)
        ax = m.plot_fit(s, fit, ax=ax, color=color, label=label, **kwargs)
    ax.legend()
    ax.set_title(title)
    return ax

def _fit_compilations(listofargs, **kwargs):
    num_plots = len(listofargs)
    fig, axes = plu.n_subplots(num_plots)

    for row in axes:
        row[0].set_ylabel('ylabel')
    for ax in axes[-1]:
        ax.set_xlabel('xlabel')
    axes = axes.flatten()[0:num_plots]
    for args, ax in zip(listofargs, axes):
        ax=_fit_compilation(*args, ax=ax, **kwargs)
    return fig, axes

def kratky_plot():
    pass

def guinier_plot():
    pass

# def _compilation(m, fixed_pars, label_pars, fit, cm, labelfunc, title, ax=None, **kwargs):
#     query = f"""
#     select id from desy_measurements where
#     """
#     for fixp in fixed pars:
#         query += f"{fixp}"

def absolutes(m):
    fig, axes = plt.subplots(nrows=3, figsize=(8,10))
    for ax in axes:
        ax.set_yscale('log')
        ax.set_xscale('log')
    axsub, axsam, axbuf = axes
    m.plot_data(ax=axsub, label=m.get_parameter('parent'))
    # axsub.legend()
    axsub.set_title(m.get_parameter('parent'))

    absdfs = m.get_absolute_dfs()[0]
    for absdf, color in zip(absdfs, plu.cm_for_l('viridis', absdfs)):
        axsam.errorbar(absdf.q, absdf.I, yerr=absdf.err_I, color=color)
    avgdf = m.get_averaged()
    axsam.errorbar(avgdf.q, avgdf.I, yerr=avgdf.err_I, color='r', label='average')
    axsam.set_title('sample absolutes')
    axsam.legend()

    b_absdfs = m.get_absolute_dfs(buf=True)[0]
    for b_absdf, color in zip(b_absdfs, plu.cm_for_l('viridis', b_absdfs)):
        axbuf.errorbar(b_absdf.q, b_absdf.I, yerr=b_absdf.err_I, color=color)
    b_avgdf = m.get_averaged(buf=True)
    axbuf.errorbar(b_avgdf.q, b_avgdf.I, yerr=b_avgdf.err_I, color='r', label='average')
    axbuf.set_title('buffer absolutes')
    axbuf.legend()



