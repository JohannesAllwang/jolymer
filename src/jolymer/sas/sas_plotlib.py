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

plt.style.use('seaborn')

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

def _compilation(m, fixed_pars, label_pars, fit, cm, labelfunc, title, ax=None, **kwargs):
    query = f"""
    select id from desy_measurements where
    """
    for fixp in fixed pars:
        query += f"{fixp}"
