from .. import database_operations as dbo
import pandas as pd
import numpy as np
import datetime as dt
import os
import matplotlib.pyplot as plt
from scipy import optimize, constants
from ..Sample import Sample
from .. import plot_utility as plu
from . import sas_plotlib as splu

from matplotlib.backends.backend_pdf import PdfPages

class Ms:

    def __init__(self, ms):
        self.ms = ms


    def fits(self, fits=True, shiftby=1, **kwargs):
        if 'ax' in kwargs:
            ax = kwargs.pop('ax')
        else:
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (5, 4), 
                                     sharex=True, sharey=True, squeeze=True)
        if 'ylim' in kwargs:
            ylim=kwargs.pop('ylim')
        else:
            ylim = [None, None]
        if 'topleft' in kwargs:
            topleft = kwargs.pop('topleft')
        else:
            topleft=''
        # shiftby = 0.07 if fits else 1
        # shift= 20**len(self.ms)
        shift = shiftby**len(self.ms)
        for m in self.ms:
            m.fit_dict, m.fit_df = m.model.fit(m, bounds=m.bounds, iqmax=m.iqmax,
                                p0=m.p0, iqmin=m.iqmin, fixed_parameters=m.fixed_pars)
            title = f'c({m.sample.PS.short_name})={m.sample.PS_gpl} c(TRY)'
            df = m.get_data(cout=False)[m.iqmin:m.iqmax]
            ax.errorbar(df.q, df.I * shift, df.err_I * shift, marker = m.marker, color=m.color,
                            linestyle='', label = m.label, elinewidth=0.2, **kwargs)
            if fits:
                ax.errorbar(m.fit_df.q, m.fit_df.fit*shift, marker='', color='black')
            shift = shift / shiftby
            m.partext=m.model.get_text(m.fit_dict)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend(fontsize = 'x-small')
        ax.set_ylim(*ylim)
        # ax.grid()
            
        ax.set_xlabel('$q$ [1/nm]')
        ax.set_ylabel('$I$ [1/cm]')
        if fits:
            ax.set_ylabel('Intensity [1/cm] shifted')
            text = '$I(q) = I_{Beaucage}(q) + I_F$'
        ax.annotate(topleft, xy=(.1, .9), xycoords='axes fraction')
    #         ax.annotate(text, xy=(.1, .1), 
    #             xycoords='axes fraction')

    def res(self, **kwargs):
        if 'ax' in kwargs:
            ax = kwargs.pop('ax')
        else:
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (5, 4), 
                                     sharex=True, sharey=True, squeeze=True)
        if 'ylim' in kwargs:
            ylim=kwargs.pop('ylim')
        else:
            ylim = [None, None]
        for m in self.ms:
            label = '{}; $\\chi^2 = {:.1f}$'.format(
                    m.label, m.fit_dict['chi2']
                    )
            df = m.fit_df
            ydata = df.res/df.err_I
            ax.plot(df.q, ydata, marker = m.marker, color=m.color,
                            linestyle='', label = label, **kwargs)
        ax.legend(fontsize = 'xx-small')
        ax.set_ylabel('Normalized Residuals')
        ax.set_xlabel('$q\\mathrm{\,[nm^{-1}]}$')

    def kratky(self, **kwargs):
        if 'ax' in kwargs:
            ax = kwargs.pop('ax')
        else:
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (5, 4), 
                                     sharex=True, sharey=True, squeeze=True)
        if 'ylim' in kwargs:
            ylim=kwargs.pop('ylim')
        else:
            ylim = [None, None]
        for m in self.ms:
            df = m.get_data(cout=False)[m.iqmin:m.iqmax]
            label = ''
            ax.errorbar(df.q, df.q**2 * df.I * 1000,  marker = m.marker, color=m.color,
                            linestyle='', label = label, elinewidth=0.2)
            ax.legend()
            #             ax.annotate(m.partext, xy=(0.0, 0.0), xycoords='axes fraction')
            # ax.grid()
            ax.set_ylim(0,5)
            ax.set_xlim(0, 2.5)
            ax.set_xlabel('$q$ [1/nm]')
            #         axes[0][1].set_xlabel('$q$ [1/nm]')
            ax.set_ylabel('$I\cdot q^2 \mathrm{\,[0.001nm^{-2}cm^{-1}]}$')
            text = 'Kratky plots'
            ax.annotate(text, xy=(.6, .8), xycoords='axes fraction')

    def dfits(self, fits=True, shiftby=1, **kwargs):
        shift = shiftby**len(self.ms)
        if 'ax' in kwargs:
            ax = kwargs.pop('ax')
        else:
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (5, 4), 
                                     sharex=True, sharey=True, squeeze=True)
        if 'ylim' in kwargs:
            ylim=kwargs.pop('ylim')
        else:
            ylim = [None, None]
        if 'topleft' in kwargs:
            topleft = kwargs.pop('topleft')
        else:
            topleft=''
        for m in self.ms:
            df = m.get_data()
            phtext = f'pH {m.sample.buffer.pH}'
            if not m.sample.tt:
                phtext='untreated'
            if fits:
                m.lowq_fitdict, m.lowq_fitdf = m.lmodel.fit(m, p0=m.lp0, qmax=m.q1, bounds=m.lbounds,
                                                         fixed_parameters=m.lfixed_pars)
                m.highq_fitdict, m.highq_fitdf = m.hmodel.fit(m, p0=m.lp0, qmin=m.q2, bounds=m.hbounds,
                                                           fixed_parameters=m.hfixed_pars)
                # label = '{}'.format(
                #     phtext, m.lowq_fitdict['chi2'], m.highq_fitdict['chi2'],
                # )
            
            ax.loglog(df.q, df.I*shift, marker=m.marker, label=m.label, color=m.color, linestyle='', **kwargs)
            if fits:
                ax.loglog(m.highq_fitdf.q, m.highq_fitdf.fit*shift, marker='', color='black')
                ax.loglog(m.lowq_fitdf.q, m.lowq_fitdf.fit*shift, marker='', color='black')
            shift = shift/shiftby
        ax.legend(fontsize='x-small')
        ax.set_xlabel('$q$ [1/nm]')
        ax.set_ylim(*ylim)
        # ax.grid()
    #     axes[0][1].set_xlabel('$q$ [1/nm]')
        ax.set_ylabel('$I$ [1/cm]')
        if fits:
            ax.set_ylabel('Intensity [1/cm] shifted')
            text = '$I(q) = I_{Beaucage}(q) + I_F$'
        ax.annotate(topleft, xy=(.1, .9), xycoords='axes fraction')

    def dres(self, **kwargs):
        if 'ax' in kwargs:
            ax = kwargs.pop('ax')
        else:
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (5, 4), 
                                     sharex=True, sharey=True, squeeze=True)
        if 'ylim' in kwargs:
            ylim=kwargs.pop('ylim')
        else:
            ylim = [None, None]
        for m in self.ms:
            label = '{}; $\\chi^2_l = {:.1f}$; $\\chi_h^2 = {:.1f}$'.format(
                    m.label, m.lowq_fitdict[ 'chi2' ], m.highq_fitdict[ 'chi2' ]
                    )
            df1 = m.lowq_fitdf
            df2 = m.highq_fitdf
            ydata1 = df1.res/df1.err_I
            ydata2 = df2.res/df2.err_I
            ax.plot(df1.q, ydata1, marker = m.marker, color=m.color,
                            linestyle='')
            ax.plot(df2.q, ydata2, marker = m.marker, color=m.color,
                            linestyle='', label = label)
        ax.legend(fontsize='xx-small')
        ax.set_ylabel('Normalized Residuals')
        ax.set_xlabel('$q\\mathrm{\,[nm^{-1}]}$')
