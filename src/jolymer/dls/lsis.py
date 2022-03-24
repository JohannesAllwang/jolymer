from .. import database_operations as dbo
import pandas as pd
import numpy as np
import datetime as dt
import os
import matplotlib.pyplot as plt
from scipy import optimize, constants
from ..Sample import Sample
from .. import plot_utility as plu
from . import dls_plotlib as dp
from . import qplot

from matplotlib.backends.backend_pdf import PdfPages

class Lsis:

    def __init__(self, ms, fit=None):
        self.ms = ms
        self.fit = fit
        # self.model = ms[0].model

    def make_plot(self, **kwargs):
        if 'ax' in kwargs:
            ax = kwargs.pop('ax')
        else:
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4),
                                   sharex=True, sharey=True, squeeze=True)
        if 'ylim' in kwargs:
            ylim = kwargs.pop('ylim')
            ax.set_ylim(*ylim)
        if 'xlim' in kwargs:
            xlim = kwargs.pop('xlim')
            ax.set_xlim(*xlim)
        return ax, kwargs

    def fits(self, phi, **kwargs):
        # philabel = dp.philabel
        # xspace = 't'
        # seq_numbers = [78, 198, 122, 66]
        ax, kwargs = self.make_plot(**kwargs)
        for m in self.ms:
            ax = dp.plot_phidls(m, phi, fit=self.fit, ax=ax, showres=False,
                                color=m.color,
                                label=m.label, marker=m.marker, linestyle='')
        ax.set_ylabel('$g_2^c-1$')
        ax.set_xlabel('$\\tau\\mathrm{\\,[s]}$')
        ax.legend()
        return ax

    def dists(self, phi, **kwargs):
        ax, kwargs = self.make_plot(**kwargs)
        for m in self.ms:
            ax = dp.plot_phidls_dist(m, phi, fit=self.fit, ax=ax,
                                     color=m.color,
                                     label=m.label, marker='')
        ax.set_ylabel('$\\tau A(\\tau)$')
        ax.set_xlabel('$\\tau\\mathrm{\\,[s]}$')
        ax.legend()
        return ax

    def rayleighs(self, **kwargs):
        ax, kwargs = self.make_plot(**kwargs)
        for m in self.ms:
            ax = qplot.plot_IminI(m, m.m_buffer, m.m_toluene, ax=ax, color=m.color,
                                  marker=m.marker, label=m.label)

        ax.legend()
        qplot.qlabel(ax, 'I')
        ax.set_ylabel('$(I_{sample} - I_{H2O})/ I_{toluene}$')
        plu.loglog(ax)

    def guiniers(self, **kwargs):
        ax, kwargs = self.make_plot(**kwargs)
        for m in self.ms:
            ax = qplot.guinier(m, m.m_buffer, m.m_toluene, ax=ax, color=m.color,
                               marker=m.marker, label=m.label)
        ax.legend()
        qplot.qqlabel(ax, 'ln($I$)')
        # ax.set_ylabel('$(I_{sample} - I_{H2O})/ I_{toluene}$')

    def Dapps(self, **kwargs):
        ax, kwargs = self.make_plot(**kwargs)
        for m in self.ms:
            ax, (xdata, ydata) = qplot.Dapp(m, ax=ax, color=m.color,
                                            fit=self.fit, marker=m.marker,
                                            label=m.label, **kwargs)
            qplot.plot_qqfit(xdata, ydata, color=m.color, marker=None, ax=ax,
                             bounds=[[0, 0], [1e9, 6e6]])
        ax.legend()
        return ax

    def Rapps(self, **kwargs):
        ax, kwargs = self.make_plot(**kwargs)
        rdict = {'label': [],
                 'rh': [],
                 'err_rh': []}
        for m in self.ms:
            ax, (xdata, ydata) = qplot.Rapp(m, ax=ax, color=m.color,
                                            fit=self.fit, marker=m.marker,
                                            label=m.label, **kwargs)
            popt, pcov = qplot.plot_qqfit(xdata, ydata, color=m.color,
                                          marker=None, ax=ax,
                                          bounds=[[-1e6, 0], [1e6, 1e3]])
            rh = popt[1]
            pstd = np.sqrt(np.diag(pcov))
            err_rh = pstd[1]
            print(m.label, ' Rh= ', rh, ' pm ', err_rh)
            rdict['label'].append(m.label)
            rdict['rh'].append(rh)
            rdict['err_rh'].append(err_rh)
        ax.legend()
        rdf = pd.DataFrame(rdict)
        return ax, rdf

    def Gammas(self, **kwargs):
        ax, kwargs = self.make_plot(**kwargs)
        for m in self.ms:
            ax, (xdata, ydata) = qplot.Gamma(m, ax=ax, color=m.color,
                                             fit=self.fit, marker=m.marker,
                                             label=m.label, **kwargs)
            qplot.plot_qqfit(xdata, ydata, color=m.color, marker=None, ax=ax)
        ax.legend()
        return ax
