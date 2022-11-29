from .. import database_operations as dbo
import pandas as pd
import numpy as np
import datetime as dt
import os
import matplotlib.pyplot as plt
import matplotlib
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
        ax.tick_params(axis="x", bottom=True, top=True, labelbottom=True, labeltop=False)
        ax.tick_params(axis="y", left=True, right=True, labelleft=True, labelright=False)
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
        ax.legend(fontsize='x-small')
        return ax

    def dists(self, phi, **kwargs):
        ax, kwargs = self.make_plot(**kwargs)
        for m in self.ms:
            ax = dp.plot_phidls_dist(m, phi, fit=self.fit, ax=ax,
                                     color=m.color, **kwargs,
                                     label=m.label, marker='')
        ax.set_ylabel('$\\tau_D A(\\tau_D)$')
        ax.set_xlabel('$\\tau_D\\mathrm{\\,[s]}$')
        ax.legend()
        ax.legend(fontsize='x-small')
        return ax

    def rayleighs(self, **kwargs):
        ax, kwargs = self.make_plot(**kwargs)
        for m in self.ms:
            ax = qplot.plot_IminI(m, m.m_buffer, m.m_toluene,
                                  ax=ax, color=m.color,
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
        ax.legend(fontsize='x-small')
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
        ax.legend(fontsize='x-small')
        rdf = pd.DataFrame(rdict)
        self.rdf = rdf
        return ax, rdf

    def Gammas(self, **kwargs):
        ax, kwargs = self.make_plot(**kwargs)
        for m in self.ms:
            ax, (xdata, ydata) = qplot.Gamma(m, ax=ax, color=m.color,
                                             fit=self.fit, marker=m.marker,
                                             label=m.label, **kwargs)
            qplot.plot_qqfit(xdata, ydata, color=m.color, marker=None, ax=ax)
        ax.legend(fontsize='x-small')
        return ax

    def Rh(self, **kwargs):
        ax, kwargs = self.make_plot(**kwargs)
        xdata = self.rdf.label
        ydata = self.rdf.rh
        err_y = self.rdf.err_rh
        ax.errorbar(xdata, ydata, yerr=err_y, **kwargs)
        ax.set_ylabel('$R_H$ [nm]')
        return ax

    def fit_to_excel(self, phi=90, onlyraw=False,
                     outpath='~/test.xlsx', **kwargs):
        datadic = {}
        fitdic = {}
        distdic = {}
        rappdict = {}
        colordic = {'labels': [m.label for m in self.ms],
                    'color': [matplotlib.colors.to_hex(
                        m.color) for m in self.ms]}
        for m in self.ms:
            df, dfs = m.get_average_g2(phi, **kwargs)
            datadic[f'{m.label}_tau'] = df.t
            # Omit label, because Origin will use the rowname as label
            datadic[f'{m.label}'] = df.g2
            datadic[f'{m.label}_err_g2'] = df.g2
            try:
                dfres = m.get_res(phi)
                fitdic[f'{m.label}_tau'] = dfres.t
                fitdic[f'{m.label}_fit'] = dfres.fit
                dfdist = m.get_Arl(phi)
                distdic[f'{m.label}_tau'] = dfdist.t
                distdic[f'{m.label}_R'] = dfdist.rapp * 1e9
                distdic[f'{m.label}_dist'] = dfdist.dist
                dfrapp = self.fit.get_phidlstable(m)
                rappdict[f'{m.label}_qq'] = dfrapp.qq * 1e-12
                rappdict[f'{m.label}_rapp'] = dfrapp.Rapp * 10**9
            except:
                pass
        df_data = pd.DataFrame(datadic)
        df_color = pd.DataFrame(colordic)
        try:
            df_fit = pd.DataFrame(fitdic)
            df_dist = pd.DataFrame(distdic)
            df_rapp = pd.DataFrame(rappdict)
            _, df_rh = self.Rapps()
        except:
            pass
        with pd.ExcelWriter(outpath) as writer:
            df_data.to_excel(writer, sheet_name='data')
            df_color.to_excel(writer, sheet_name='color')
            try:
                df_fit.to_excel(writer, sheet_name='fit')
                df_dist.to_excel(writer, sheet_name='dist')
                df_rapp.to_excel(writer, sheet_name='rapp')
                df_rh.to_excel(writer, sheet_name='rh')
            except:
                pass

    def dist_to_excel(self, angle=90, **kwargs):
        pass
