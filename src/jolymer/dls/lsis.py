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
                                color=m.color, label=m.label, marker=m.marker,
                                linestyle='', **kwargs)
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
        if 'rapp' in kwargs:
            if kwargs['rapp'] is True:
                ax.set_ylabel('$R_H A(R_H)$')
                ax.set_xlabel('$R_H\\mathrm{\\,[nm]}$')
        ax.legend()
        ax.legend(fontsize='x-small')
        return ax

    def rayleighs(self, **kwargs):
        ax, kwargs = self.make_plot(**kwargs)
        for m in self.ms:
            ax = qplot.rayley(m, m.m_buffer, m.m_toluene,
                              ax=ax, color=m.color,
                              marker=m.marker, label=m.label,
                              **kwargs)

        ax.legend()
        qplot.qlabel(ax, 'I')
        ax.set_ylabel('$(I_{sample} - I_{H2O})/ I_{toluene}$')
        # plu.loglog(ax)
        return ax

    def guiniers(self, **kwargs):
        ax, kwargs = self.make_plot(**kwargs)
        for m in self.ms:
            ax = qplot.rayley(m, m.m_buffer, m.m_toluene, ax=ax, color=m.color,
                               marker=m.marker, label=m.label, useq2=True)
        ax.legend()
        qplot.qqlabel(ax, 'ln($I$)')
        ax.set_yscale('log')
        # ax.set_ylabel('$(I_{sample} - I_{H2O})/ I_{toluene}$')
        return ax

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
                 'err_rh': [],
                 'slope': []}
        rapp_fit_dic = {}
        for m in self.ms:
            ax, (xdata, ydata) = qplot.Rapp(m, ax=ax, color=m.color,
                                            fit=self.fit, marker=m.marker,
                                            label=m.label, **kwargs)
            popt, pcov = qplot.plot_qqfit(xdata, ydata, color=m.color,
                                          marker=None, ax=ax,
                                          bounds=[[-1e6, 0], [1e6, 1e3]])
            Rslope = popt[0]
            rh = popt[1]
            pstd = np.sqrt(np.diag(pcov))
            err_rh = pstd[1]
            print(m.label, ' Rh= ', rh, ' pm ', err_rh)
            rdict['label'].append(m.label)
            rdict['rh'].append(rh)
            rdict['err_rh'].append(err_rh)
            rdict['slope'].append(Rslope)
            rapp_fit_dic[f'{m.label}_qq'] = [0, max(xdata)]
            rapp_fit_dic[f'{m.label}_rapp_fit'] = [rh, rh+Rslope*max(xdata)]

        print(rapp_fit_dic)
        ax.legend()
        ax.legend(fontsize='x-small')
        rdf = pd.DataFrame(rdict)
        rapp_fit_df = pd.DataFrame(rapp_fit_dic)
        self.rdf = rdf
        self.rapp_fit_df = rapp_fit_df
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
        qqlist = [self.ms[0].qq(phi) for qhi in self.ms[0].angles]
        theqq = np.array(qqlist)
        for m in self.ms:
            df, dfs = m.get_average_g2(phi, **kwargs)
            datadic[f'{m.label}_tau'] = df.t
            # Omit label, because Origin will use the rowname as label
            datadic[f'{m.label}'] = df.g2
            datadic[f'{m.label}_err_g2'] = df.err_g2
            try:
                dfres = m.get_res(phi)
                fitdic[f'{m.label}_tau'] = dfres.t
                fitdic[f'{m.label}_fit'] = dfres.fit
                dfdist = m.get_Arl(phi)
                thedistt = dfdist.t
                distdic[f'{m.label}_tau'] = dfdist.t
                distdic[f'{m.label}_R'] = dfdist.rapp * 1e9
                distdic[f'{m.label}_dist'] = dfdist.dist
                dfrapp = self.fit.get_phidlstable(m)
                theqq = df.rapp.qq * 1e-12
                rappdict[f'{m.label}_qq'] = dfrapp.qq * 1e-12
                rappdict[f'{m.label}_rapp'] = dfrapp.Rapp * 10**9
                rappdict[f'{m.label}_Gamma'] = dfrapp.Gamma
                rappdict[f'{m.label}_Dapp'] = dfrapp.Dapp * 1e6
            except Exception as e:
                print(e)

        df_data = pd.DataFrame(datadic)
        df_color = pd.DataFrame(colordic)
        for lol in fitdic:
            print(len(fitdic[lol]))
        try:
            df_fit = pd.DataFrame(fitdic)
            df_dist = pd.DataFrame(distdic)
            df_rapp = pd.DataFrame(rappdict)
            _, df_rh = self.Rapps()
            df_rapp_fit = self.rapp_fit_df
        except:
            pass
        print(outpath)
        with pd.ExcelWriter(outpath) as writer:
            df_data.to_excel(writer, sheet_name='data')
            df_color.to_excel(writer, sheet_name='color')
            try:
                df_fit.to_excel(writer, sheet_name='fit')
                df_dist.to_excel(writer, sheet_name='dist')
                df_rapp.to_excel(writer, sheet_name='rapp')
                print('test')
                df_rapp_fit.to_excel(writer, sheet_name='rapp_fit')
                df_rh.to_excel(writer, sheet_name='rh')
            except:
                pass

    def dist_to_excel(self, angle=90, **kwargs):
        pass
