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
from . import sas_plotlib as splu
from .beaucage import beaucage

from matplotlib.backends.backend_pdf import PdfPages


def binArray(data, axis, binstep, binsize, func=np.nanmean):
    data = np.array(data)
    dims = np.array(data.shape)
    argdims = np.arange(data.ndim)
    argdims[0], argdims[axis] = argdims[axis], argdims[0]
    data = data.transpose(argdims)
    data = [func(np.take(data, np.arange(int(i*binstep),
            int(i*binstep+binsize)), 0), 0)
            for i in np.arange(dims[axis]//binstep)]
    data = np.array(data).transpose(argdims)
    return data


class Ms:

    def __init__(self, ms):
        self.ms = ms
        self.model = ms[0].model

    def __index__(self, index):
        return self.ms[index]

    def __len__(self):
        return len(self.ms)

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

    def fits(self, fits=True, shiftby=1, dataqmin=0, dataqmax=2300, **kwargs):
        ax, kwargs = self.make_plot(**kwargs)
        shift = shiftby**len(self.ms)
        for m in self.ms:
            # print(m.label)
            dfall = m.get_data(cout=False)
            df = dfall[m.dataqmin: m.dataqmax]
            ax.errorbar(df.q, df.I*shift, df.err_I*shift,
                        marker=m.marker, color=m.color,
                        linestyle='', label=m.label,
                        elinewidth=0.2, **kwargs)
            df = dfall[m.iqmin: m.iqmax]
            if fits:
                m.fit_dict, m.fit_df = m.model.fit(m, bounds=m.bounds,
                                                   iqmax=m.iqmax, p0=m.p0,
                                                   iqmin=m.iqmin,
                                                   fixed_parameters=m.fixed_pars)
                ax.errorbar(m.fit_df.q, m.fit_df.fit*shift, marker='',
                            color='black')
                m.partext = m.model.get_text(m.fit_dict)
            shift = shift / shiftby
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend(fontsize='xx-small')
        ax.set_xlabel('$q \\,\\mathrm{[nm^{-1}]}$')
        ax.set_ylabel('$I \\,\\mathrm{[cm^{-1}]}$')

    def data(self, shiftby=1, **kwargs):
        return self.fits(fits=False, shiftby=shiftby, **kwargs)

    def res(self, **kwargs):
        if 'ax' in kwargs:
            ax = kwargs.pop('ax')
        else:
            fig, ax = plt.subplots(nrows=1, ncols=1,
                                   figsize = (5, 4),
                                   sharex=True, sharey=True,
                                   squeeze=True)
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
        ax.legend(fontsize='xx-small')
        ax.set_ylabel('Normalized Residuals')
        ax.set_xlabel('$q\\mathrm{\\,[nm^{-1}]}$')

    def kratky(self, **kwargs):
        ax, kwargs = self.make_plot(**kwargs)
        for m in self.ms:
            df = m.get_data(cout=False)[m.iqmin:m.iqmax]
            xdata = df.q
            ydata = 1000 * df.I * df.q**2
            if 'binnum' in kwargs:
                binnum = kwargs.pop('binnum')
                xdata = binArray(df.q, 0, binnum, binnum)
                ydata = 1000 * binArray(df.I, 0, binnum, binnum) * xdata**2
            label = m.label
            ax.errorbar(xdata, ydata, marker=m.marker, color=m.color,
                        linestyle='', label=label, elinewidth=0.2, **kwargs)
            # ax.legend()
            ax.legend(fontsize='xx-small')
            #             ax.annotate(m.partext, xy=(0.0, 0.0), xycoords='axes fraction')
            # ax.grid()
            ax.set_ylim(0,5)
            ax.set_xlim(0, 2.5)
            ax.set_xlabel('$q \\,\\mathrm{[nm^{-1}]}$')
            #         axes[0][1].set_xlabel('$q$ [1/nm]')
            ax.set_ylabel('$I\\cdot q^2 \\mathrm{\\,[0.001nm^{-2}cm^{-1}]}$')

    def q3I(self, **kwargs):
        ax, kwargs = self.make_plot(**kwargs)
        for m in self.ms:
            df = m.get_data(cout=False)[m.iqmin:m.iqmax]
            label = m.label
            ax.errorbar(df.q, df.q**3 * df.I * 10_000,  marker=m.marker, color=m.color,
                        linestyle='', label=label, elinewidth=0.2, **kwargs)
            #             ax.annotate(m.partext, xy=(0.0, 0.0), xycoords='axes fraction')
            # ax.grid()
            ax.set_ylim(0,0.0025)
            ax.set_xlim(0, 0.6)
            ax.set_xlabel('$q \\,\\mathrm{[nm^{-1}]}$')
            #         axes[0][1].set_xlabel('$q$ [1/nm]')
            ax.set_ylabel('$I\\cdot q^3 \\mathrm{\\,[nm^{-3}cm^{-1}]}$')
            ax.set_ylabel('$I\\cdot q^3 \\mathrm{\\,[A.U.]}$')
        ax.legend(fontsize='xx-small')

    def debye(self, **kwargs):
        pass

    def get_results(self):
        par_dict = {}
        for par in self.model.parameters:
            par_dict[par] = []
            par_dict[f'std_{par}'] = []
            for m in self.ms:
                par_dict[par].append(m.fit_dict[par])
                par_dict[f'std_{par}'].append(m.fit_dict['std_'+par])
        self.df = pd.DataFrame(par_dict)
        return self.df

    def save_fit_results(self, path, **kwargs):
        self.df = self.get_results()
        self.df.to_csv(path)
        return self.df

    def plot_par(self, par, **kwargs):
        ax, kwargs = self.make_plot(**kwargs)
        df = self.get_results()
        xdata = [m.label for m in self.ms]
        ydata = df[par]
        err_ydata = df[f'std_{par}']
        ax.errorbar(xdata, ydata, err_ydata, **kwargs)
        return ax

    def markdown_table(self, fixed_pars=[], e_pars=[], **kwargs):
        out = '| parameter |'
        for m in self.ms:
            out += f'{m.label} |'
        out += '\n | --- |'
        for m in self.ms:
            out += f'--- |'
        out += '\n'
        for par in self.model.parameters:
            out += f'| {par} |'
            for m in self.ms:
                print(m.fit_dict[f'std_{par}'])
                if m.fit_dict[f'std_{par}'] == 'fixed':
                    out += '{:.2e} fix |'.format(m.fit_dict[par])
                elif par in e_pars:
                    print(m.fit_dict[f'std_{par}'])
                    out += '{:.2e} |'.format(m.fit_dict[par], m.fit_dict['std_'+par])
                else:
                    out += '{:.2f} ± {:.2f} |'.format(m.fit_dict[par], m.fit_dict['std_'+par])
            out += '\n'
        if 'beaucage_scale' in self.model.parameters:
            out += f'| beaucage_C |'
            for m in self.ms:
                out += '{:.2e} |'.format(beaucage.get_C(m.fit_dict['beaucage_scale'],
                                                        m.fit_dict['beaucage_rg'],
                                                        m.fit_dict['beaucage_exp']))
            out += '\n'
        return out

    def markdown_sasview(self):
        pass

    def markdown_latex_table(self, **kwargs):
        out = '| parameter |'
        for m in self.ms:
            out += f'{m.label} |'
        out += '\n | --- |'
        for m in self.ms:
            out += f' --- |'
        out += '\n'
        for par in self.model.parameters:
            out += f'| ${self.model.pardict["par"]["latex"]}$ |'
            for m in self.ms:
                out += '{:.2f} \\pm {:.2f} |'.format(m.fit_dict[par], m.fit_dict['std_'+par])
            out += '\n'
        return out

    def to_excel(self, fits=True, onlyraw=False,
                 outpath='~/test_saxs_to_excel.xlsx',
                 shiftby=10, **kwargs):
        datadic = {}
        fitdic = {}
        colordic = {'labels': [m.label for m in self.ms],
                    'color': [matplotlib.colors.to_hex(
                        m.color) for m in self.ms]}
        shift = shiftby ** len(self.ms)
        for m in self.ms:
            dfall = m.get_data(cout=False)
            df = dfall[m.dataqmin: m.dataqmax]
            datadic[f'{m.label}_q'] = df.q
            datadic[f'{m.label}'] = df.I
            datadic[f'{m.label}_err_I'] = df.err_I
            m.fit_dict, m.fit_df = m.model.fit(
                    m, bounds=m.bounds,
                    iqmax=m.iqmax, p0=m.p0,
                    iqmin=m.iqmin,
                    fixed_parameters=m.fixed_pars)
            fitdic[f'{m.label}_q'] = df.q
            fitdic[f'{m.label}'] = df.I * shift
            fitdic[f'{m.label}_err_I'] = df.err_I * shift
            fitdic[f'{m.label}_q'] = m.fit_df.q * shift
            fitdic[f'{m.label}_fit'] = m.fit_df.fit * shift
            shift = shift / shiftby
        df_data = pd.DataFrame(datadic)
        df_color = pd.DataFrame(colordic)
        df_fit = pd.DataFrame(fitdic)
        df_par = self.get_results()
        with pd.ExcelWriter(outpath) as writer:
            df_data.to_excel(writer, sheet_name='data')
            df_color.to_excel(writer, sheet_name='color')
            df_fit.to_excel(writer, sheet_name='fit')
            df_par.to_excel(writer, sheet_name='par')
