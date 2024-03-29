# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 16:26:57 2020

@author: xcill
"""
import pandas as pd
import numpy as np
import datetime as dt
import os
from os.path  import join
import matplotlib.pyplot as plt
from scipy import optimize, constants, signal

from ..Sample import Sample
from .. import plot_utility as plu
# from .DLS_Measurement import DL
from .. import os_utility as osu
from ..Measurement import Measurement

from matplotlib.backends.backend_pdf import PdfPages

# plt.style.use('seaborn')

# Here come some functions used to generate labels:


def seqlabel(m, fit, seq_number):
    return str(seq_number)


def philabel(m, fit, seq_number):
    phi = m.phifromseq(seq_number)
    return f'{phi} $^\\circ$'


def _dat_compilation(m, seq_numbers, cm, labelfunc, ax=None, **kwargs):
    citer = plu.cm_for_l(cm, seq_numbers)
    if 'title' in kwargs:
        title = kwargs.pop('title')
    for s, color in zip(seq_numbers, citer):
        label = labelfunc(m, None, s)
        ax = m.plot_data(s, ax=ax, color=color, label=label, **kwargs)
    ax.legend()
    ax.set_title(title)
    return ax


def seq_dat(m, seq_numbers, **kwargs):
    return _dat_compilation(m,  seq_numbers, 'viridis', seqlabel, **kwargs)


def phi_dat(m, seq_numbers, kwargs):
    return _dat_compilation(m,  seq_numbers, 'viridis', philabel, **kwargs)


def _fit_compilation(m, seq_numbers, cm, labelfunc, title, ax=None, fit=None,
                     showlegend=True, legendargs={}, **kwargs):
    citer = plu.cm_for_l(cm, seq_numbers)
    for s, color in zip(seq_numbers, citer):
        label = labelfunc(m, fit, s)
        ax = m.plot_fit(s, fit, ax=ax, color=color, label=label, **kwargs)
    if showlegend:
        ax.legend(**legendargs)
    ax.set_title(title)
    return ax

def plot_res(m, phi, ax=None, fit=None, rapp=False, color='black', **kwargs):
    if ax is None:
        fig, ax = plt.subplots()
    if fit is None:
        df, dfs = m.get_average_g2(phi)
    elif fit.name == 'repes':
        df, dfs = m.get_average_g2(phi)
        dfres = m.get_res(phi)
        xdata = df.t
        if rapp:
            xdata = df.rapp
    else:
        df = dfres = fit.get_phidlsfit(m, phi)
    ax.plot(dfres.t, dfres.res, color=color, **kwargs)
    ax.set_ylabel('Relative residuals')
    ax.set_ylim(-0.05, 0.05)
    return ax

def plot_phidls(m, phi, fitcolor='black', ax=None, showres=True, fit=None,
                rapp=False, **kwargs):
    if ax is None:
        fig, ax = plt.subplots()
    if fit is None:
        df, dfs = m.get_average_g2(phi)
        ax.errorbar(df.t, df.g2, df.err_g2, **kwargs)
        showres = False
    elif fit.name == 'repes':
        df, dfs = m.get_average_g2(phi)
        dfres = m.get_res(phi)
        xdata = df.t
        if rapp:
            xdata = df.rapp
        ax.errorbar(xdata, df.g2, df.err_g2, **kwargs)
        ax.errorbar(dfres.t, dfres.fit, color=fitcolor, linestyle='-', marker='')
    else:
        df = dfres = fit.get_phidlsfit(m, phi)
        ax.errorbar(df.t, df.g2, df.err_g2, **kwargs)
        ax.errorbar(df.t, df.fit, color=fitcolor)
    ax.set_xscale('log')
    rescolor = fitcolor
    if 'color' in kwargs:
        rescolor = kwargs['color']
    if showres:
        axres = ax.twinx()
        axres.plot(dfres.t, dfres.res, alpha=0.2, color=rescolor)
        axres.set_ylabel('Relative Residuals')
        axres.set_ylim(-0.05, 0.05)
    return ax


def plot_phidls_dist(m, phi, ax=None, fit='repes', max_is_one=False,
                     rapp=False, A='A', **kwargs):
    if ax is None:
        fig, ax = plt.subplots()
    # if fit == 'contin':
    #     ax.errorbar(df.t, df.g2, df.err_g2, **kwargs)
    #     showres = False
    elif fit.name == 'repes':
        df = m.get_Arl(phi, A=A)
        if rapp:
            xdata = df.rapp * 1e9
        else:
            xdata = df.t
        if max_is_one:
            ydata = df.dist / np.max(df.dist)
        else:
            ydata = df.dist
        ax.plot(xdata, ydata, **kwargs)
    else:
        df = fit.get_phidlsfit(m, phi)
        ax.errorbar(df.t, df.g2, df.err_g2, **kwargs)
    ax.set_xscale('log')
    return ax

def res_compilation(m, phis, cm, ax=None, fit=None,
                    showlegend=True, legendargs={}, **kwargs):
    citer = plu.cm_for_l(cm, phis)
    for phi, color in zip(phis, citer):
        label = f'{phi} °'
        ax = plot_res(m, phi, fit=fit, ax=ax,
                      color=color, label=label, **kwargs)
    if showlegend:
        ax.legend(**legendargs)
    return ax

def phidls_compilation(m, phis, cm, ax=None, fit=None,
                       showlegend=True, legendargs={}, **kwargs):
    citer = plu.cm_for_l(cm, phis)
    for phi, color in zip(phis, citer):
        label = f'{phi} °'
        ax = plot_phidls(m, phi, fit=fit, ax=ax,
                         color=color, label=label, **kwargs)
    if showlegend:
        ax.legend(**legendargs)
    # ax.set_title(title)
    return ax


def phidls_dist_compilation(m, phis, cm, ax=None, fit=None,
                            showlegend=True, legendargs={},
                            **kwargs):
    citer = plu.cm_for_l(cm, phis)
    for phi, color in zip(phis, citer):
        label = f'{phi} °'
        ax = plot_phidls_dist(m, phi, fit=fit, ax=ax,
                              color=color, label=label, **kwargs)
    if showlegend:
        ax.legend(**legendargs)
    # ax.set_title(title)
    return ax


def seq_compilation(m, seq_numbers, fit=None, **kwargs):
    print(fit)
    if 'labelfunc' in kwargs:
        labelfunc = kwargs.pop('labelfunc')
    else:
        labelfunc = seqlabel
    return _fit_compilation(m, seq_numbers, 'viridis', labelfunc, None, fit=fit,
                            **kwargs)


def phi_compilation(m, seq_numbers, fit=None):
    return _fit_compilation(m, seq_numbers, 'viridis', philabel, None, fit=fit)


def _fit_compilations(listofargs, **kwargs):
    num_plots = len(listofargs)
    fig, axes = plu.n_subplots(num_plots)

    for row in axes:
        row[0].set_ylabel('$g_2-1$')
    for ax in axes[-1]:
        ax.set_xlabel('$\\tau$ [s]')
    axes = axes.flatten()[0:num_plots]
    for args, ax in zip(listofargs, axes):
        ax = _fit_compilation(*args, ax=ax, **kwargs)
    return fig, axes


def seq_compilations(m, seqll, fit=None, cm='viridis', title=None, **kwargs):
    listofargs = []
    for seq_numbers in seqll:
        args = [m, seq_numbers, cm, seqlabel, title]
        listofargs.append(args)
    fig, axes = _fit_compilations(listofargs, fit=fit, **kwargs)
    return fig, axes


def _rawdata_page(m, phi, fit=None, per_plot=5, **kwargs):
    seqll = []
    seq_numbers = m.phirange(phi)
    if len(seq_numbers) > 10:
        nax = 4
    else:
        nax = 2
    per_plot = int(np.ceil(len(seq_numbers)/nax))
    for i in range(nax):
        seql = []
        for n in range(i*per_plot, (i+1)*per_plot):
            try:
                seql.append(seq_numbers[n])
            except:
                pass
        seqll.append(seql)

    fig, axes = seq_compilations(m, seqll, fit=fit, title=None)
    fig.suptitle(f'$2\\Theta = $ {phi} $^\\circ$')
    plt.tight_layout()


def rawdata_pdf(m, fit=None, filename=None):
    # filename = os.path.join(m.path, 'rawdata.pdf')
    osu.create_path(m.figures_path)
    # fitname = 'nofit' if fit is None else fit.name
    if filename is None:
        filename = join(Measurement.figures_path, 'temp.pdf')
    else:
        filename = join(m.figures_path, filename)
    with PdfPages(filename) as pdf:
        for phi in m.angles:
            print(phi)
            _rawdata_page(m, phi, fit=fit)
            pdf.savefig()
            plt.close()
    return filename


def _dist_compilation(m, seq_numbers, cm, labelfunc, title, ax=None, fit=None,
                      xspace='t', showlegend=True, legendargs={}, **kwargs):
    citer = plu.cm_for_l(cm, seq_numbers)
    for s, color in zip(seq_numbers, citer):
        label = labelfunc(m, fit, s)
        ax = m.plot_dist(s, fit, ax=ax, color=color,
                         xspace=xspace, label=label, **kwargs)
    if showlegend:
        ax.legend(**legendargs)
    if 'xlabel' in kwargs:
        ax.set_xlabel(kwargs['xlabel'])
    ax.set_title(title)
    return ax


def dist_seq(m, seq_numbers, fit=None, xspace='t', ax=None,
             **kwargs):
    if ax is None:
        fig, (axf, axd) = plt.subplots(nrows=2, figsize=(10, 10))
    else:
        axd=ax
    axd = _dist_compilation(m, seq_numbers, 'viridis', seqlabel, None, fit=fit,
                            ax=axd, xspace=xspace)
    # axf = _fit_compilation(m, seq_numbers, 'viridis', seqlabel, None, fit=fit,
    #                        ax=axf)
    return axd


def dist_phi(m, seq_numbers, xspace='t', figsize=(10, 10), xlim=None, fit=None,
             data_marker='.'):
    fig, (axf, axd) = plt.subplots(nrows=2, figsize=figsize)
    axd = _dist_compilation(m, seq_numbers, 'viridis', philabel, None, fit=fit,
                            ax=axd, xspace=xspace)
    axd.set_ylabel('$\\tau A(\\tau)$ [s]')
    axd.set_xlabel('$\\tau$ [s]')
    axd.set_xlim(xlim)
    axf = _fit_compilation(m, seq_numbers, 'viridis', philabel, None, fit=fit,
                           ax=axf, marker=data_marker)
    axf.set_ylabel('$g_2 - 1$')
    axf.set_xlim(xlim)
    return fig, (axf, axd)


def rh_seq(*args, **kwargs):
    return dist_seq(*args, xspace='rh', **kwargs)


def rh_phi(*args, **kwargs):
    return dist_phi(*args, xspace='rh', **kwargs)


def _raw_contin_page(m, phi, fit=None, per_plot=5, **kwargs):
    seqll = []
    seq_numbers = m.phirange(phi)
    for i in range(4):
        seql = []
        for n in range(i*per_plot, (i+1)*per_plot):
            try:
                seql.append(seq_numbers[n])
            except:
                pass
        seqll.append(seql)

    fig, ((axf1, axf2), (axd1, axd2)) = plt.subplots(nrows=2, ncols=2,
                                                     figsize=(12, 8))
    axd1 = _dist_compilation(m, m.phirange(phi)[0:5], 'viridis',
                             seqlabel, None, fit=fit, ax=axd1)
    axd1.set_xlabel('$\\tau $ [s]')
    axd1.set_ylabel('$\\tau \\cdot G(\\tau) $ [s]')
    axd2 = _dist_compilation(m, m.phirange(phi)[5::], 'viridis', seqlabel,
                             None, fit=fit, ax=axd2)
    axd2.set_xlabel('$\\tau $ [s]')
    axf1 = _fit_compilation(m, m.phirange(phi)[0:5], 'viridis', seqlabel,
                            None, fit=fit, ax=axf1)
    axf1.set_ylabel('$g_2 - 1$')
    axf2 = _fit_compilation(m, m.phirange(phi)[5::], 'viridis', seqlabel,
                            None, fit=fit, ax=axf2)

    fig.suptitle(f'$2\\Theta = $ {phi} $^\\circ$ \t $T = $ {m.TC} $^\\circ$C')
    plt.tight_layout()


def contin_pdf(m, fit=None, filename=None):
    if filename is None:
        filename = join(Measurement.figures_path, 'temp.pdf')
    osu.create_path(m.figures_path)
    with PdfPages(filename) as pdf:
        for phi in m.angles:
            print(phi)
            _raw_contin_page(m, phi, fit=fit)
            pdf.savefig()
            plt.close()


def get_full_df(m, fit, par):
    df = fit.get_phitable(m)
    if par == 'Gamma':
        df['Gamma'] = df.qq * df.Dapp
        constant_fit = [0, df.loc[df.phi == 'constant'].Dapp * max(df.qq)]
        y_intercept = float(df.loc[df.phi == 'y_intercept'][par])
        slope = float(df.loc[df.phi == 'slope'][par])
        linear_fit = [0, max(df.qq) * (y_intercept + max(df.qq)*slope)]
    else:
        constant_fit = [df.loc[df.phi == 'constant'][par] for x in [1, 1]]
        y_intercept = float(df.loc[df.phi == 'y_intercept'][par])
        slope = float(df.loc[df.phi == 'slope'][par])
        linear_fit = [y_intercept, y_intercept + max(df.qq)*slope]
    return df, constant_fit, linear_fit


def qplot(m, fit, par, ax=None, plot_constant=False, plot_linear=False,
          **kwargs):
    df, constant_fit, linear_fit = get_full_df(m, fit, par)
    dfp = df[df.qq > 0]
    if ax is None:
        fig, ax = plt.subplots()
    if 'xlabel' in kwargs:
        ax.set_xlabel(kwargs.pop('xlabel'))
    if 'ylabel' in kwargs:
        ax.set_ylabel(kwargs.pop('ylabel'))
    if 'title' in kwargs:
        ax.set_title(kwargs.pop('title'))
    ax.errorbar(dfp.qq, dfp[par], dfp[f'err_{par}'], fmt='o', **kwargs)
    if plot_constant:
        ax.errorbar([0, max(dfp.qq)], constant_fit, fmt='-', **kwargs)
    if plot_linear:
        # print(linear_fit)
        ax.errorbar([0, max(dfp.qq)], linear_fit, fmt='-', **kwargs)
    return df, ax


def _qplot_mfitlist(mfitlist, parameter, figure=None, **kwargs):
    num_plots = len(mfitlist)
    if figure is None:
        fig, axes = plu.n_subplots(num_plots)
    else:
        fig, axes = figure
    axlist = axes.flatten()
    for fitm, ax in zip(mfitlist, axlist):
        m, fit = fitm
        ax = m.qplot(fit, [parameter], ax=ax, **kwargs)

# def qplot_dapp(mlist, )


def plot_Dapp(m, fit, **kwargs):
    pass


def compare_prob_rej(self, c, meas_num, method='one5'):
    ""
    print("Care!! This only works for one5")
    uloz_list = []
    max_list = []
    phi = "TODO"
    prob_rej = [['B', '.999'], ['C', '.99'], ['D', '.9'],
                ['A', '.5'], ['E', '.1'], ['F', '.01'],
                ['G', '.001'], ['H', '.0001']]
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))
    ax, ax2 = axes
    for i, color in zip(prob_rej, iter(plt.cm.viridis(np.linspace(0, .9, 9)))):
        df = pd.read_csv(self.Arl_file(c, meas_num, one='_one5', A=i[0]),
                         header=None, names=['logtime', 'msdist'])
        df['t'] = 0.001 * 10**df.logtime
        df['dist'] = 0.001 * df.msdist
        ax.plot(df.t, df.dist, color=color, label=i[1])
        uloz, peak_split = self.moa_reader(self.moA_file(c,
                                                         meas_num,
                                                         A=i[0],
                                                         one='_one5'))
        peaks, _ = signal.find_peaks(df.dist, height=.000005)
        maximum = df.loc[peaks]
        maximum.reset_index(drop=True, inplace=True)

        for j in uloz:
            ax.plot([0.001 * j[1]], [-0.000_001], '*', color=color)
            if j[1] > 5 and j[1] < 10:
                uloz_list.append(j[1])

        for j in peak_split:
            ax.plot([0.001 * j[1]], [-0.000_003], "p", color=color)
            if j[1] > .005 and j[1] < 10:
                uloz_list.append(j[1])

        ax.plot(maximum.t, maximum.dist, 'X', color=color)
        for j in maximum.t:
            if j > .005 and j < 0.010:

                max_list.append(1000 * j)

    ax.set_xscale('log')
    ax.legend()
    ax.set_title("Measurement " + str(meas_num) + ' Angle: ' + str(phi)+'°')
    ax.set_xlabel("$\\tau$ [s]")
    ax.set_ylabel("$\\tau_Dw(\\tau_D)$ [s]")

    ax2.plot(['.999',  '.99', '.9', '.5', '.1', '.01', '.001', '.0001'],
             uloz_list, label='$\\left< \\tau_D \\right>$ [ms]')
    ax2.plot(['.999',  '.99', '.9', '.5', '.1', '.01', '.001', '.0001'],
             max_list, label='maximum [ms]')
    ax2.legend()
    ax2.set_xlabel('ProbRej.')
    ax2.set_ylabel('$\\tau_D$ [ms]')
    plt.tight_layout()


def gamma_star_plot(self, method, figure=None, sup_titeling = 'name', labeling = 'method', colormap = plt.cm.ocean(np.linspace(0,.9, 10)),
                    fitname='linear', min_fit=0, max_fi=0):
    max_fit = len(self.angles) if max_fi==0 else max_fi
    if figure == None:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (10, 10))
    else:
        fig, ax = figure
    df_csv = pd.DataFrame()
    qrh_min, qrh_max = [100,0]
    fit_data = pd.read_csv(self.dzrhch_file(method, fitname, min_fit, max_fit))
    for index, c, color in zip(range(self.amc), self.c_string, iter(colormap)):
        df = pd.read_csv(self.avg_par_file(c, method))
        gamma_star =  self.visc * df.decay_rate / (kB * self.T * df.qq**1.5) * 10**-18
        qrh = np.sqrt(df.qq) * fit_data.Rh[index] * 0.001
        label = str(self.c_float[index]) + ' g/l'
        if labeling == 'tag':
            label += self.tag
        elif labeling == 'name_tag':
            label += self.name +self.tag
        ax.plot(qrh, gamma_star, 'o', color = color, label = label)
        for j in qrh:
            if j < qrh_min:
                qrh_min=j
            if j > qrh_max:
                qrh_max=j
        df_csv[c + 'qrh'] = qrh
        df_csv[c + 'gamma_star'] = gamma_star
    df_csv.to_csv(self.dir_dat_method(method) + self.name + '_gamma_star_' + method + '_' + fitname +'_fit_' + str(min_fit+1) +'to'+ str(max_fit) +'.csv' )
    qrh_range = np.arange(qrh_min, qrh_max, (qrh_max - qrh_min)/100)
    ax.plot(qrh_range, (6*np.pi*qrh_range)**-1, label = 'hard sphere')
    ax.set_ylabel("$\\Gamma^*$")
    ax.set_xlabel("$qR_H$")
    ax.legend()
    title = self.name
    if sup_titeling=='method':
        title+= '_' + method
    elif sup_titeling == 'name_tag_method':
        title += '_' + self.tag + '_' + method
    fig.suptitle(title)
    fig.tight_layout()
    return fig, ax
