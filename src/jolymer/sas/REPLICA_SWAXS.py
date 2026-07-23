from scipy import optimize
import numpy as np
import re
from os.path import join
from matplotlib.colors import LogNorm, Normalize
from matplotlib.cm import ScalarMappable
from pylab import cm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd

import shutil
import os

import sasmodels
from sasmodels import data as sasmodels_data
import pyFAI
import fabio
from pathlib import Path

import MDAnalysis as mda

from .. import os_utility as osu
from .GROMACS_SWAXS import GROMACS_SWAXS
from .. import jocolors

from dataclasses import dataclass, field

def _colorbar(mappable):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib.pyplot as plt
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax)
    plt.sca(last_axes)
    return cbar

def _make_grid():
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(7, 4))
    gs = GridSpec(4, 3, figure=fig)
    ax_large = fig.add_subplot(gs[0:4, 1:3])
    ax_small1 = fig.add_subplot(gs[0, 0])
    ax_small2 = fig.add_subplot(gs[1, 0])
    ax_small3 = fig.add_subplot(gs[2, 0])
    ax_small4 = fig.add_subplot(gs[3, 0])
    axes = np.array([ax_small1, ax_small2, ax_small3, ax_small4])
    return fig, axes, ax_large

def _wsum(q, I_list, weights):
    """
    Weighted sum of spectra.

    I_list: list/array of I_i(q)
    weights: array of w_i
    """
    I = np.array(I_list)
    w = np.array(weights)

    return np.sum(w[:, None] * I, axis=0) / np.sum(w)

@dataclass
class REPLICA_SWAXS(GROMACS_SWAXS):

    gss: list[GROMACS_SWAXS] = field(default_factory=list)
    name: str = ""

    def __post_init__(self):
        for i, gs in enumerate(gss):
            gs.NAME = f'R{i}_{gs.NAME}'

    def __len__(self):
        return len(self.gss)

    def __iter__(self):
        return iter(self.gss)

    @property
    def nreplicas(self):
        return len(self.gss)

    def get_chi2(self):
        rows = []
        for irep, gs in enumerate(self.gss):
            df = gs.get_all()
            if "chi2" not in df.columns:
                continue
            tmp = df[["time", "chi2"]].copy()
            tmp["replica"] = irep
            rows.append(tmp)
        if not rows:
            return pd.DataFrame()
        return pd.concat(rows, ignore_index=True)

    def best(self, n=10):
        chi2df = self.get_chi2()
        return (
            chi2df
            .sort_values("chi2")
            .head(n)
            .reset_index(drop=True)
        )

    def plot_replicas(
        self,
        figsize=(7,5),
        every_n=1,
        annotate=True,
        **kwargs
        ):
        n = len(self.gss)
        ncols = int(np.ceil(np.sqrt(n)))
        nrows = int(np.ceil(n / ncols))
        fig, axes = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            figsize=figsize
        )
        axes = np.atleast_1d(axes).flatten()
        outdicts = []
        for i, (ax, gs) in enumerate(zip(axes, self.gss)):
            outdict = gs.plot_spectra(
                ax=ax,
                every_n=every_n,
                **kwargs
            )
            if annotate:
                ax.annotate(
                    f"{chr(97+i)})",
                    xy=(-0.1, 1.0),
                    xycoords="axes fraction"
                )
            outdicts.append(outdict)
        return fig, axes, outdicts

    def collect_spectra(self):
        dfs = []
        for irep, gs in enumerate(self.gss):
            outdict = gs.plot_spectra(plot=False)
            for idf, df in enumerate(outdict["df"]):
                if idf < gs.min_index:
                    continue
                tmp = df.copy()
                tmp["replica"] = irep
                dfs.append(tmp)
        return dfs

    def select_best_spectra(self, n=50):
        selected = []
        best = self.best(n)
        for _, row in best.iterrows():
            gs = self.gss[int(row.replica)]
            df = gs.get_all()
            nearest = np.argmin(
                np.abs(df.time - row.time)
            )
            selected.append(
                gs.load_spectrum(nearest)
            )
        return selected

    def get_rg_dataframe(self, qmax=0.08):
        rows = []
        for irep, gs in enumerate(self.gss):
            outdict = gs.plot_spectra(plot=False)
            # print('outdict', outdict)
            for ispec, (df, chi2) in enumerate(
                    zip(outdict["df"], outdict["chi2"])):
                if ispec < gs.min_index:
                    continue
                try:
                    rgdict = SAXS_Measurement.get_rg(self, df=df, qmax=qmax)
                    rows.append({
                        "replica": irep,
                        "spectrum": ispec,
                        "Rg": rgdict['Rg'],
                        "chi2": chi2,
                        "df": df
                    })
                except Exception as e:
                    print(e)
                    print('get_rg failed')
                    print('df', df)
        return pd.DataFrame(rows)

    def plot_rg_colored_spectra(
        self,
        ax=None,
        qmax_rg=0.08,
        highlight="best_chi2",
        cmap="plasma"):
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure
        rgdf = self.get_rg_dataframe(qmax=qmax_rg)
        rgdf = rgdf[np.isfinite(rgdf.Rg)]
        norm = Normalize(
            rgdf.Rg.min(),
            rgdf.Rg.max()
        )
        cmap_obj = plt.get_cmap(cmap)
        for _, row in rgdf.sort_values("Rg").iterrows():
            ax.plot(
                row.df.q,
                row.df.I,
                color=cmap_obj(norm(row.Rg)),
                alpha=0.5,
                lw=1,
            )
        if highlight == "best_chi2":
            best = rgdf.loc[rgdf.chi2.idxmin()]
            ax.plot(
                best.df.q,
                best.df.I,
                color="green",
                lw=3,
                label=(
                    fr"$R_g={best.Rg:.2f}\,\AA$"
                    "\n"
                    fr"$\chi^2={best.chi2:.2f}$"
                )
            )
        ax.set_xscale("log")
        ax.set_yscale("log")
        sm = ScalarMappable(
            norm=norm,
            cmap=cmap_obj
        )
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax)
        cbar.set_label(r"$R_g$ ($\AA$)")
        self.gss[0].plot_data(color='black', marker='', ax=ax, unit='A')
        return rgdf

    def plot_rg_histogram(
            self,
            ax=None,
            qmax_rg=0.08,
            bins=30,
            weighted=False):
        if ax is None:
            fig, ax = plt.subplots()
        rgdf = self.get_rg_dataframe(qmax=qmax_rg)
        rg = rgdf.Rg.values
        chi2 = rgdf.chi2.values
        mask = np.isfinite(rg)
        rg = rg[mask]
        chi2 = chi2[mask]
        if weighted:
            weights = np.exp(-chi2)
        else:
            weights = None
        ax.hist(
            rg,
            bins=bins,
            weights=weights
        )
        ax.set_xlabel(r"$R_g$ ($\AA$)")
        ax.set_ylabel("Counts")
        return rgdf

    def ensemble_spectrum(
            self,
            chi2_scale=5.0,
            min_q=None,
            max_q=None,
            weight_mode="exp"):
        """
        Build weighted ensemble SAXS curve from all replicas.
        """
        all_I = []
        all_weights = []
        q_ref = None
        for gs in self.gss:
            out = gs.get_all()
            for _, row in out.iterrows():
                df = row["df"]
                chi2 = row.get("chi2", np.nan)
                if np.isnan(chi2):
                    continue
                # define weights
                if weight_mode == "exp":
                    w = np.exp(-(chi2 - out["chi2"].min()) / chi2_scale)
                elif weight_mode == "inv":
                    w = 1.0 / (chi2 + 1e-6)
                else:
                    w = 1.0
                q = df["q"].values
                I = df["I"].values
                if q_ref is None:
                    q_ref = q
                else:
                    # safety check: enforce same grid
                    if len(q) != len(q_ref):
                        continue
                all_I.append(I)
                all_weights.append(w)
        all_I = np.array(all_I)
        all_weights = np.array(all_weights)
        weights = all_weights[:, None]
        I_ens = np.sum(weights * all_I, axis=0) / np.sum(weights)
        return q_ref, I_ens

    def plot_best_fit(self, exp_df, ax=None, chi2_scale=5.0):
        """
        Compare weighted ensemble to experimental SAXS curve.
        """
        if ax is None:
            fig, ax = plt.subplots()
        q_sim, I_sim = self.ensemble_spectrum(
            chi2_scale=chi2_scale,
            weight_mode="exp"
        )
        ax.plot(
            exp_df["q"],
            exp_df["I"],
            "k-",
            label="experiment"
        )
        ax.plot(
            q_sim,
            I_sim,
            "r-",
            lw=2,
            label="weighted ensemble"
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$q$")
        ax.set_ylabel(r"$I(q)$")
        ax.legend()
        return ax

    def get_dfs(self):
        dfs = []
        chi2s = []
        for gs in self.gss:
            outdict = gs.plot_spectra( plot=False,
                              every_n=1)
            for odf in outdict['df']:
                dfs.append(odf)
            for ochi2 in outdict['chi2']:
                chi2s.append(ochi2)
        return dfs, np.array(chi2s)

    def return_wsfit(self):
        from scipy import optimize

        df_exp = self.gss[0].get_data()
        df_exp.q = df_exp.q/10
        dfs, _ = self.get_dfs()
        df_exp = df_exp[df_exp.q>=dfs[0].q.min()]
        df_exp = df_exp[df_exp.q<dfs[0].q.max()]
        def combomodel(q, *weights):
            outI = np.zeros_like(q, dtype=float)
            for weight, df in zip(weights, dfs):
                # outI += weight * df.I.values[0:326]
                outI += weight * df.I.values[0:len(df_exp)]
            return outI
        n = len(dfs)
        p0 = np.ones(n) / n   # initial guess: equal weights
        popt, pcov = optimize.curve_fit(
            combomodel,
            df_exp.q, df_exp.I, sigma=df_exp.err_I,
            p0=p0, bounds=(0,1)
        )
        fig, ax = plt.subplots()
        gs = self.gss[0]
        gs.plot_data(ax=ax, unit='A', linestyle='', marker='o', color=jocolors.tstum1)
        # popt_equal = 1/len(dfs) * np.ones(len(dfs))
        popt_equal = popt
        df_fit = pd.DataFrame({'q': df_exp.q, 'fit': combomodel(df_exp.q, *popt_equal)})
        chi2 = gs.get_chi2(df_exp, df_fit)[1]
        label=f'fit; $\\chi^2={chi2:.2f}$'
        ax.errorbar(df_exp.q, combomodel(df_exp.q, *popt_equal),
                    linestyle='-', marker='', color=jocolors.tstum7, label=label)
        ax.legend()
        return popt

    def run_analysis(self, analysis_name, *args, **kwargs):
        for gs in self.gss:
            results, aux_results = gs.run_analysis(analysis_name)
            gs.analysis_results[analysis_name] = [results, aux_results]

    def load_analysis(self, analysis_name):
        for gs in self.gss:
            results, aux_results = gs.load_analysis(analysis_name)
            gs.analysis_results[analysis_name] = [results, aux_results]

