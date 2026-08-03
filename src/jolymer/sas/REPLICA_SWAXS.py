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
from .SAXS_Measurement import SAXS_Measurement
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
    subfolders: bool=True

    def __post_init__(self):
        print('replica from:')
        self.NAME = f'R_{self.gss[0].NAME}'
        for irep, gs in enumerate(self.gss):
            if self.subfolders:
                gs.npt_filename = f'../{gs.npt_filename}'
            else:
                gs.npt_filename = f'{gs.npt_filename}'
            gs.irep = irep
            gs.NAME = f'R{irep}_{gs.NAME}'
            print(gs.NAME)

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
            outdict = gs.plot_spectra(plot=False,
                                      get_Rg=True)
            # print('outdict', outdict)
            for ispec, (df, chi2, Rg, err_Rg) in enumerate(
                    zip(outdict["df"], outdict["chi2"],
                        outdict["Rg"], outdict["err_Rg"])):
                if ispec < gs.min_index:
                    continue
                rgdict = SAXS_Measurement.get_rg(self, df=df, qmax=qmax)
                try:
                    rows.append({
                        "replica": irep,
                        "spectrum": ispec,
                        "Rg": rgdict['Rg'],
                        "err_Rg": rgdict['err_Rg'],
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
            fig, ax = plt.subplots(figsize=(4,3))
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

    def plot_average( self,
                    qmax_rg=0.08,
                     ax=None, fit_color=jocolors.tstum7):
        if ax is None:
            fig, ax = plt.subplots(figsize=(4,3))
        else:
            fig = ax.figure
        self.gss[0].plot_data(color=jocolors.tstum1, marker='o',
                              linestyle='',
                              ax=ax, unit='A')
        rgdf = self.get_rg_dataframe(qmax=qmax_rg)
        dfs = []
        for _, row in rgdf.sort_values("Rg").iterrows():
            dfs.append(row.df)
        average_df = dfs[0].copy()
        n = min(len(df) for df in dfs)
        average_df['I'] = np.vstack([df.I.iloc[:n] for df in dfs]).mean(axis=0)
        ax.errorbar(average_df.q, average_df.I, marker='', linestyle='-', color=jocolors.tstum7)
        return ax


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
        dfs = []
        for irep, gs in enumerate(self.gss):
            results, aux_results = gs.run_analysis(
                analysis_name, *args, **kwargs
            )
            gs.analysis_results[analysis_name] = [results, aux_results]
            results = results.copy()
            results["irep"] = irep
            dfs.append(results)
        return pd.concat(dfs, ignore_index=True)

    def load_analysis(self, analysis_name):
        dfs = []
        for irep, gs in enumerate(self.gss):
            results, aux_results = gs.load_analysis(analysis_name)
            gs.analysis_results[analysis_name] = [results, aux_results]
            results = results.copy()
            results["irep"] = irep
            dfs.append(results)
        return pd.concat(dfs, ignore_index=True)

    def load_analysis(self, analysis_name):
        for irep, gs in enumerate(self.gss):
            results, aux_results = gs.load_analysis(analysis_name)
            gs.analysis_results[analysis_name] = [results, aux_results]

    def get_occupancy_histogramm(self, df, parameters=['Rg'], bins=10):
        X = df[parameters].to_numpy()
        H, edges = np.histogramdd(X, bins=bins)
        bin_indices = []
        for i, edge in enumerate(edges):
            idx = np.digitize(X[:, i], edge) - 1
            idx = np.clip(idx, 0, len(edge)-2)
            bin_indices.append(idx)
            print(idx)
        name = "_".join(parameters)
        df = df.copy()
        df[f"bin_{name}"] = list(zip(*bin_indices))
        return df, H, edges

    def get_bins_outdir(self):
        if self.subfolders:
            return Path(self.gss[0].mdpath).parent / self.NAME
        else:
            return Path(self.gss[0].mdpath) / self.NAME

    def save_bins(self, df, parameters, save=True):
        bin_column = '_'.join(parameters)
        bin_column = f'bin_{bin_column}'
        def load_u(gs):
            gs.to_sol()
            return gs.get_u()
        universes = [load_u(gs) for gs in self.gss]
        outdir = self.get_bins_outdir()
        osu.create_path(outdir)
        lookup = []
        for bin in sorted(df[bin_column].unique()):
            bin_name = "_".join(str(int(x)) for x in bin)
            outname = outdir / f"{bin_column}_{bin_name}.xtc"
            with mda.Writer(str(outname), n_atoms=universes[0].atoms.n_atoms) as W:
                subset = df[df[bin_column] == bin]
                if len(subset) < 100:
                    continue
                else:
                    print('bin', outname, 'N', len(subset))
                lookup_row = {
                    "file": outname,
                    "N": len(subset),
                }
                for p in parameters:
                    lookup_row[p] = subset[p].mean()
                lookup.append(lookup_row)
                if save:
                    for _, row in subset.iterrows():
                        u = universes[int(row["irep"])]
                        u.trajectory[int(row["frame"])]
                        W.write(u.atoms)
        lookup_df = pd.DataFrame(lookup)
        lookup_df.to_csv(
            outdir / "lookup.dat",
            sep="\t",
            index=False
        )
        return lookup_df

    def get_bins_final(self, parameters):
        bin_column = '_'.join(parameters)
        bin_column = f'bin_{bin_column}'
        gs = self.gss[0]
        outdir = self.get_bins_outdir()
        osu.create_path(outdir)
        outdict = {'dfs': [],
                   'bins': [],
                   'par_values': []}
        lookup_df = pd.read_csv(
                outdir / "lookup.dat",
                sep="\t")
        lookup_df['bin'] = (
            lookup_df['file']
            .str.split(f'{parameters[-1]}_').str[1]
            .str.split('.xtc').str[0]
            .str.split('_')
        )
        for _, row in lookup_df.iterrows():
            bin_name = row["file"].split(f"{parameters[-1]}_")[1].split(".xtc")[0]
            dfd = gs.get_data()
            dfg = gs.get_gromacs(
                path=outdir,
                filename=f"rerun_{bin_column}_{bin_name}_final.xvg"
            )
            dfgs = gs.scale_and_offset_fit(dfd, dfg)["df"]
            par_values = row[parameters].tolist()
            outdict["dfs"].append(dfgs)
            outdict["bins"].append(bin_name)
            outdict["par_values"].append(par_values)
        return outdict

    def plot_bins_final(self, parameters,
                        ranges={},
                        ax=None):
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        pardfdict = self.get_bins_final(parameters)
        gs = self.gss[0]
        print(parameters.index('Rg'))
        for df, bin, par_values in zip(
                pardfdict['dfs'],
                pardfdict['bins'],
                pardfdict['par_values']):
            ok = True
            for par, value in zip(parameters, par_values):
                if par in ranges:
                    low, high = ranges[par]
                    if not (low <= value <= high):
                        ok = False
                        break
            if ok:
                label = ''
                for par, parv in zip(parameters, par_values):
                    label += f'{par}={parv:.1f};'
                print(bin)
                print(label)
                df['err_I'] = 1
                rgdict = SAXS_Measurement.get_rg(gs, df=df, qmax=0.1)
                print('RG', rgdict['Rg'])
                ax.errorbar(df.q, df.I, label=label, marker='')

    def plot_bins(self, parameters, H, edges, axes=None):
        from itertools import combinations
        import copy
        pairs = list(combinations(range(len(parameters)), 2))
        if axes is None:
            nrows = len(pairs)
            ncols = 1
            if len(pairs)>3:
                nrows = len(pairs) // 2
                ncols = 2
            fig, axes = plt.subplots(
                nrows=ncols,
                ncols=nrows,
                squeeze=False,
                figsize=(3*nrows, 2.5*ncols))
            axes = axes.flatten()
        else:
            fig = axes[0].get_figure()
        for ax, (i, j) in zip(axes, pairs):
            H2 = H
            sum_axes = sorted(
                set(range(H.ndim)) - {i, j},
                reverse=True
            )
            for axis in sum_axes:
                H2 = H2.sum(axis=axis)
            cmap = copy.copy(plt.cm.viridis)
            cmap.set_bad("white")
            Hplot = np.ma.masked_less(H2, 100)
            im = ax.pcolormesh(
                edges[i],
                edges[j],
                Hplot.T,
                shading="auto",
                cmap=cmap
            )
            fig.colorbar(im, ax=ax, label="Frames")
            ax.set_xlabel(parameters[i])
            ax.set_ylabel(parameters[j])
        return fig

    @staticmethod
    def get_weights(x, log_w0):
        """Softmax parametrization: guarantees w > 0 and sum(w) = 1."""
        z = log_w0 + x
        z = z - z.max()
        w = np.exp(z)
        w /= w.sum()
        return w


    def maxent_reweight(self, i_exp, calc_curves, w0=None, sigma=None, theta=1.0,
                         maxiter=5000, tol=1e-14):
        """
        Parameters
        ----------
        i_exp : (n_q,) array
            Target/experimental curve.
        calc_curves : (n_bins, n_q) array
            Calculated curve per bin, on the same q-grid as i_exp.
        w0 : (n_bins,) array, optional
            Prior weights. Defaults to uniform (1/n_bins each). If your bins
            came from unevenly populated sampling (not flat/umbrella sampling),
            pass the actual population fractions instead of uniform.
        sigma : (n_q,) array, optional
            Experimental uncertainty per q-point. Defaults to ones (unweighted
            least squares in chi2 terms).
        theta : float
            Entropy regularization strength. Larger = closer to prior.
        Returns
        -------
        w_opt : (n_bins,) array, posterior weights
        n_eff : float, effective number of bins = exp(S_rel)
        chi2  : float, reduced chi2 of the reweighted fit
        res   : scipy OptimizeResult
        """
        calc_curves = np.asarray(calc_curves, dtype=float)
        i_exp = np.asarray(i_exp, dtype=float)
        n_bins, n_q = calc_curves.shape
        if w0 is None:
            w0 = np.ones(n_bins) / n_bins
        else:
            w0 = np.asarray(w0, dtype=float)
            w0 = w0 / w0.sum()
        if sigma is None:
            sigma = np.ones(n_q)
        else:
            sigma = np.asarray(sigma, dtype=float)
        log_w0 = np.log(w0)
        def objective(x):
            w = get_weights(x, log_w0)
            calc_avg = w @ calc_curves
            chi2 = np.sum(((i_exp - calc_avg) / sigma) ** 2) / n_q
            # avoid log(0) for numerically-zero weights
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.where(w > 0, w / w0, 1.0)
                s_rel = -np.sum(w * np.log(ratio))
            return chi2 - theta * s_rel
        x0 = np.zeros(n_bins)
        res = minimize(objective, x0, method="L-BFGS-B",
                        options={"maxiter": maxiter, "ftol": tol})
        w_opt = get_weights(res.x, log_w0)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(w_opt > 0, w_opt / w0, 1.0)
            s_rel = -np.sum(w_opt * np.log(ratio))
        n_eff = np.exp(s_rel)
        calc_avg = w_opt @ calc_curves
        chi2 = np.sum(((i_exp - calc_avg) / sigma) ** 2) / n_q
        return w_opt, n_eff, chi2, res


    def theta_scan(i_exp, calc_curves, thetas, w0=None, sigma=None):
        """
        Scan theta and return chi2, S_rel, n_eff for each value, so you can pick
        theta at the 'elbow' of the chi2 vs S_rel L-curve (standard BME practice:
        don't just pick the theta with lowest chi2, that overfits).
        """
        rows = []
        for theta in thetas:
            w_opt, n_eff, chi2, _ = maxent_reweight(
                i_exp, calc_curves, w0=w0, sigma=sigma, theta=theta
            )
            n_bins = calc_curves.shape[0]
            w0_arr = w0 if w0 is not None else np.ones(n_bins) / n_bins
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.where(w_opt > 0, w_opt / w0_arr, 1.0)
                s_rel = -np.sum(w_opt * np.log(ratio))
            rows.append({"theta": theta, "chi2": chi2, "S_rel": s_rel, "n_eff": n_eff})
        return pd.DataFrame(rows)

    def build_reweighting_inputs(self, parameters):
        """
        Turns self.get_bins_final(parameters) + self.gss[0].get_data() into the
        arrays maxent_reweight() needs.
        Assumes each per-bin dataframe (outdict['dfs'][k]) and the target
        dataframe share column names q_col / i_col (and optionally sigma_col).
        If bins are not on the exact same q-grid as the target (they may not be,
        since dfgs comes from a separate GROMACS rerun), each curve is linearly
        interpolated onto the target's q-grid.
        """
        outdict = self.get_bins_final(parameters)
        target_df = self.gss[0].get_data()
        q_target = target_df[q_col].to_numpy()
        i_target = target_df[i_col].to_numpy()
        if sigma_col in target_df.columns:
            sigma = target_df[sigma_col].to_numpy()
        else:
            sigma = np.ones_like(i_target)
        calc_curves = []
        for dfgs in outdict["dfs"]:
            q_bin = dfgs[q_col].to_numpy()
            i_bin = dfgs[i_col].to_numpy()
            i_interp = np.interp(q_target, q_bin, i_bin)
            calc_curves.append(i_interp)
        calc_curves = np.array(calc_curves)
        par_values = np.array(outdict["par_values"], dtype=float)  # (n_bins, n_params)
        bin_names = outdict["bins"]
        return {
            "q_target": q_target,
            "i_target": i_target,
            "sigma": sigma,
            "calc_curves": calc_curves,
            "par_values": par_values,
            "bin_names": bin_names,
            "parameters": parameters,
        }


    def reweight_parameters(self, parameters, theta=1.0, w0=None):
        """
        End-to-end: build inputs, reweight, and return a tidy DataFrame with
        prior vs. posterior weights and per-parameter reweighted averages.
        Example
        -------
        result, w_opt, n_eff = reweight_parameters(self, ['Ree', 'Rg', 'fbs'], theta=1.0)
        print(result['summary'])
        """
        inputs = self.build_reweighting_inputs(self, parameters)
        n_bins = inputs["calc_curves"].shape[0]
        if w0 is None:
            w0 = np.ones(n_bins) / n_bins  # assume evenly sampled bins by default
        w_opt, n_eff, chi2, res = maxent_reweight(
            inputs["i_target"], inputs["calc_curves"],
            w0=w0, sigma=inputs["sigma"], theta=theta,
        )
        per_bin = pd.DataFrame(inputs["par_values"], columns=parameters)
        per_bin["bin"] = inputs["bin_names"]
        per_bin["w0"] = w0
        per_bin["w_opt"] = w_opt
        prior_avg = per_bin[parameters].mul(w0, axis=0).sum()
        post_avg = per_bin[parameters].mul(w_opt, axis=0).sum()
        summary = pd.DataFrame({"prior_avg": prior_avg, "reweighted_avg": post_avg})
        result = {
            "per_bin": per_bin,
            "summary": summary,
            "n_eff": n_eff,
            "chi2": chi2,
            "theta": theta,
        }
        return result, w_opt, n_eff

        def run_maxEnt(self):
            n_bins, n_q = 20, 50
            q = np.linspace(0.01, 0.5, n_q)
            true_w = rng.dirichlet(np.ones(n_bins))
            curves = rng.normal(loc=1.0, scale=0.2, size=(n_bins, n_q)) * np.exp(-5 * q)
            w_opt, n_eff, chi2, res = self.maxent_reweight(i_exp, curves, theta=0.1)
            print("converged:", res.success, "chi2:", chi2, "n_eff:", n_eff)

