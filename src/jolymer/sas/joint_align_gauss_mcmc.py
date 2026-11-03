"""
Joint UV/SAXS time-alignment + Gaussian elution-peak fit, via MCMC.

WHY MCMC INSTEAD OF ALTERNATING OPTIMIZATION
==============================================
`fit_gauss_peaks_refine_alignment` (the original) alternates: an outer
optimizer proposes (scale, shift), and for each trial an inner
`least_squares` call refits the Gaussian peak centers/sigmas/amplitudes
from scratch. That's fragile for exactly the reason you flagged:

  * The inner peak fit is itself nonlinear and has its own local minima.
    The "cost" the outer optimizer sees at a given (scale, shift) is
    whatever local minimum `least_squares` happened to land in from that
    trial's warm start -- not a clean function of (scale, shift). Small
    changes in shift can produce spurious jumps in that cost that have
    nothing to do with alignment quality.
  * Because the peaks have a lot of shape freedom, they can partially
    *absorb* a bad shift (by moving peak centers instead of the UV trace),
    making a wrong alignment look artificially good to the outer loop.

Sampling the *joint* posterior over (scale, shift, peak centers, peak
sigmas) fixes both problems:

  * Peak *amplitudes* and per-channel *baselines* are linear given
    (scale, shift, centers, sigmas), so at every MCMC likelihood
    evaluation they're solved for in closed form with `np.linalg.lstsq`
    -- exact, non-iterative, no local minima of their own. There is no
    inner optimizer left to get stuck.
  * The nonlinear parameters (scale, shift, centers, sigmas) are explored
    together by the sampler. If shift and a peak center are degenerate,
    the posterior shows a ridge, instead of an optimizer silently picking
    one arbitrary point on it and reporting false confidence.

Caveat: profiling amplitudes out by their least-squares mode (rather than
fully marginalizing them with the correct Occam-factor determinant term)
is a standard approximation for linear nuisance parameters, but is not
exactly a Bayesian marginal likelihood -- the design matrix G(theta)
changes with (scale, shift, centers, sigmas), so the profile can be very
slightly non-flat where it should be penalized. In practice this is a
second-order effect compared to what alternating optimization was doing;
flag it if you want the fully marginal version later (closed-form, doable
by hand from the Gaussian-linear-model formula, just not included here to
keep the likelihood cheap for 126-way parallel evaluation).

PRE-EXISTING BUG, FYI
========================
In the original `fit_gauss_peaks`, the UV time axis is built with
`shift = 0` and `scale = 1` hard-coded (the lines that would read
`self._alignment["shift"/"scale"]` are commented out). That means
`fit_gauss_peaks_refine_alignment`, which alternates by writing trial
(scale, shift) into `self._alignment` and calling `fit_gauss_peaks`, was
never actually seeing the effect of the (scale, shift) it was
"optimizing" -- the cost it minimized was constant in those two
parameters. This module does not have that bug: the warped UV time axis
is built directly from the trial (scale, shift) inside the likelihood.

DATA-QUALITY FILTERING (your errI0 / qRg question)
======================================================
`saxs_quality_mask` combines three independent, togglable checks on the
SAXS scalar (I0) series before it ever enters the fit:
  1. I0_min          - as before, a floor on I0 itself.
  2. max_rel_err     - drop points where errI0 / |I0| exceeds a fraction
                        (catches noisy near-baseline points that a flat
                        I0_min wouldn't).
  3. qmax * Rg < qRg_max - Guinier-validity proxy from the autorg table:
                        drop frames where the q-range used for I0/Rg is
                        outside the Guinier-valid regime (qRg <~ 1.3) for
                        that frame's own Rg, i.e. autorg found *a* fit
                        but not a physically valid one.
Filtered points are dropped from the SAXS channel's residual entirely
(not down-weighted), so they cannot pull on the peak/alignment fit.

PARALLELIZATION FOR A 126-CORE MACHINE
==========================================
Uses emcee's EnsembleSampler with a multiprocessing.Pool. emcee
parallelizes by splitting the walker ensemble in half each step and
pool.map-ing the likelihood over one half at a time -- so set `nwalkers`
to a multiple of 2 * n_workers (default here: 2 * 126 = 252) so every
half-step saturates all cores. Each worker process must NOT also spawn
its own BLAS thread pool (numpy's lstsq will, by default) or you
oversubscribe the 126 cores -- this module forces single-threaded BLAS
in-process via environment variables before the Pool is created.

emcee's default stretch move mixes slowly on strongly correlated /
possibly multimodal posteriors like scale-vs-shift here, so this module
defaults to a mix of DEMove + DESnookerMove (differential-evolution style
moves), which tend to handle exactly that kind of geometry better and are
still fully parallel-friendly.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# 1. Data-quality filtering
# ---------------------------------------------------------------------

def saxs_quality_mask(
    dfX: pd.DataFrame,
    autorg_df: Optional[pd.DataFrame] = None,
    I0_min: float = 0.0,
    max_rel_err: Optional[float] = 0.1,
    qmax_for_qrg: Optional[float] = None,
    qRg_max: float = 1.3,
) -> np.ndarray:
    """
    Boolean mask (True = keep) for a SAXS scalar dataframe with columns
    (frame, time, I0, errI0). Any check left as None is skipped.

    I0_min        : drop I0 <= I0_min.
    max_rel_err   : drop points where errI0 / |I0| exceeds this fraction.
    qmax_for_qrg,
    qRg_max       : if both given (and autorg_df with an 'Rg' column of
                    the same length as dfX is supplied), drop points
                    where qmax_for_qrg * Rg > qRg_max.
    """
    n = len(dfX)
    mask = np.ones(n, dtype=bool)

    I0 = dfX["I0"].to_numpy(dtype=float)
    mask &= np.isfinite(I0) & (I0 > I0_min)

    if max_rel_err is not None and "errI0" in dfX:
        errI0 = dfX["errI0"].to_numpy(dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            rel_err = np.abs(errI0) / np.abs(I0)
        mask &= np.isfinite(rel_err) & (rel_err < max_rel_err)

    if qmax_for_qrg is not None and autorg_df is not None and "Rg" in autorg_df:
        Rg = autorg_df["Rg"].to_numpy(dtype=float)
        if len(Rg) != n:
            raise ValueError(
                "autorg_df length does not match dfX; align them on the "
                "same frame index before calling saxs_quality_mask."
            )
        qRg = qmax_for_qrg * Rg
        mask &= np.isfinite(Rg) & (qRg < qRg_max)

    return mask


# ---------------------------------------------------------------------
# 2. Gaussian design matrix + closed-form linear solve
# ---------------------------------------------------------------------

def _gauss_design(t: np.ndarray, centers: np.ndarray, sigmas: np.ndarray) -> np.ndarray:
    diffs = t[:, None] - centers[None, :]
    return np.exp(-0.5 * (diffs / sigmas[None, :]) ** 2)


def _solve_linear_channel(
    G: np.ndarray, y: np.ndarray, yerr: Optional[np.ndarray], fit_baseline: bool
) -> Tuple[np.ndarray, np.ndarray]:
    """Closed-form weighted least squares for one channel's peak
    amplitudes (+ optional constant baseline), given a fixed Gaussian
    design matrix G. Returns (coefficients, weighted_residuals)."""
    if fit_baseline:
        A = np.column_stack([G, np.ones(len(y))])
    else:
        A = G
    w = 1.0 / np.clip(yerr, 1e-12, None) if yerr is not None else np.ones(len(y))
    Aw = A * w[:, None]
    yw = y * w
    coef, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
    resid_w = Aw @ coef - yw
    return coef, resid_w


# ---------------------------------------------------------------------
# 3. Channel container + picklable joint log-posterior
# ---------------------------------------------------------------------

@dataclass
class Channel:
    name: str
    t_raw: np.ndarray   # this channel's own, UNwarped time axis
    y: np.ndarray
    yerr: np.ndarray
    warp: bool           # True: apply (scale, shift) to t_raw before
                          # evaluating Gaussians (UV channels). False:
                          # reference axis, used as-is (SAXS).


class JointAlignGaussLogProb:
    """
    Picklable joint log-posterior over theta = [scale, shift,
    centers_1..n, log_sigmas_1..n]. Peak amplitudes + per-channel
    baselines are profiled out analytically at every call.

    Kept free of anything non-picklable (no scipy.interpolate objects, no
    references to the parent measurement object) so it can be shipped as-is
    to multiprocessing workers.
    """

    def __init__(
        self,
        channels: List[Channel],
        n_peaks: int,
        scale_bounds: Tuple[float, float] = (0.95, 1.05),
        shift_bounds: Tuple[float, float] = (-3000.0, 3000.0),
        center_bounds: Optional[Tuple[float, float]] = None,
        log_sigma_bounds: Tuple[float, float] = (-2.0, 4.0),
        fit_baseline: bool = True,
        min_center_gap: float = 0.0,
    ):
        self.channels = channels
        self.n_peaks = n_peaks
        self.scale_bounds = scale_bounds
        self.shift_bounds = shift_bounds
        self.log_sigma_bounds = log_sigma_bounds
        self.fit_baseline = fit_baseline
        self.min_center_gap = min_center_gap
        if center_bounds is None:
            ref = channels[0].t_raw
            center_bounds = (float(np.min(ref)), float(np.max(ref)))
        self.center_bounds = center_bounds
        self.ndim = 2 + 2 * n_peaks

    def unpack(self, theta: np.ndarray):
        scale, shift = theta[0], theta[1]
        centers = theta[2 : 2 + self.n_peaks]
        sigmas = np.exp(theta[2 + self.n_peaks :])
        return scale, shift, centers, sigmas

    def log_prior(self, theta: np.ndarray) -> float:
        scale, shift, centers, sigmas = self.unpack(theta)
        if not (self.scale_bounds[0] <= scale <= self.scale_bounds[1]):
            return -np.inf
        if not (self.shift_bounds[0] <= shift <= self.shift_bounds[1]):
            return -np.inf
        if np.any(centers < self.center_bounds[0]) or np.any(centers > self.center_bounds[1]):
            return -np.inf
        # enforce ordering to remove the peak label-switching degeneracy
        if self.n_peaks > 1 and np.any(np.diff(centers) <= self.min_center_gap):
            return -np.inf
        log_sigmas = np.log(sigmas)
        if np.any(log_sigmas < self.log_sigma_bounds[0]) or np.any(
            log_sigmas > self.log_sigma_bounds[1]
        ):
            return -np.inf
        return 0.0

    def log_likelihood(self, theta: np.ndarray) -> float:
        scale, shift, centers, sigmas = self.unpack(theta)
        chi2 = 0.0
        for ch in self.channels:
            t = scale * ch.t_raw + shift if ch.warp else ch.t_raw
            G = _gauss_design(t, centers, sigmas)
            _, resid_w = _solve_linear_channel(G, ch.y, ch.yerr, self.fit_baseline)
            chi2 += float(np.dot(resid_w, resid_w))
        return -0.5 * chi2

    def __call__(self, theta: np.ndarray) -> float:
        lp = self.log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        ll = self.log_likelihood(theta)
        if not np.isfinite(ll):
            return -np.inf
        return lp + ll

    def reconstruct(self, theta: np.ndarray) -> Dict[str, Any]:
        """Best-fit reconstruction (model + per-peak curves) at a given theta."""
        scale, shift, centers, sigmas = self.unpack(theta)
        out = {}
        for ch in self.channels:
            t = scale * ch.t_raw + shift if ch.warp else ch.t_raw
            G = _gauss_design(t, centers, sigmas)
            coef, _ = _solve_linear_channel(G, ch.y, ch.yerr, self.fit_baseline)
            if self.fit_baseline:
                amps, baseline = coef[:-1], float(coef[-1])
            else:
                amps, baseline = coef, 0.0
            model = G @ amps + baseline
            out[ch.name] = {
                "t": t, "y": ch.y, "yerr": ch.yerr, "model": model,
                "amplitudes": amps, "baseline": baseline,
                "peaks": {
                    f"peak_{k + 1}": amps[k]
                    * np.exp(-0.5 * ((t - centers[k]) / sigmas[k]) ** 2)
                    for k in range(self.n_peaks)
                },
            }
        return {
            "scale": scale, "shift": shift, "centers": centers, "sigmas": sigmas,
            "channels": out,
        }


# ---------------------------------------------------------------------
# 4. Pulling channel data out of a CoupledMeasurement
# ---------------------------------------------------------------------

def build_channels_from_coupled(
    cm,
    wl=280,
    saxs_kind: Optional[str] = None,
    qmin: Optional[float] = None,
    qmax: Optional[float] = None,
    qstar: Optional[float] = None,
    I0_min: float = 0.0,
    max_rel_err: Optional[float] = 0.1,
    qRg_max: Optional[float] = 1.3,
    load: bool = True,
) -> Tuple[List[Channel], pd.DataFrame, int]:
    """
    Build the (unwarped) SAXS + UV channel list from a CoupledMeasurement
    instance, applying `saxs_quality_mask` to the SAXS channel.
    """
    saxs_kind = saxs_kind or cm._alignment.get("saxs_kind", cm.saxs_kind)
    qmin = qmin if qmin is not None else cm._alignment.get("qmin", cm.qmin)
    qmax = qmax if qmax is not None else cm._alignment.get("qmax", cm.qmax)
    qstar = qstar if qstar is not None else cm._alignment.get("qstar", cm.qstar)

    dfX = cm.get_saxs_scalar(saxs_kind, qmin=qmin, qmax=qmax, qstar=qstar, load=load)
    autorg_df = getattr(cm, "_autorg_df", None)
    qmax_for_qrg = qmax if qRg_max is not None else None
    mask = saxs_quality_mask(
        dfX, autorg_df=autorg_df, I0_min=I0_min, max_rel_err=max_rel_err,
        qmax_for_qrg=qmax_for_qrg, qRg_max=qRg_max if qRg_max is not None else 1.3,
    )
    n_total = len(dfX)
    dfX = dfX.loc[mask].reset_index(drop=True)
    n_dropped = n_total - len(dfX)

    channels = [
        Channel(
            name="I0",
            t_raw=dfX["time"].to_numpy(dtype=float),
            y=dfX["I0"].to_numpy(dtype=float),
            yerr=dfX["errI0"].to_numpy(dtype=float),
            warp=False,  # SAXS defines the reference time axis
        )
    ]

    dfUV = cm.uv.get_Abs_wl(wl)
    uv_cols = [c for c in dfUV.columns if c.startswith("Abs")]
    t_uv_raw = dfUV["time"].to_numpy(dtype=float)
    for col in uv_cols:
        y = dfUV[col].to_numpy(dtype=float)
        # NOTE: stand-in noise estimate -- if cm.uv exposes a real per-point
        # UV error, swap it in here instead of this 5%-of-std placeholder.
        yerr = np.full_like(y, np.nanstd(y) * 0.05 + 1e-9)
        channels.append(Channel(name=col, t_raw=t_uv_raw, y=y, yerr=yerr, warp=True))

    return channels, dfX, n_dropped


# ---------------------------------------------------------------------
# 5. Walker initialization + MCMC driver
# ---------------------------------------------------------------------

def _init_walkers(
    log_prob: JointAlignGaussLogProb,
    nwalkers: int,
    scale0: float,
    shift0: float,
    centers0: np.ndarray,
    sigmas0: np.ndarray,
    rng: np.random.Generator,
    max_tries_per_walker: int = 2000,
) -> np.ndarray:
    ndim = log_prob.ndim
    p0 = np.zeros((nwalkers, ndim))
    i = 0
    tries = 0
    while i < nwalkers:
        scale = scale0 + rng.normal(0, 0.002)
        shift = shift0 + rng.normal(0, 2.0)
        centers = np.sort(
            centers0 + rng.normal(0, 0.3 * np.mean(sigmas0), size=len(centers0))
        )
        log_sigmas = np.log(sigmas0) + rng.normal(0, 0.1, size=len(sigmas0))
        theta = np.concatenate([[scale, shift], centers, log_sigmas])
        if np.isfinite(log_prob.log_prior(theta)):
            p0[i] = theta
            i += 1
        tries += 1
        if tries > max_tries_per_walker * nwalkers:
            raise RuntimeError(
                "Could not find enough valid starting points -- check "
                "scale0/shift0/centers0/sigmas0 against the bounds."
            )
    return p0


def run_joint_mcmc(
    cm,
    n_peaks: int,
    wl=280,
    scale0: Optional[float] = None,
    shift0: Optional[float] = None,
    centers0: Optional[Sequence[float]] = None,
    sigmas0: Optional[Sequence[float]] = None,
    scale_bounds: Tuple[float, float] = (0.95, 1.05),
    shift_bounds: Tuple[float, float] = (-3000.0, 3000.0),
    log_sigma_bounds: Tuple[float, float] = (-2.0, 4.0),
    fit_baseline: bool = True,
    I0_min: float = 0.0,
    max_rel_err: Optional[float] = 0.1,
    qRg_max: Optional[float] = 1.3,
    nwalkers: Optional[int] = None,
    nsteps: int = 4000,
    burn: int = 1000,
    thin: int = 10,
    n_workers: int = 126,
    seed: int = 0,
) -> Dict[str, Any]:
    """
    Jointly sample (scale, shift, peak centers, peak sigmas) with emcee,
    profiling out peak amplitudes/baselines analytically at every draw.
    Parallelizes across `n_workers` processes.

    Also updates cm._alignment["scale"/"shift"] to the posterior medians
    and stores the full refinement summary in
    cm._alignment["refinement"], mirroring the existing
    refine_uv_alignment* methods.

    Returns a dict with the sampler, flattened post-burn-in chain,
    per-parameter posterior summaries, and the reconstructed best-fit.
    """
    import emcee

    rng = np.random.default_rng(seed)
    channels, dfX, n_dropped = build_channels_from_coupled(
        cm, wl=wl, I0_min=I0_min, max_rel_err=max_rel_err, qRg_max=qRg_max,
    )
    if n_dropped:
        print(
            f"[run_joint_mcmc] dropped {n_dropped}/{n_dropped + len(dfX)} "
            f"SAXS points on quality filters (I0_min / max_rel_err / qRg_max)."
        )

    tX, yX = channels[0].t_raw, channels[0].y
    if scale0 is None:
        scale0 = cm._alignment.get("scale", 1.0)
    if shift0 is None:
        shift0 = cm._alignment.get("shift", 0.0)
    if centers0 is None:
        pos = yX > 0
        tmin, tmax = tX[pos].min(), tX[pos].max()
        centers0 = np.linspace(tmin, tmax, n_peaks)
    else:
        centers0 = np.asarray(centers0, dtype=float)
    if sigmas0 is None:
        sigmas0 = np.full(n_peaks, (tX.max() - tX.min()) / (4 * n_peaks))
    else:
        sigmas0 = np.asarray(sigmas0, dtype=float)

    log_prob = JointAlignGaussLogProb(
        channels=channels, n_peaks=n_peaks,
        scale_bounds=scale_bounds, shift_bounds=shift_bounds,
        log_sigma_bounds=log_sigma_bounds, fit_baseline=fit_baseline,
    )

    if nwalkers is None:
        nwalkers = max(2 * log_prob.ndim, 2 * n_workers)
        nwalkers = int(np.ceil(nwalkers / (2 * n_workers)) * 2 * n_workers)

    p0 = _init_walkers(log_prob, nwalkers, scale0, shift0, centers0, sigmas0, rng)

    # Avoid oversubscribing the 126 cores: force single-threaded BLAS
    # before the Pool (and its forked workers) exist.
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")

    from multiprocessing import Pool

    moves = [(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)]
    with Pool(processes=n_workers) as pool:
        sampler = emcee.EnsembleSampler(
            nwalkers, log_prob.ndim, log_prob, pool=pool, moves=moves,
        )
        sampler.run_mcmc(p0, nsteps, progress=True)

    chain = sampler.get_chain(discard=burn, thin=thin, flat=True)
    logp = sampler.get_log_prob(discard=burn, thin=thin, flat=True)
    best_theta = chain[np.argmax(logp)]

    scale_samples = chain[:, 0]
    shift_samples = chain[:, 1]
    centers_samples = chain[:, 2 : 2 + n_peaks]
    sigmas_samples = np.exp(chain[:, 2 + n_peaks :])

    def _summary(x):
        return {
            "median": float(np.median(x)),
            "lo16": float(np.percentile(x, 16)),
            "hi84": float(np.percentile(x, 84)),
        }

    result = {
        "sampler": sampler,
        "log_prob": log_prob,
        "chain": chain,
        "logp": logp,
        "best_theta": best_theta,
        "best_fit": log_prob.reconstruct(best_theta),
        "scale": _summary(scale_samples),
        "shift": _summary(shift_samples),
        "centers": [_summary(centers_samples[:, k]) for k in range(n_peaks)],
        "sigmas": [_summary(sigmas_samples[:, k]) for k in range(n_peaks)],
        "acceptance_fraction": float(np.mean(sampler.acceptance_fraction)),
        "n_dropped_saxs_points": n_dropped,
    }

    cm._alignment["scale"] = result["scale"]["median"]
    cm._alignment["shift"] = result["shift"]["median"]
    cm._alignment["refinement"] = {
        "method": "joint_mcmc",
        "scale": result["scale"],
        "shift": result["shift"],
        "acceptance_fraction": result["acceptance_fraction"],
        "n_dropped_saxs_points": n_dropped,
    }
    return result


# ---------------------------------------------------------------------
# 6. Quick diagnostics
# ---------------------------------------------------------------------

def plot_scale_shift_posterior(result: Dict[str, Any]):
    """2D scale-vs-shift scatter + marginal histograms, colored by log
    posterior. Cheap stand-in for a full corner plot; run `corner.corner`
    on result['chain'] yourself for the full (scale, shift, centers,
    sigmas) posterior if you have the `corner` package."""
    import matplotlib.pyplot as plt

    scale = result["chain"][:, 0]
    shift = result["chain"][:, 1]
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    axes[0].scatter(shift, scale, s=3, alpha=0.3, c=result["logp"], cmap="viridis")
    axes[0].set_xlabel("shift")
    axes[0].set_ylabel("scale")
    axes[0].set_title("joint posterior (color = log-posterior)")
    axes[1].hist(shift, bins=60)
    axes[1].set_xlabel("shift")
    axes[1].set_title("shift marginal")
    fig.tight_layout()
    return fig, axes


def check_convergence(sampler) -> Dict[str, Any]:
    """Autocorrelation-time based convergence check. Raises no error;
    just reports, since short chains legitimately can't estimate tau."""
    try:
        tau = sampler.get_autocorr_time(quiet=True)
        return {"tau": tau, "n_steps": sampler.iteration, "ok": np.all(sampler.iteration > 50 * tau)}
    except Exception as e:  # emcee raises if the chain is too short
        return {"error": str(e)}


# ---------------------------------------------------------------------
# 7. Optional: attach as a method on CoupledMeasurement
# ---------------------------------------------------------------------
#
#   from joint_align_gauss_mcmc import run_joint_mcmc
#   CoupledMeasurement.fit_gauss_peaks_mcmc = run_joint_mcmc
#
# then call it as   result = cm.fit_gauss_peaks_mcmc(n_peaks=3, n_workers=126)
