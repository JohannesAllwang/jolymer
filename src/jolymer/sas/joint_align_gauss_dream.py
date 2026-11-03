"""
Bumps/DREAM backend for the joint UV/SAXS alignment + Gaussian
elution-peak fit.

Same physical model and same reasons for going Bayesian/joint as in
`joint_align_gauss_mcmc.py` (see that module's docstring for the full
rationale on why alternating optimization gives a misleading picture, the
pre-existing shift/scale bug in the original `fit_gauss_peaks`, and the
errI0/qRg quality filtering). This module reuses that module's data
plumbing (`Channel`, `saxs_quality_mask`, `build_channels_from_coupled`,
the closed-form linear amplitude solve) and swaps emcee for `bumps`
+ DREAM as the sampler.

WHY BUMPS/DREAM AS AN ALTERNATIVE TO EMCEE
==============================================
Bumps (Kienzle et al., used throughout the SasView/refl1d/reflectometry
ecosystem you're already adjacent to via ATSAS) is built specifically for
"physical parameters with bounds -> nonlinear least squares or Bayesian
fit" problems, and its DREAM sampler (DiffeRential Evolution Adaptive
Metropolis) was designed for exactly this failure mode: strongly
correlated / possibly multimodal posteriors from a handful of physical
parameters, which is a good description of (scale, shift) vs. peak
shape here. Practically, compared to the emcee version:

  * Parameter bounds ARE the prior -- you declare
    `Parameter(value=..., bounds=(lo, hi))` once and don't hand-roll a
    `log_prior` (the peak-ordering constraint is the one thing still done
    by hand, as a penalty term in `nllf`, same as before).
  * `bumps.fitters.fit(problem, method="dream", parallel=N, ...)` handles
    the process pool for you (`bumps.mapper.MPMapper`), including
    cloudpickling the problem to workers -- more forgiving than emcee's
    plain-pickle Pool if you want to define the Fitness/model in a
    notebook or a closure rather than strictly at module level.
  * You get DREAM's built-in Gelman-Rubin-style convergence check for
    free, plus `bumps.dream.state` tools for trace/corner plots.

Trade-off: bumps' `Parameter` objects carry mutable `.value` state that
`nllf()` reads off; this module still forbids any non-numpy state in the
Fitness object's hot path so cloudpickling stays cheap, but it is a
slightly different style than emcee's purely-functional `theta -> logp`.

Requires: `pip install bumps`

Note: bumps 1.0.5's chain "outlier" removal (`outliers="iqr"`, the
default in the CLI) hits an internal indexing bug on some
pop/thin/chain-length combinations (`IndexError` in
`dream/state.py:_replace_outlier`), reproduced during testing of this
module. `run_joint_dream` defaults to `outliers="none"` to sidestep it;
pass `outliers="iqr"` yourself if you want DREAM's outlier-chain
detection and your bumps version doesn't hit the bug.
"""

from __future__ import annotations

import os
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from jolymer.sas.joint_align_gauss_mcmc import (  # reuse the shared model/data plumbing
    Channel,
    saxs_quality_mask,
    build_channels_from_coupled,
    _gauss_design,
    _solve_linear_channel,
)


# ---------------------------------------------------------------------
# Fitness: bumps' equivalent of JointAlignGaussLogProb
# ---------------------------------------------------------------------

class JointAlignGaussFitness:
    """
    bumps `Fitness` implementation of the joint (scale, shift, peak
    centers, peak log-sigmas) model. Peak amplitudes + per-channel
    baselines are profiled out analytically inside `nllf()`, exactly as
    in the emcee version -- no nested nonlinear optimizer.
    """

    def __init__(
        self,
        channels: List[Channel],
        n_peaks: int,
        centers0: Sequence[float],
        sigmas0: Sequence[float],
        scale0: float = 1.0,
        shift0: float = 0.0,
        scale_bounds: Tuple[float, float] = (0.95, 1.05),
        shift_bounds: Tuple[float, float] = (-3000.0, 3000.0),
        center_bounds: Optional[Tuple[float, float]] = None,
        log_sigma_bounds: Tuple[float, float] = (-2.0, 4.0),
        fit_baseline: bool = True,
        min_center_gap: float = 0.0,
        ordering_penalty: float = 1e8,
    ):
        from bumps.parameter import Parameter

        self.channels = channels
        self.n_peaks = n_peaks
        self.fit_baseline = fit_baseline
        self.min_center_gap = min_center_gap
        self.ordering_penalty = ordering_penalty

        if center_bounds is None:
            ref = channels[0].t_raw
            center_bounds = (float(np.min(ref)), float(np.max(ref)))
        self.center_bounds = center_bounds

        self.scale = Parameter(value=scale0, bounds=scale_bounds, name="scale")
        self.shift = Parameter(value=shift0, bounds=shift_bounds, name="shift")
        self.centers = [
            Parameter(value=float(c), bounds=center_bounds, name=f"center{i}")
            for i, c in enumerate(centers0)
        ]
        self.log_sigmas = [
            Parameter(value=float(np.log(s)), bounds=log_sigma_bounds, name=f"log_sigma{i}")
            for i, s in enumerate(sigmas0)
        ]

    # -- required Fitness protocol methods --------------------------
    def parameters(self):
        return {
            "scale": self.scale,
            "shift": self.shift,
            "centers": self.centers,
            "log_sigmas": self.log_sigmas,
        }

    def update(self):
        pass

    def numpoints(self):
        return sum(len(ch.y) for ch in self.channels)

    def nllf(self) -> float:
        centers = np.array([c.value for c in self.centers])
        sigmas = np.exp(np.array([s.value for s in self.log_sigmas]))
        scale, shift = self.scale.value, self.shift.value
        # ordering constraint: kills the peak label-switching degeneracy,
        # the one thing not expressible purely via per-parameter bounds
        if self.n_peaks > 1 and np.any(np.diff(centers) <= self.min_center_gap):
            return self.ordering_penalty
        chi2 = 0.0
        for ch in self.channels:
            t = scale * ch.t_raw + shift if ch.warp else ch.t_raw
            G = _gauss_design(t, centers, sigmas)
            _, resid_w = _solve_linear_channel(G, ch.y, ch.yerr, self.fit_baseline)
            chi2 += float(np.dot(resid_w, resid_w))
        return 0.5 * chi2

    def residuals(self) -> np.ndarray:
        centers = np.array([c.value for c in self.centers])
        sigmas = np.exp(np.array([s.value for s in self.log_sigmas]))
        scale, shift = self.scale.value, self.shift.value
        out = []
        for ch in self.channels:
            t = scale * ch.t_raw + shift if ch.warp else ch.t_raw
            G = _gauss_design(t, centers, sigmas)
            _, resid_w = _solve_linear_channel(G, ch.y, ch.yerr, self.fit_baseline)
            out.append(resid_w)
        return np.concatenate(out)

    def resynth_data(self):
        pass

    def restore_data(self):
        pass

    def save(self, basename):
        pass

    def plot(self, view="linear"):
        pass

    # -- convenience: best-fit reconstruction at current parameter values --
    def reconstruct(self) -> Dict[str, Any]:
        centers = np.array([c.value for c in self.centers])
        sigmas = np.exp(np.array([s.value for s in self.log_sigmas]))
        scale, shift = self.scale.value, self.shift.value
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
        return {"scale": scale, "shift": shift, "centers": centers, "sigmas": sigmas, "channels": out}


# ---------------------------------------------------------------------
# Driver: build channels from a CoupledMeasurement, run DREAM
# ---------------------------------------------------------------------

def run_joint_dream(
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
    samples: int = 200_000,
    burn: int = 20_000,
    thin: int = 5,
    pop: int = 20,
    n_workers: int = 126,
    seed: int = 0,
    outliers: str = "none",
) -> Dict[str, Any]:
    """
    Jointly fit (scale, shift, peak centers, peak sigmas) with bumps'
    DREAM sampler, profiling out amplitudes/baselines in closed form at
    every evaluation (see JointAlignGaussFitness.nllf). Parallelizes the
    population evaluation across `n_workers` processes via bumps'
    MPMapper.

    `samples` is the total number of posterior draws across the whole
    population (bumps' convention), so with population size
    `pop * nvars` the number of DREAM generations is
    samples / (pop * nvars) -- for 126 workers you generally want a
    population at least ~2-3x n_workers so each generation's proposals
    keep the pool busy; `pop` here is *per free parameter*, matching
    bumps' own `--pop` CLI option.

    Updates cm._alignment["scale"/"shift"] to the posterior medians and
    stores the refinement summary in cm._alignment["refinement"], as the
    other alignment methods do.
    """
    from bumps.fitproblem import FitProblem
    import bumps.fitters as fitters

    rng = np.random.default_rng(seed)
    channels, dfX, n_dropped = build_channels_from_coupled(
        cm, wl=wl, I0_min=I0_min, max_rel_err=max_rel_err, qRg_max=qRg_max,
    )
    if n_dropped:
        print(
            f"[run_joint_dream] dropped {n_dropped}/{n_dropped + len(dfX)} "
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

    fitness = JointAlignGaussFitness(
        channels=channels, n_peaks=n_peaks, centers0=centers0, sigmas0=sigmas0,
        scale0=scale0, shift0=shift0,
        scale_bounds=scale_bounds, shift_bounds=shift_bounds,
        log_sigma_bounds=log_sigma_bounds, fit_baseline=fit_baseline,
    )
    problem = FitProblem(fitness)

    # avoid oversubscribing the 126 cores: each worker's numpy lstsq must
    # not also spawn its own BLAS thread pool.
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")

    result = fitters.fit(
        problem, method="dream", parallel=n_workers,
        samples=samples, burn=burn, thin=thin, pop=pop,
        outliers=outliers, trim=(outliers != "none"), verbose=True,
    )

    draw = result.state.draw()
    labels = draw.labels
    points = draw.points  # (n_draws, n_free_params)
    logp = draw.logp

    idx = {name: labels.index(name) for name in labels}

    def _summary(colname):
        x = points[:, idx[colname]]
        return {"median": float(np.median(x)), "lo16": float(np.percentile(x, 16)),
                "hi84": float(np.percentile(x, 84))}

    scale_summary = _summary("scale")
    shift_summary = _summary("shift")
    center_summaries = [_summary(f"center{k}") for k in range(n_peaks)]
    sigma_summaries = []
    for k in range(n_peaks):
        x = np.exp(points[:, idx[f"log_sigma{k}"]])
        sigma_summaries.append({
            "median": float(np.median(x)), "lo16": float(np.percentile(x, 16)),
            "hi84": float(np.percentile(x, 84)),
        })

    # set the Fitness parameters to the best point and reconstruct
    best_x, best_logp = result.state.best()
    for name, val in zip(labels, best_x):
        if name == "scale":
            fitness.scale.value = val
        elif name == "shift":
            fitness.shift.value = val
        elif name.startswith("center"):
            fitness.centers[int(name[len("center"):])].value = val
        elif name.startswith("log_sigma"):
            fitness.log_sigmas[int(name[len("log_sigma"):])].value = val

    out = {
        "problem": problem,
        "fitness": fitness,
        "fit_result": result,
        "state": result.state,
        "draw": draw,
        "best_fit": fitness.reconstruct(),
        "scale": scale_summary,
        "shift": shift_summary,
        "centers": center_summaries,
        "sigmas": sigma_summaries,
        "n_dropped_saxs_points": n_dropped,
    }

    cm._alignment["scale"] = scale_summary["median"]
    cm._alignment["shift"] = shift_summary["median"]
    cm._alignment["refinement"] = {
        "method": "joint_dream",
        "scale": scale_summary,
        "shift": shift_summary,
        "n_dropped_saxs_points": n_dropped,
    }
    return out


# ---------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------

def plot_scale_shift_posterior(result: Dict[str, Any]):
    """Cheap scale-vs-shift diagnostic. For the full posterior use
    `bumps.dream.views.plot_all(result['state'])` (needs matplotlib) or
    hand result['draw'].points / .labels to `corner.corner`."""
    import matplotlib.pyplot as plt

    labels = result["draw"].labels
    points = result["draw"].points
    scale = points[:, labels.index("scale")]
    shift = points[:, labels.index("shift")]
    logp = result["draw"].logp

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    axes[0].scatter(shift, scale, s=3, alpha=0.3, c=logp, cmap="viridis")
    axes[0].set_xlabel("shift")
    axes[0].set_ylabel("scale")
    axes[0].set_title("joint posterior (color = -nllf)")
    axes[1].hist(shift, bins=60)
    axes[1].set_xlabel("shift")
    axes[1].set_title("shift marginal")
    fig.tight_layout()
    return fig, axes


# ---------------------------------------------------------------------
# Optional: attach as a method on CoupledMeasurement
# ---------------------------------------------------------------------
#
#   from joint_align_gauss_dream import run_joint_dream
#   CoupledMeasurement.fit_gauss_peaks_dream = run_joint_dream
#
# then call it as   result = cm.fit_gauss_peaks_dream(n_peaks=3, n_workers=126)
