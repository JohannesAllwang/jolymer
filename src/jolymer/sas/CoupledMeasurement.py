import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.signal import correlate
from pathlib import Path
from dataclasses import dataclass
import re

from .SAXS_Measurement import SAXS_Measurement
from .ms import Ms
from jolymer.samples.bioMOLECULE import *
from jolymer.uv.onlineUV import onlineUV

def ms_from_folder(path='workdirs/ACGT6c80_waxs', file_prefix='waxs',
                   min_seqi=0, max_seqi=100000, q_beamstop=0.006,
                   exclude_seqi=[], frame_time=2.1,
                   angular_unit='A'):
    import os
    import re
    SRC_DIR = Path(path)
    pattern = re.compile(rf"^{file_prefix}_(\d{{1,3}})\.dat$")
    files = []
    for f in os.listdir(SRC_DIR):
        match = pattern.match(f)
        if match:
            num = int(match.group(1))
            files.append((num, f))
    files = sorted(files, key=lambda x: x[0])
    sorted_filenames = [f for num, f in files]
    I_list = []
    sigma_list = []
    q_values = None
    mlist = []
    for i, filename in enumerate(sorted_filenames):
        seqi = filename.split(f'{file_prefix}_')[1].split('.dat')[0]
        seqi = int(seqi)
        if seqi>max_seqi:
            break
        if seqi<min_seqi:
            continue
        if seqi in exclude_seqi:
            continue
        m = SAXS_Measurement(path=SRC_DIR, filename=filename,
                             qmin=q_beamstop, angular_unit=angular_unit)
        m.time = frame_time * float(seqi)
        m.seqi = seqi
        mlist.append(m)
    print("got number of ms:", len(mlist))
    out =  Ms(mlist)
    out.name = file_prefix
    out.min_seqi = min_seqi
    out.max_seqi = max_seqi
    return out

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import differential_evolution

@dataclass
class AlignSAXS_WAXS:

    ms_saxs: Ms
    ms_waxs: Ms

    def interpolate_matrix(
            self,
            qmin=None,
            qmax=None,
    ):
        q_s, I_s, sigma_s, _ = self.ms_saxs.build_saxs_matrix()
        q_w, I_w, sigma_w, _ = self.ms_waxs.build_saxs_matrix()
        outdict = {'q_s': q_s,
                   'q_w': q_w,
                   'I_s': I_s,
                   'I_w': I_w,
                   "sigma_s": sigma_s,
                   "sigma_w": sigma_w}
        qmin_auto = max(q_s.min(), q_w.min())
        qmax_auto = min(q_s.max(), q_w.max())
        print('auto q:', qmin_auto, qmax_auto)
        qmin = qmin if qmin is not None else qmin_auto
        qmax = qmax if qmax is not None else qmax_auto
        mask_w = (q_w > qmin) & (q_w < qmax)
        if mask_w.sum() < 3:
            raise ValueError("Overlap region too small")
        outdict['num_frames'] = I_w.shape[1]
        outdict['qw_ov'] = q_w[mask_w]
        outdict['Iw_ov'] = I_w[mask_w, :]
        interp = interp1d(q_s, I_s, axis=0, bounds_error=False,
                          fill_value=np.nan, kind='linear')
        Is_ov = interp(outdict['qw_ov'])
        sigma_w_ov = sigma_w[mask_w, :]
        interpw = interp1d(q_s, sigma_s, axis=0, bounds_error=False,
                           fill_value=np.nan, kind='linear')
        sigma_s_ov = interpw(outdict['qw_ov'])
        outdict['Is_ov'] = Is_ov
        outdict['sigma_s_ov'] = sigma_s_ov
        outdict['sigma_w_ov'] = sigma_w_ov
        return outdict

    def align_saxs_waxs_global(
        self,
        qmin=None,
        qmax=None,
        bg_bounds=0.5,
        use_errors=True,
        global_offset=False
    ):
        """
        Global alignment of WAXS to SAXS using all datasets simultaneously.
        Fits:
            I_waxs(q) = a * I_saxs(q) + b
        SAXS is unchanged; WAXS is scaled.
        Returns:
            a, b, diagnostics
        """
        matrix_dict = self.interpolate_matrix(qmin=qmin,
                                              qmax=qmax)
        def loss(params):
            a = params[0]
            offsets = params[1:]  # shape: (n_frames,)
            model = a * matrix_dict['Is_ov'] + offsets
            sigma2 = (
                matrix_dict['sigma_w_ov']**2 +
                (a * matrix_dict['sigma_s_ov'])**2
            )
            resid = matrix_dict['Iw_ov'] - model
            chi2 = np.sum((resid**2) / sigma2)
            return chi2
        mean_int = np.mean(matrix_dict['Iw_ov'])
        bounds = [(0.0, 10.0)] + [(-bg_bounds*np.abs(mean_int), bg_bounds*np.abs(mean_int))] *\
                matrix_dict['num_frames']
        if global_offset:
            bounds = [(0.0, 10.0)] + [(-bg_bounds*np.abs(mean_int), bg_bounds*np.abs(mean_int))]
        result = differential_evolution(loss, bounds,
                                        strategy='best1bin',
                                        recombination=0.7)
        return result.x[0], result.x[1:], result

    def rebin_overlap(
        self,
        a,
        b,
        qmin=None,
        qmax=None,
        minus_df=None,
        n_bins=200,
        prefix='merge'
    ):
        """
        Merge SAXS and aligned WAXS in the overlap region and
        store merged result in self.ms_saxs.ms.
        Parameters
        ----------
        a : float
            Global SAXS->WAXS scale factor
        b : ndarray
            Frame-wise offsets (n_frames,)
        """
        # --- unpack ---
        q_s, I_s, s_s, _ = self.ms_saxs.build_saxs_matrix()
        q_w, I_w, s_w, _ = self.ms_waxs.build_saxs_matrix()
        n_qs, n_frames = I_s.shape
        # --- apply alignment ---
        I_w_aligned = (I_w - b) / a
        s_w_aligned = s_w / abs(a)
        # --- overlap definition ---
        qmin_auto = max(q_s.min(), q_w.min())
        qmax_auto = min(q_s.max(), q_w.max())
        qmin = qmin if qmin is not None else qmin_auto
        qmax = qmax if qmax is not None else qmax_auto
        mask_s_ov = (q_s >= qmin) & (q_s <= qmax)
        mask_w_ov = (q_w >= qmin) & (q_w <= qmax)
        # --- rebin ---
        # --- define common q grid in overlap ---
        q_bins = np.linspace(qmin, qmax, n_bins + 1)
        q_centers = 0.5 * (q_bins[:-1] + q_bins[1:])
        # --- allocate merged arrays ---
        I_merge = np.full((n_bins, n_frames), np.nan)
        s_merge = np.full((n_bins, n_frames), np.nan)
        # --- frame-wise merge ---
        for t in range(n_frames):
            q_ov = np.concatenate([q_s[mask_s_ov], q_w[mask_w_ov]])
            I_ov = np.concatenate([
                I_s[mask_s_ov, t],
                I_w_aligned[mask_w_ov, t]
            ])
            s_ov = np.concatenate([
                s_s[mask_s_ov, t],
                s_w_aligned[mask_w_ov, t]
            ])
            for i in range(n_bins):
                m = (q_ov >= q_bins[i]) & (q_ov < q_bins[i + 1])
                if not np.any(m):
                    continue
                w = 1.0 / (s_ov[m] ** 2)
                I_merge[i, t] = np.sum(I_ov[m] * w) / np.sum(w)
                s_merge[i, t] = np.sqrt(1.0 / np.sum(w))
        # --- prepend / append pure SAXS (frame-wise) ---
        mask_lo = q_s < qmin
        mask_hi = q_w > qmax
        valid_bins = np.any(~np.isnan(I_merge), axis=1)
        q_centers = q_centers[valid_bins]
        I_merge = I_merge[valid_bins, :]
        s_merge = s_merge[valid_bins, :]
        q_final = np.concatenate([
            q_s[mask_lo],
            q_centers,
            q_w[mask_hi],
        ])
        I_final = np.vstack([
            I_s[mask_lo, :],
            I_merge,
            I_w_aligned[mask_hi, :],
        ])
        s_final = np.vstack([
            s_s[mask_lo, :],
            s_merge,
            s_w_aligned[mask_hi, :],
        ])
        # --- save per frame ---
        minusI = 0
        if not minus_df is None:
            minusI = minus_df.I
        for t, m in enumerate(self.ms_saxs.ms):
            path = Path(m.path)
            filename = f"{prefix}.{m.filename}"
            outdf = pd.DataFrame({
                "q": q_final,
                "I": I_final[:, t] - minusI,
                "err_I": s_final[:, t],
            })
            m.save_data(path / filename, df=outdf)
        return q_final, I_final, s_final

    def plot_alignment_diagnostics(
        self,
        a,
        b,
        qmin=None,
        qmax=None,
        chi2_ylim=(0, 10),
        use_mean=False
    ):
        """
        Diagnostic plot for SAXS/WAXS alignment.
        Shows per-frame reduced chi^2 and fitted offsets b_t.
        Parameters
        ----------
        a : float
            Global SAXS->WAXS scale factor
        b : ndarray
            Frame-wise offsets (n_frames,)
        """
        import matplotlib.pyplot as plt
        md = self.interpolate_matrix(qmin=qmin, qmax=qmax)
        Iw = md['Iw_ov']
        Is = md['Is_ov']
        sw = md['sigma_w_ov']
        ss = md['sigma_s_ov']
        n_q, n_frames = Iw.shape
        chi2_per_frame = np.zeros(n_frames)
        for i in range(n_frames):
            # model = a * Is[:, i] + b.mean()
            if use_mean:
                model = a * Is[:, i] + b.mean()
            else:
                model = a * Is[:, i] + b[i]
            resid = Iw[:, i] - model
            sigma2 = sw[:, i]**2 + (a * ss[:, i])**2
            chi2_per_frame[i] = np.mean(resid**2 / sigma2)
        # --- plotting ---
        fig, ax1 = plt.subplots(figsize=(6, 2.5))
        ax1.plot(
            chi2_per_frame,
            color="tab:red",
            lw=1.5,
            label=r"$\chi^2_\mathrm{red}$",
        )
        ax1.set_xlabel("Frame number")
        ax1.set_ylabel(r"Reduced $\chi^2$", color="tab:red")
        ax1.tick_params(axis="y", labelcolor="tab:red")
        ax1.set_ylim(*chi2_ylim)
        ax2 = ax1.twinx()
        ax2.plot(
            b,
            color="tab:blue",
            alpha=0.5,
            lw=1.2,
            label=r"Offset $b_t$",
        )
        ax2.set_ylabel("Background offset $b$", color="tab:blue")
        ax2.tick_params(axis="y", labelcolor="tab:blue")
        plt.title("SAXS–WAXS Global Alignment Diagnostics")
        fig.tight_layout()
        return chi2_per_frame


_AUTORG_LINE_RE = re.compile(
    r"^\s*(?P<Rg>[\d.]+)\s+(?P<Rg_stdev_pct>\d+)%\s+"
    r"(?P<I0>[\d.eE+-]+)\s+"
    r"(?P<g_first>\d+)\s*-\s*(?P<g_last>\d+)\s*\(\s*(?P<g_npoints>\d+)\s*\)\s+"
    r"(?P<quality_pct>\d+)%\s*(?P<flag>[A-Za-z]?)\s+"
    r"(?P<file>\S+)\s*$"
)


def parse_autorg_table(path):
    """
    Parse an ATSAS `autorg` batch-summary file into a DataFrame.
    Source columns: Rg, stdev(%), I(0), Guinier points, Quality(%), File.
    There is no direct I(0) error in this table - only a *relative* error
    on Rg (stdev, %) and an overall fit Quality (%). `Rg_err` is derived as
    Rg * Rg_stdev_pct / 100 for convenience; `I0_err` is left as NaN.
    Important: the `File` column is frequently NOT just the bare filename -
    ATSAS includes the parent folder (e.g. "B3/B3_0.dat") whenever it needs
    to disambiguate or fit the column width. `filename` here is the
    basename only, and is what should be used to match against your own
    measurement objects' filenames (see `load_autorg` below) - matching on
    the full `file` string will silently fail to find any rows whenever
    that prefix is present.
    """
    rows = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            m = _AUTORG_LINE_RE.match(line)
            if m is None:
                continue  # header line or anything unparseable - skip it
            d = m.groupdict()
            rows.append({
                "Rg": float(d["Rg"]),
                "Rg_stdev_pct": float(d["Rg_stdev_pct"]),
                "I0": float(d["I0"]),
                "guinier_first": int(d["g_first"]),
                "guinier_last": int(d["g_last"]),
                "guinier_npoints": int(d["g_npoints"]),
                "quality_pct": float(d["quality_pct"]),
                "aggregation_flag": d["flag"].upper() == "A",
                "file": d["file"],
                "filename": Path(d["file"]).name,
            })
    if not rows:
        raise ValueError(f"No parseable autorg rows found in {path}")
    df = pd.DataFrame(rows)
    df["Rg_err"] = df["Rg"] * df["Rg_stdev_pct"] / 100.0
    df["I0_err"] = np.nan
    return df

@dataclass
class CoupledMeasurement:
    """
    Couples SAXS and online UV measurements recorded at the same spot.
    Handles alignment, interpolation, and joint diagnostics.
    """
    saxs_list: Ms
    uv: onlineUV
    sample: bioMOLECULE
    qmin: float = 0.001
    qmax: float = 0.5
    qstar: float = 0.1
    saxs_kind: str = 'autorg'
    autorg_filename: str = 'autorg.dat'

    def __post_init__(self):
        self._I = None
        self._sigma = None
        self._q = None
        self._x = None
        self._alignment = {}
        self._uv_on_saxs = None
        self._alignment["saxs_kind"] = self.saxs_kind
        self._alignment["qstar"] = self.qstar
        self._alignment["qmin"] = self.qmin
        self._alignment["qmax"] = self.qmax
        self._alignment["shift"] = 0
        self._alignment["scale"] = 1

    def build_saxs_matrix(self):
        I_list = []
        sigma_list = []
        q_values = None
        for m in self.saxs_list:
            df = m.get_data()
            # df = df.iloc[0:]  # your skip logic
            if q_values is None:
                q_values = df.q.to_numpy()
            I_list.append(df.I.to_numpy())
            sigma_list.append(df.err_I.to_numpy())
        self._I = np.column_stack(I_list)
        self._sigma = np.column_stack(sigma_list)
        self._q = q_values
        self._x = np.arange(self._I.shape[1]) * 2.1
        return self._q, self._I, self._sigma, self._x

    def get_saxs_scalar(self, kind="autorg", qstar=0.1,
                        autorg_filename=None,
                        qmin=0.04, qmax=0.4, load=False):
        """
        Returns frame, time, I0 (+ optionally errI0) for the saxs series.
        New: kind="autorg" sources I0/Rg from ATSAS `autorg` output
        (`self._autorg_df`, built by `load_autorg`) instead of the internal
        Guinier-region fit in `self.saxs_list.get_rgs(...)`. This is much
        faster (autorg already ran once from the GUI) and matches whatever
        autorg itself decided the Guinier region was, rather than the fixed
        (qmin, qmax) bounds / bounded fit used by kind="I0".
        """
        if autorg_filename is None:
            autorg_filename = self.autorg_filename
        q = self._q if self._q is not None else self.build_saxs_matrix()[0]
        if kind == "autorg":
            autorg_path = Path(self.get_saxs_path()) / autorg_filename
            if not autorg_path.exists():
                self.run_autorg()
            autorg_df = self.load_autorg(autorg_filename=autorg_filename)
            frames = [getattr(m, "seqi", i) for i, m in enumerate(self.saxs_list)]
            autorg_df.insert(0, "frame", frames)
            return autorg_df[["frame", "time", "I0", "errI0"]]
        elif kind == "I0":
            frames = []
            times = []
            rgdf = self.saxs_list.get_rgs(qmin=qmin, qmax=qmax,
                                          load=load, I00=0.01, Rg0=12,
                                          bounds=((10, 0), (20, 1.0)))
            for frame, m in enumerate(self.saxs_list):
                frames.append(m.seqi)
                times.append(m.time)
            return pd.DataFrame({
                'frame': frames,
                'time': times,
                'I0': rgdf.I0,
                'errI0': rgdf.err_I0})
        elif kind == "Iq":
            idx = np.argmin(np.abs(q - qstar))
            frames = []
            times = []
            I0s = []
            errI0s = []
            for frame, m in enumerate(self.saxs_list):
                frames.append(m.seqi)
                times.append(m.time)
                I0s.append(float(m.get_data()['I'][idx]))
                errI0s.append(float(m.get_data()['err_I'][idx]))
            return pd.DataFrame({
                'frame': frames,
                'time': times,
                'I0': I0s,
                'errI0': errI0s})
        else:
            raise ValueError(f"Unknown SAXS scalar: {kind}")

    def reload_saxs(self, autorg_filename=None):
        """
        Permanently move the UV<->SAXS time alignment onto the SAXS
        measurements themselves: each m.time is remapped from raw SAXS
        time into UV time units, and self._alignment is reset to the
        identity (scale=1, shift=0).
        Only the time axis is touched -- I, sigma, Rg, I0 etc. on each
        measurement are left exactly as they were. Call this once you're
        happy with an alignment (e.g. after align_uv_to_saxs /
        refine_uv_alignment*), so that afterwards SAXS and UV share a
        single, UV-derived time axis and no further scale/shift is needed
        (fit_gauss_peaks, get_saxs_scalar, etc. all just work in that
        common frame).
        Returns the (now identity) alignment dict.
        """
        if autorg_filename is None:
            autorg_filename = self.autorg_filename
        scale = self._alignment["scale"]
        shift = self._alignment["shift"]
        if scale == 1.0 and shift == 0.0:
            # nothing to do -- avoid silently "reloading" for no reason
            return self._alignment
        for m in self.saxs_list:
            m.time = (m.time - shift) / scale
        # anything cached against the old time axis is now stale
        self._uv_on_saxs = None
        self._autorg_df = None
        self._alignment["scale"] = 1.0
        self._alignment["shift"] = 0.0
        self._alignment["reloaded_from"] = {"scale": scale, "shift": shift}
        # refresh the cached autorg table (frame/time columns) against
        # the new m.time values, if that's the scalar source in use
        if self._alignment.get("saxs_kind", self.saxs_kind) == "autorg":
            self.load_autorg(autorg_filename=autorg_filename)
        return self._alignment


    def preprocess_uv(dfUV, mode="identity", **kwargs):
        """
        Preprocess UV dataframe before alignment.
        Extend this later with deconvolution, R^6 scaling, etc.
        """
        df = dfUV.copy()
        if mode == "identity":
            return df
        elif mode == "custom":
            # placeholder for future logic
            # e.g. df['Abs'] = deconvolve(df['Abs'], ...)
            # e.g. df['Abs'] *= df['R']**6
            raise NotImplementedError("Custom UV preprocessing not implemented")
        else:
            raise ValueError(f"Unknown UV preprocessing mode: {mode}")

    def get_dfUV(self):
        dfUV = self.uv.get_scaled_Abs(refwl=self.uv.refwl,
                                      outwl=self.uv.outwl,
                                      alignment_time=self.uv.alignment_time)
        return dfUV

    def reset_times(self):
        if not hasattr(self, "_t_raw"):
            self._t_raw = [m.time for m in self.saxs_list]
        for m, t in zip(self.saxs_list, self._t_raw):
            m.time = t
        self._uv_on_saxs = None
        self._autorg_df = None
        self._alignment.update(scale=1.0, shift=0.0)
        return self._alignment

    @staticmethod
    def align_uv_df_to_saxs_df(dfUV, dfSAXS, uv_col="Abs", saxs_col="I0", min_t_saxs=0, max_t_saxs=np.inf,
                               I0_min=0, scale0shift0=(2.2, -1400.0), scale_rel=0.2, shift_window=600.0,
                               n_scale=81, n_shift=241, **_):
        """Model: t_saxs = scale * t_uv + shift. Deterministic grid search around scale0shift0, then Nelder-Mead polish."""
        from scipy.optimize import minimize
        ts = dfSAXS["time"].to_numpy(float)
        ys = dfSAXS[saxs_col].to_numpy(float)
        m = np.isfinite(ys) & (ts >= min_t_saxs) & (ts <= max_t_saxs) & (ys >= I0_min)
        ts, ys = ts[m], ys[m]
        ys = (ys - ys.mean()) / ys.std()
        tu = dfUV["time"].to_numpy(float)
        yu = dfUV[uv_col].to_numpy(float)
        o = np.argsort(tu)
        tu, yu = tu[o], yu[o]
        def score(p):
            if p[0] <= 0:
                return -2.0
            yi = np.interp(ts, p[0] * tu + p[1], yu, left=np.nan, right=np.nan)
            ok = np.isfinite(yi)
            if ok.sum() < 10 or np.std(yi[ok]) == 0:
                return -2.0
            return np.corrcoef(ys[ok], yi[ok])[0, 1]
        a0, b0 = scale0shift0
        best = (-np.inf, a0, b0)
        for a in np.linspace(a0 * (1 - scale_rel), a0 * (1 + scale_rel), n_scale):
            for b in np.linspace(b0 - shift_window, b0 + shift_window, n_shift):
                s = score((a, b))
                if s > best[0]:
                    best = (s, a, b)
        _, a, b = best
        res = minimize(lambda p: -score(p), [a, b], method="Nelder-Mead",
                       options=dict(initial_simplex=[[a, b], [a * 1.002, b], [a, b + 2.0]], xatol=1e-6, fatol=1e-9))
        a, b = res.x
        return {"method": "affine", "scale": float(a), "shift": float(b), "score": float(-res.fun), "success": bool(res.success)}

    def align_uv_to_saxs(self, saxs_kind="autorg", qstar=0.1, qmin=0.04, qmax=0.4, load=True,
                         min_t_saxs=0, max_t_saxs=np.inf, I0_min=0, scale0shift0=(2.2, -1400.0), **kwargs):
        self.reset_times()
        dfUV = self.get_dfUV()
        self.uv.df = dfUV
        dfSAXS = self.get_saxs_scalar(kind=saxs_kind, qstar=qstar, qmin=qmin, qmax=qmax, load=load)
        alignment = self.align_uv_df_to_saxs_df(dfUV, dfSAXS, min_t_saxs=min_t_saxs, max_t_saxs=max_t_saxs,
                                                I0_min=I0_min, scale0shift0=scale0shift0, **kwargs)
        alignment.update(saxs_kind=saxs_kind, qstar=qstar, qmin=qmin, qmax=qmax)
        self._alignment = alignment
        return self._alignment

    def get_saxs_path(self):
        saxs_path = self.saxs_list.ms[0].path
        return Path(saxs_path)

    def run_autorg(self, autorg_name="autorg.dat", pattern="GR_*.dat",
                   directory=None):
        """
        Run ATSAS autorg on all matching files and write the combined
        output to a CSV file.
        Parameters
        ----------
        directory : str or pathlib.Path
            Directory containing the GR_*.dat files.
        output : str or pathlib.Path
            Output CSV filename.
        pattern : str
            Input file pattern.
        """
        import subprocess
        if directory is None:
            saxs_path = self.get_saxs_path()
        else:
            saxs_path = directory
        output = directory / autorg_name
        files = sorted(directory.glob(pattern))
        if not files:
            raise FileNotFoundError(
                f"No files matching {pattern!r} found in {directory}"
            )
        subprocess.run(
            ["autorg", *map(str, files), "-o", str(output)],
            check=True,
        )
        return output

    def load_autorg(self, autorg_filename=None, directory=None,
                    wavelength=[]):
        """
        Read an autorg.dat produced by `run_autorg` and attach Rg/I0 to each
        measurement in `self.saxs_list`.
        Parses the whole table once with `parse_autorg_table` (see above) and
        matches rows to measurements by **basename**, since the table's `File`
        column often includes a parent-folder prefix (e.g. "B3/B3_0.dat") that
        won't equal a bare `m.filename` - matching on the full string was
        silently dropping every row (Rg/I0 all NaN) whenever that prefix was
        present.
        """
        if autorg_filename is None:
            autorg_filename = self.autorg_filename
        saxs_path = Path(directory) if directory is not None else self.get_saxs_path()
        autorg_path = Path(saxs_path, autorg_filename)
        table = parse_autorg_table(autorg_path).set_index("filename")
        autorg_dict = {"time": [], "Rg": [], "I0": [], "errRg": [], "errI0": []}
        for m in self.saxs_list:
            key = Path(m.filename).name  # normalize in case m.filename ever carries a folder too
            if key in table.index:
                row = table.loc[key]
                if isinstance(row, pd.DataFrame):
                    raise ValueError(
                        f"Ambiguous autorg match for {key!r}: multiple rows "
                        f"share this basename ({list(row['file'])}). Use full "
                        f"relative paths for m.filename, or dedupe the input files."
                    )
                m.autorg = [row.Rg, row.I0]
                err_rg = row.Rg_err
                err_i0 = row.I0_err
            else:
                m.autorg = [np.nan, np.nan]
                err_rg = np.nan
                err_i0 = np.nan
            autorg_dict["time"].append(m.time)
            autorg_dict["Rg"].append(m.autorg[0])
            autorg_dict["I0"].append(m.autorg[1])
            autorg_dict["errRg"].append(err_rg)
            autorg_dict["errI0"].append(err_i0)
        self._autorg_df = pd.DataFrame(autorg_dict)
        return self._autorg_df

    def interpolate_uv(self):
        if self._alignment is None:
            raise RuntimeError("UV–SAXS alignment not computed")
        shift = self._alignment["shift"]
        scale = self._alignment["scale"]
        uv_time = scale*(self.uv.df.time) + shift
        uv_abs = self.uv.df['Abs']
        interp = interp1d(
            uv_time,
            uv_abs,
            bounds_error=False,
            fill_value=0.0,
            kind='cubic',
        )
        saxs_time = np.array([m.time for m in self.saxs_list])
        uv_on_saxs = interp(saxs_time)
        self._uv_on_saxs = pd.DataFrame({
            "frame": self._x,
            "time": saxs_time,
            "Abs": uv_on_saxs,
        })
        return self._uv_on_saxs

    def refine_uv_alignment(
        self,
        saxs_kind="autorg",
        qstar=0.1,
        qmin=0.04,
        qmax=0.4,
        load=True,
        min_t_saxs=0,
        max_t_saxs=np.inf,
        I0_min=0,
        shift_bounds=(-100, 100),
    ):
        dfUV = self.uv.get_scaled_Abs(
            refwl=self.uv.refwl,
            outwl=self.uv.outwl,
            alignment_time=self.uv.alignment_time,
        )
        self.uv.df = dfUV
        dfSAXS = self.get_saxs_scalar(
            kind=saxs_kind,
            qstar=qstar,
            qmin=qmin,
            qmax=qmax,
            load=load,
        )
        refinement = self.refine_uv_shift(
            dfUV,
            dfSAXS,
            min_t_saxs=min_t_saxs,
            max_t_saxs=max_t_saxs,
            I0_min=I0_min,
            shift_bounds=shift_bounds,
        )
        # Keep the existing alignment information
        self._alignment["shift"] = refinement["shift"]
        self._alignment["scale"] = 1.0
        self._alignment["refinement"] = refinement
        return self._alignment

    def interpolate_uv_matrix(self, uv_wide, dt=0.946, mintime=None, maxtime=None):
        if self._alignment is None:
            raise RuntimeError("UV–SAXS alignment not computed")
        shift = self._alignment["shift"]
        scale = self._alignment["scale"]
        wl = uv_wide.iloc[:,0].values
        abs_matrix = uv_wide.iloc[:,1:].values
        ntime = abs_matrix.shape[1]
        raw_time = np.linspace(0, ntime, ntime) * dt
        mask = np.ones(ntime, dtype=bool)
        if mintime is not None:
            mask &= raw_time > mintime
        if maxtime is not None:
            mask &= raw_time < maxtime
        raw_time = raw_time[mask]
        abs_matrix = abs_matrix[:,mask]
        uv_time = scale*raw_time + shift
        interp = interp1d(uv_time, abs_matrix, axis=1, bounds_error=False, fill_value=0.0,
                          kind='cubic')
        saxs_time = np.array([m.time for m in self.saxs_list])
        uv_on_saxs = interp(saxs_time)
        cols = ["wl"] + [f"Abs{i}" for i in range(len(saxs_time))]
        return pd.DataFrame(np.column_stack([wl, uv_on_saxs]), columns=cols)

    def _scan_shift(self, theta, n, blocks, bounds, baseline, num):
        t = theta.copy()
        best = (np.inf, theta[2 * n])
        for sh in np.linspace(*bounds, num):
            t[2 * n] = sh
            cost = np.sum(self._eval_blocks(t, n, blocks[:1], True, baseline)[0][2] ** 2)
            if cost < best[0]:
                best = (cost, sh)
        return best[1]

    @staticmethod
    def _block_report(b, G, A, n, baseline):
        cols = [f"peak_{k + 1}" for k in range(n)]
        base = A[n] - A[n + 1] if baseline else np.zeros(A.shape[1])
        return dict(
            time=b["t"],
            data=pd.DataFrame(b["Y"], index=b["t"], columns=b["names"]),
            fit=pd.DataFrame(G @ A, index=b["t"], columns=b["names"]),
            amplitudes=pd.DataFrame(A[:n].T, index=b["names"], columns=cols),
            baselines=pd.Series(base, index=b["names"]),
            peaks={nm: {cols[k]: G[:, k] * A[k, j] for k in range(n)} for j, nm in enumerate(b["names"])})

    def fit_gauss_peaks(self, n_peaks, wl=280, p0=None, baseline=False, fit_alignment=True,
                        fit_scale=False, shift_window=100.0, scale_tol=0.05, uv_weight=1.0,
                        scan_points=401, apply=False, sigma_bounds=None, **ls_kwargs):
        """
        Joint fit of SAXS I0 and UV absorbance (one or many wavelengths) as a sum of
        Gaussians on their native time grids, no interpolation. UV centers/sigmas are
        shared; the SAXS block sees them as c' = scale*c + shift, s' = scale*s.
        shift/scale bounds are relative to the current self._alignment.
        apply=True stores the fitted (scale, shift) in self._alignment; call
        reload_saxs() afterwards to move SAXS onto the UV time axis.
        """
        from scipy.optimize import least_squares
        n = n_peaks
        a = self._alignment
        scale0, shift0 = a.get("scale", 1.0), a.get("shift", 0.0)
        blocks = self._build_blocks(wl, uv_weight)
        theta0, (lo, hi) = self._theta_start_bounds(
            blocks, n, True, (scale0 * (1 - scale_tol), scale0 * (1 + scale_tol)),
            (shift0 - shift_window, shift0 + shift_window), sigma_bounds)
        if p0 is not None:
            theta0[:n] = p0["centers"]
            theta0[n:2 * n] = np.log(p0["sigmas"])
        theta0 = np.clip(theta0, lo, hi)
        free = np.ones(len(theta0), bool)
        if not fit_alignment:
            free[2 * n:] = False
        elif not fit_scale:
            free[2 * n + 1] = False
        if fit_alignment:
            theta0[2 * n] = self._scan_shift(theta0, n, blocks, (lo[2 * n], hi[2 * n]), baseline, scan_points)
        def fun(x):
            th = theta0.copy()
            th[free] = x
            return self._joint_residual(th, n, blocks, True, baseline)
        res = least_squares(fun, theta0[free], bounds=(lo[free], hi[free]), x_scale="jac", **ls_kwargs)
        theta = theta0.copy()
        theta[free] = res.x
        c, s, scale, shift = self._unpack_theta(theta, n, True)
        ev = self._eval_blocks(theta, n, blocks, True, baseline)
        saxs, uv = (self._block_report(b, G, A, n, baseline) for b, (G, A, _) in zip(blocks, ev))
        wl_arr = np.atleast_1d(np.asarray(list(wl) if isinstance(wl, range) else wl, float))
        uv["wl"] = wl_arr if len(wl_arr) == len(uv["data"].columns) else np.arange(len(uv["data"].columns))
        if apply and fit_alignment:
            self._store_alignment(scale, shift, 2 * res.cost, res.success)
        return dict(result=res, cost=float(2 * res.cost), centers=c, sigmas=s, scale=float(scale),
                    shift=float(shift), centers_saxs=scale * c + shift, sigmas_saxs=scale * s,
                    saxs=saxs, uv=uv)

    def plot_gauss_peaks(self, fit, max_uv=3, channels=None):
        import matplotlib.pyplot as plt
        names = list(fit["uv"]["data"].columns)
        if channels is None:
            channels = [names[i] for i in np.unique(np.linspace(0, len(names) - 1, min(max_uv, len(names))).astype(int))]
        panels = [("I0", fit["saxs"], lambda t: (t - fit["shift"]) / fit["scale"])]
        panels += [(c, fit["uv"], lambda t: t) for c in channels]
        fig, axes = plt.subplots(len(panels), 1, figsize=(7, 3 * len(panels)), sharex=True, squeeze=False)
        for ax, (ch, rep, tf) in zip(axes[:, 0], panels):
            t = tf(rep["time"])
            ax.plot(t, rep["data"][ch], label=ch)
            for k, v in rep["peaks"][ch].items():
                ax.plot(t, v, ":", label=k)
            ax.plot(t, rep["fit"][ch], label="fit")
            ax.set_ylabel(ch)
            ax.legend()
        axes[-1, 0].set_xlabel("UV time (s)")
        fig.tight_layout()
        return fig, axes[:, 0]

    def plot_gauss_spectra(self, fit):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 3))
        A = fit["uv"]["amplitudes"]
        for k in A.columns:
            yval = A[k].to_numpy()
            yval = yval/yval.max()
            ax.plot(fit["uv"]["wl"], yval, label=k)
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Peak amplitude")
        ax.legend()
        return fig, ax


    def refine_uv_alignment_scale_shift(
        self,
        saxs_kind="autorg",
        qstar=0.1,
        qmin=0.04,
        qmax=0.4,
        load=True,
        min_t_saxs=0,
        max_t_saxs=np.inf,
        I0_min=0,
        scale_bounds=(0.98, 1.02),
        shift_bounds=(-10.0, 10.0),
    ):
        """Instance wrapper for refine_uv_scale_shift, mirroring
        refine_uv_alignment's wrapping of refine_uv_shift."""
        dfUV = self.get_dfUV()
        self.uv.df = dfUV
        dfSAXS = self.get_saxs_scalar(
            kind=saxs_kind, qstar=qstar, qmin=qmin, qmax=qmax, load=load,
        )
        refinement = self.refine_uv_scale_shift(
            dfUV,
            dfSAXS,
            min_t_saxs=min_t_saxs,
            max_t_saxs=max_t_saxs,
            I0_min=I0_min,
            scale0=self._alignment.get("scale", 1.0),
            shift0=self._alignment.get("shift", 0.0),
            scale_bounds=scale_bounds,
            shift_bounds=shift_bounds,
        )
        self._alignment["shift"] = refinement["shift"]
        self._alignment["scale"] = refinement["scale"]
        self._alignment["refinement"] = refinement
        return self._alignment

    def to_regals(self):
        dfX = self.get_saxs_scalar(self._alignment["saxs_kind"],
                                      qmin=self._alignment["qmin"],
                                      qmax=self._alignment["qmax"],
                                      qstar=self._alignment["qstar"],
                                   load=True
                                      )
        uv_time = self._uv_on_saxs["time"]
        uv_y = self._uv_on_saxs["Abs"]
        outdict = {
                'x': np.array(dfX.time),
                'idx': np.array(dfX.frame),
                'q': self._q,
                'I': self._I,
                'sigma': self._sigma,
                'uv_meas': uv_y}
        return outdict

    def plot_saxs_scalar(self):
        import matplotlib.pyplot as plt
        dfX = self.get_saxs_scalar(self._alignment["saxs_kind"],
                                      qmin=self._alignment["qmin"],
                                      qmax=self._alignment["qmax"],
                                      qstar=self._alignment["qstar"],
                                   load=True
                                      )
        fig, ax = plt.subplots()
        ax.plot(dfX.time, dfX.I0, label="SAXS")
        ax.legend()

    def plot_alignment(self, ax=None):
        import matplotlib.pyplot as plt
        dfX = self.get_saxs_scalar(self._alignment["saxs_kind"],
                                      qmin=self._alignment["qmin"],
                                      qmax=self._alignment["qmax"],
                                      qstar=self._alignment["qstar"],
                                   load=True
                                      )
        uv_time = self._uv_on_saxs["time"]
        uv_y = self._uv_on_saxs["Abs"]
        if ax is None:
            fig, ax = plt.subplots()
        ax.plot(dfX.time, dfX.I0 / dfX.I0.max(), label="SAXS")
        ax.plot(dfX.time, uv_y / uv_y.max(), label="UV")
        ax.legend()
        return ax

    def plot_uv_autorg(
        self,
        autorg_name=None,
        pattern="GR_*.dat",
        wl_min=260,
        wl_max=290,
        refwl=280,
        xlim=None,
        ylim_rg=None,
        ylim_i0=None,
        ylim_uv=None,
        name=None,
        minutes=False,
        ax=None,
    ):
        if autorg_name is None:
            autorg_name = self.autorg_filename
        if ax is None:
            fig, ax = plt.subplots(figsize=(7, 3))
        autorg_path = Path(self.get_saxs_path()) / autorg_name
        if not autorg_path.exists():
            self.run_autorg(autorg_name=autorg_name, pattern=pattern)
        autorg_df = self.load_autorg(autorg_filename=autorg_name)
        autorg_df.time = autorg_df.time / self._alignment['scale'] - self._alignment['shift'] / self._alignment['scale']
        axU = ax.twinx()
        self.uv.plot_full_wavelength_map(
            wl_min=wl_min,
            wl_max=wl_max,
            refwl=refwl,
            ax=axU,
            alignment_time=self.uv.alignment_time,
        )
        if ylim_uv is not None:
            axU.set_ylim(*ylim_uv)
        ax.plot(autorg_df.time, autorg_df.Rg, color="green")
        ax.set_ylabel(r"$R_g$ (nm)", color="green")
        ax.tick_params(axis="y", colors="green")
        if ylim_rg is not None:
            ax.set_ylim(*ylim_rg)
        axI = ax.twinx()
        axI.spines["right"].set_position(("outward", 50))
        axI.plot(autorg_df.time, autorg_df.I0, color="black")
        axI.set_ylabel(r"$I(0)$", color="black")
        axI.tick_params(axis="y", colors="black")
        if ylim_i0 is not None:
            axI.set_ylim(*ylim_i0)
        ax.set_xlabel("Time (s)")
        if xlim is not None:
            ax.set_xlim(*xlim)
        if name is not None:
            ax.annotate(name, xy=(0.5, 0.8), xycoords="axes fraction")
        if minutes:
            ax.xaxis.set_major_formatter(
                matplotlib.ticker.FuncFormatter(lambda x, pos: f"{x / 60:.0f}")
            )
            ax.set_xlabel("Time (min)")
        return ax, autorg_df


    @staticmethod
    def refine_uv_shift(
        dfUV,
        dfSAXS,
        uv_col="Abs",
        saxs_col="I0",
        min_t_saxs=0,
        max_t_saxs=np.inf,
        I0_min=0,
        shift0=0.0,
        shift_bounds=(-10.0, 10.0),
    ):
        """
        Refine an already approximately aligned UV/SAXS time axis
        by optimizing only a small time shift.
        The UV time axis is transformed as:
            t_uv_corrected = t_uv + shift
        No scale correction is performed.
        """
        from scipy.optimize import minimize_scalar
        t_saxs = dfSAXS["time"].to_numpy()
        y_saxs = dfSAXS[saxs_col].to_numpy()
        mask = (
            (t_saxs >= min_t_saxs) &
            (t_saxs <= max_t_saxs) &
            (y_saxs >= I0_min)
        )
        t_saxs = t_saxs[mask]
        y_saxs = y_saxs[mask]
        t_uv = dfUV["time"].to_numpy()
        y_uv = dfUV[uv_col].to_numpy()
        # Normalize SAXS once
        y_saxs_n = (
            (y_saxs - y_saxs.mean()) /
            y_saxs.std()
        )
        def loss(shift):
            f = interp1d(
                t_uv + shift,
                y_uv,
                kind="cubic",
                bounds_error=False,
                fill_value=0,
            )
            y_uv_i = f(t_saxs)
            if np.std(y_uv_i) == 0:
                return np.inf
            y_uv_n = (
                (y_uv_i - y_uv_i.mean()) /
                y_uv_i.std()
            )
            return -np.sum(y_saxs_n * y_uv_n)
        res = minimize_scalar(
            loss,
            bounds=shift_bounds,
            method="bounded",
            options={"xatol": 0.01},
        )
        shift = float(res.x)
        score = -res.fun
        return {
            "method": "shift_refinement",
            "scale": 1.0,
            "shift": shift,
            "score": score,
            "success": res.success,
        }

    @staticmethod
    def refine_uv_scale_shift(
        dfUV,
        dfSAXS,
        uv_col="Abs",
        saxs_col="I0",
        min_t_saxs=0,
        max_t_saxs=np.inf,
        I0_min=0,
        scale0=1.0,
        shift0=0.0,
        scale_bounds=(0.98, 1.02),
        shift_bounds=(-10.0, 10.0),
    ):
        """
        Local refinement of (scale, shift) together. Same correlation
        objective as align_uv_df_to_saxs_df, but started from
        (scale0, shift0) and optimized with L-BFGS-B instead of
        differential_evolution -- appropriate once the alignment is
        already close and you just want a small, fast, reproducible
        correction rather than a full global search.
        """
        from scipy.optimize import minimize

        t_saxs = dfSAXS["time"].to_numpy()
        y_saxs = dfSAXS[saxs_col].to_numpy()
        mask = (
            (t_saxs >= min_t_saxs) &
            (t_saxs <= max_t_saxs) &
            (y_saxs >= I0_min)
        )
        t_saxs_m = t_saxs[mask]
        y_saxs_m = y_saxs[mask]
        t_uv = dfUV["time"].to_numpy()
        y_uv = dfUV[uv_col].to_numpy()

        y_saxs_n = (y_saxs_m - y_saxs_m.mean()) / y_saxs_m.std()

        def loss(p):
            a, b = p
            f = interp1d(
                a * t_uv + b,
                y_uv,
                kind="cubic",
                bounds_error=False,
                fill_value=0,
            )
            y_uv_i = f(t_saxs_m)
            if np.std(y_uv_i) == 0:
                return np.inf
            y_uv_n = (y_uv_i - y_uv_i.mean()) / y_uv_i.std()
            return -np.sum(y_saxs_n * y_uv_n)

        res = minimize(
            loss,
            x0=[scale0, shift0],
            bounds=[scale_bounds, shift_bounds],
            method="L-BFGS-B",
        )
        a, b = res.x
        score = -res.fun
        return {
            "method": "scale_shift_refinement",
            "scale": float(a),
            "shift": float(b),
            "score": score,
            "success": res.success,
        }

    def _gauss_peaks_cost_at_alignment(self, scale, shift, n_peaks, wl, p0, fit_baseline):
        """Run fit_gauss_peaks at a trial (scale, shift), restoring the
        previous alignment afterwards. Returns (cost, fit_dict)."""
        saved_scale = self._alignment["scale"]
        saved_shift = self._alignment["shift"]
        self._alignment["scale"] = scale
        self._alignment["shift"] = shift
        try:
            fit = self.fit_gauss_peaks(n_peaks, wl=wl, p0=p0, fit_baseline=fit_baseline)
        finally:
            self._alignment["scale"] = saved_scale
            self._alignment["shift"] = saved_shift
        cost = float(np.sum(fit["result"].fun ** 2))
        return cost, fit

    def run_joint_dream(self, *args, **kwargs):
        from jolymer.sas.joint_align_gauss_dream import run_joint_dream
        return run_joint_dream(self, *args, **kwargs)

    @staticmethod
    def _gauss_design(t, centers, sigmas, baseline=False):
        G = np.exp(-(t[:, None] - centers[None, :])**2 / (2 * sigmas[None, :]**2))
        return np.column_stack([G, np.ones_like(t), -np.ones_like(t)]) if baseline else G

    @staticmethod
    def _solve_amplitudes(G, Y):
        from scipy.optimize import nnls
        if np.isfinite(Y).all():
            A = np.linalg.lstsq(G, Y, rcond=None)[0]
            bad = np.where((A < 0).any(axis=0))[0]
        else:
            A = np.zeros((G.shape[1], Y.shape[1]))
            bad = range(Y.shape[1])
        for j in bad:
            m = np.isfinite(Y[:, j])
            A[:, j] = nnls(G[m], Y[m, j])[0]
        return A

    @staticmethod
    def _unpack_theta(theta, n, fit_scale):
        scale = theta[2 * n + 1] if fit_scale else 1.0
        return theta[:n], np.exp(theta[n:2 * n]), scale, theta[2 * n]

    @staticmethod
    def _guess_peaks(t, y, n):
        from scipy.signal import find_peaks, peak_widths
        y = np.nan_to_num(np.asarray(y, float))
        dt = np.median(np.diff(t))
        idx, prop = find_peaks(y, prominence=0.02 * np.ptp(y))
        idx = idx[np.argsort(prop["prominences"])[::-1][:n]]
        sig = peak_widths(y, idx, rel_height=0.5)[0] * dt / 2.355 if len(idx) else []
        c = list(t[idx])
        s = list(np.maximum(sig, dt))
        k = n - len(c)
        if k > 0:
            c += list(np.linspace(t[0], t[-1], k + 2)[1:-1])
            s += [(t[-1] - t[0]) / (4 * n)] * k
        o = np.argsort(c)
        return np.array(c)[o], np.array(s)[o]

    def _build_blocks(self, wl=280, uv_weight=1.0):
        a = self._alignment
        dfX = self.get_saxs_scalar(a["saxs_kind"], qmin=a["qmin"], qmax=a["qmax"], qstar=a["qstar"], load=True)
        dfUV = self.uv.get_Abs_wl(list(wl) if isinstance(wl, range) else wl)
        cols = [c for c in dfUV.columns if c.startswith("Abs")]
        tX = dfX["time"].to_numpy(float)
        yX = dfX["I0"].to_numpy(float)[:, None]
        tU = dfUV["time"].to_numpy(float)
        YU = dfUV[cols].to_numpy(float)
        sx = np.nanmax(np.abs(yX))
        su = np.nanmax(np.abs(YU)) * np.sqrt(len(cols)) / uv_weight
        return [dict(t=tX, Y=yX, cw=np.array([1 / sx]), affine=True, names=["I0"]),
                dict(t=tU, Y=YU, cw=np.full(len(cols), 1 / su), affine=False, names=cols)]

    def _eval_blocks(self, theta, n, blocks, fit_scale=True, baseline=False):
        c, s, scale, shift = self._unpack_theta(theta, n, fit_scale)
        out = []
        for b in blocks:
            cc, ss = (scale * c + shift, scale * s) if b["affine"] else (c, s)
            G = self._gauss_design(b["t"], cc, ss, baseline)
            A = self._solve_amplitudes(G, b["Y"])
            d = G @ A - b["Y"]
            out.append((G, A, np.where(np.isfinite(d), d, 0.0) * b["cw"]))
        return out

    def _joint_residual(self, theta, n, blocks, fit_scale=True, baseline=False):
        return np.concatenate([r.ravel() for _, _, r in self._eval_blocks(theta, n, blocks, fit_scale, baseline)])

    def _theta_start_bounds(self, blocks, n, fit_scale=True, scale_bounds=(0.95, 1.05), shift_bounds=(-100.0, 100.0), sigma_bounds=None):
        tU = blocks[1]["t"]
        c0, s0 = self._guess_peaks(tU, np.nansum(blocks[1]["Y"], axis=1), n)
        dt = np.median(np.diff(tU))
        smin, smax = (dt / 2, tU.max() - tU.min()) if sigma_bounds is None else sigma_bounds
        s0 = np.clip(s0, smin * 1.01, smax * 0.99)
        x0 = [*c0, *np.log(s0), self._alignment.get("shift", 0.0)]
        lo = [*np.full(n, tU.min()), *np.full(n, np.log(smin)), shift_bounds[0]]
        hi = [*np.full(n, tU.max()), *np.full(n, np.log(smax)), shift_bounds[1]]
        if fit_scale:
            x0.append(self._alignment.get("scale", 1.0))
            lo.append(scale_bounds[0])
            hi.append(scale_bounds[1])
        return np.array(x0, float), (np.array(lo), np.array(hi))

    def _store_alignment(self, scale, shift, cost, success):
        self._alignment["scale"] = float(scale)
        self._alignment["shift"] = float(shift)
        self._alignment["refinement"] = {"method": "gauss_peak_joint", "scale": float(scale), "shift": float(shift), "cost": float(cost), "success": bool(success)}
        return self._alignment

#

