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

    def get_saxs_scalar(self, kind="I0", qstar=0.1,
                     qmin=0.04, qmax=0.4, load=False,
                        autorg_filename="autorg.dat"):
        """
        Returns frame, time, I0 (+ optionally errI0) for the saxs series.
        New: kind="autorg" sources I0/Rg from ATSAS `autorg` output
        (`self._autorg_df`, built by `load_autorg`) instead of the internal
        Guinier-region fit in `self.saxs_list.get_rgs(...)`. This is much
        faster (autorg already ran once from the GUI) and matches whatever
        autorg itself decided the Guinier region was, rather than the fixed
        (qmin, qmax) bounds / bounded fit used by kind="I0".
        """
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

    def align_uv_to_saxs(
        self,
        saxs_kind="I0",
        qstar=0.1,
        qmin=0.04, qmax=0.4,
        allow_fallback=True,
        load=True,
        min_t_saxs=0,
        max_t_saxs=np.inf,
        I0_min=0,
        scale0shift0=[2.2, -1400],
    ):
        dfUV = self.uv.get_scaled_Abs(refwl=self.uv.refwl,
                                      outwl=self.uv.outwl,
                                      alignment_time=self.uv.alignment_time)
        self.uv.df = dfUV
        dfSAXS = self.get_saxs_scalar(kind=saxs_kind, qstar=qstar, qmin=qmin, qmax=qmax,
                                      load=load)
        alignment = self.align_uv_df_to_saxs_df(
            dfUV,
            dfSAXS,
            min_t_saxs=min_t_saxs,
            max_t_saxs=max_t_saxs,
            I0_min=I0_min,
            preprocess_uv_mode="identity",
            scale0shift0=scale0shift0,
        )
        self._alignment = alignment
        self._alignment["saxs_kind"] = saxs_kind
        self._alignment["qstar"] = qstar
        self._alignment["qmin"] = qmin
        self._alignment["qmax"] = qmax
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

    def load_autorg(self, autorg_filename="autorg.dat", directory=None):
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
        autorg_name="autorg.dat",
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
    def align_uv_df_to_saxs_df(
        dfUV,
        dfSAXS,
        uv_col="Abs",
        saxs_col="I0",
        min_t_saxs=0,
        max_t_saxs=np.inf,
        I0_min=0,
        preprocess_uv_mode="identity",
        preprocess_uv_kwargs=None,
        scale_bounds=(0.99, 1.01),
        shift_bounds=(-3000.0, 3000.0),
        scale0shift0=[2.2, -1400],
    ):
        """
        Affine-align UV dataframe to SAXS dataframe in time.
        Returns:
            alignment dict with scale, shift, and score
        """
        from scipy.optimize import minimize
        from scipy.optimize import differential_evolution
        if preprocess_uv_kwargs is None:
            preprocess_uv_kwargs = {}
        # --- preprocess UV (hook for later deconvolution etc.)
        dfUVp = CoupledMeasurement.preprocess_uv(dfUV, mode=preprocess_uv_mode,
                               **preprocess_uv_kwargs)
        # --- extract arrays (but keep DF as source of truth)
        t_saxs = dfSAXS["time"].to_numpy()
        y_saxs = dfSAXS[saxs_col].to_numpy()
        mask = (
            (t_saxs >= min_t_saxs) &
            (t_saxs <= max_t_saxs) &
            (y_saxs >= I0_min)
        )
        t_saxs_m    = t_saxs[mask]
        y_saxs_m   = y_saxs[mask]
        # w = y_saxs  / (np.abs(dfSAXS[f'err{saxs_col}'].to_numpy()**2))
        # w = w  / np.sum(w)
        t_uv = dfUVp["time"].to_numpy()
        y_uv = dfUVp[uv_col].to_numpy()
        # normalize SAXS once
        y_saxs_n = (y_saxs_m - y_saxs_m.mean()) / y_saxs_m.std()
        def loss(p):
            a, b = p
            t_uv_warped = a * t_uv + b
            f = interp1d(
                t_uv_warped,
                y_uv,
                kind="cubic",
                # bounds_error=True,
                bounds_error=False,
                fill_value=0,
            )
            y_uv_i = f(t_saxs_m)
            # normalize UV after interpolation
            if np.std(y_uv_i) == 0:
                return np.inf
            y_uv_n = (y_uv_i - y_uv_i.mean()) / y_uv_i.std()
            # negative correlation (we minimize)
            return -np.sum(y_saxs_n * y_uv_n)
        # initial guess: no scaling, no shift
        # res = minimize(
        #     loss,
        #     x0=scale0shift0,
        #     bounds=[scale_bounds, shift_bounds],
        #     method="L-BFGS-B",
        # )
        res = differential_evolution(loss,
            x0=scale0shift0,
            bounds=[scale_bounds, shift_bounds])
        print("t_saxs:", t_saxs[:5], "...", t_saxs[-5:])
        print("t_uv:", t_uv[:5], "...", t_uv[-5:])
        a, b = res.x
        score = -res.fun
        return {
            "method": "affine",
            "scale": a,
            "shift": b,
            "score": score,
            "success": res.success,
        }

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
