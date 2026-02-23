import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.signal import correlate
from pathlib import Path
from dataclasses import dataclass

from .SAXS_Measurement import SAXS_Measurement
from .ms import Ms
from jolymer.samples.bioMOLECULE import *
from jolymer.uv.onlineUV import onlineUV

def ms_from_folder(path='workdirs/ACGT6c80_waxs', file_prefix='waxs',
                   min_seqi=0, max_seqi=100000, q_beamstop=0.006,
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
    time=0
    for i, filename in enumerate(sorted_filenames):
        if i>max_seqi:
            break
        if i<min_seqi:
            continue
        m = SAXS_Measurement(path=SRC_DIR, filename=filename,
                             qmin=q_beamstop, angular_unit=angular_unit)
        m.time = time
        mlist.append(m)
        time += 2.1
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
        qmin = qmin if qmin is not None else qmin_auto
        qmax = qmax if qmax is not None else qmax_auto
        mask_w = (q_w > qmin) & (q_w < qmax)
        if mask_w.sum() < 5:
            raise ValueError("Overlap region too small")
        outdict['num_frames'] = I_w.shape[1]
        outdict['qw_ov'] = q_w[mask_w]
        outdict['Iw_ov'] = I_w[mask_w, :]
        interp = interp1d(q_s, I_s, axis=0, bounds_error=False,
                          fill_value=np.nan)
        Is_ov = interp(outdict['qw_ov'])
        sigma_w_ov = sigma_w[mask_w, :]
        interpw = interp1d(q_s, sigma_s, axis=0, bounds_error=False,
                           fill_value=np.nan)
        sigma_s_ov = interpw(outdict['qw_ov'])
        outdict['Is_ov'] = Is_ov
        outdict['sigma_s_ov'] = sigma_s_ov
        outdict['sigma_w_ov'] = sigma_w_ov
        return outdict

    def align_saxs_waxs_global(
        self,
        qmin=None,
        qmax=None,
        use_errors=True,
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
        # def loss(params):
        #     a = params[0]
        #     offsets = params[1:]
        #     model = a * matrix_dict['Is_ov'] + offsets
        #     total_sigma = np.sqrt(matrix_dict['sigma_w_ov']**2 +\
        #             (a * matrix_dict['sigma_s_ov'])**2)
        #     weighted_mae = np.mean(np.abs(matrix_dict['Iw_ov'] - model) / total_sigma)
        #     penalty = 0
        #     if np.any(model < 0):
        #         penalty = np.sum(np.abs(model[model < 0])) * 1e6
        #     return weighted_mae + penalty
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
            # Optional hard constraint: no negative modeled intensity
            # if np.any(model < 0):
            #     chi2 += 1e12
            return chi2
        mean_int = np.mean(matrix_dict['Iw_ov'])
        bounds = [(0.1, 10.0)] + [(-20*mean_int, 20*mean_int)] *\
                matrix_dict['num_frames']
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
        n_bins=200,
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
        for t, m in enumerate(self.ms_saxs.ms):
            path = Path(m.path)
            filename = f"merge.{m.filename}"
            outdf = pd.DataFrame({
                "q": q_final,
                "I": I_final[:, t],
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
        # --- rebuild overlap (guarantees consistency) ---
        md = self.interpolate_matrix(qmin=qmin, qmax=qmax)
        Iw = md['Iw_ov']
        Is = md['Is_ov']
        sw = md['sigma_w_ov']
        ss = md['sigma_s_ov']
        n_q, n_frames = Iw.shape
        chi2_per_frame = np.zeros(n_frames)
        for i in range(n_frames):
            # model = a * Is[:, i] + b.mean()
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
    saxs_kind: str = 'I0'

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
                        qmin=0.04, qmax=0.4, load=False):
        """
        Returns frame, time, I0 for the saxs
        """
        q = self._q if self._q is not None else self.build_saxs_matrix()[0]
        if kind == "I0":
            frames = []
            times = []
            rgdf = self.saxs_list.get_rgs(qmin=qmin, qmax=qmax,
                                          load=load, I00=0.01, Rg0=12,
                                          bounds=((10, 0), (20, 1.0)))
            for frame, m in enumerate(self.saxs_list):
                frames.append(frame+1)
                times.append(m.time)
            return pd.DataFrame({
                'frame': frames,
                'time': times,
                'I0': rgdf.I0,
                'errI0': rgdf.err_I0})
        elif kind == "Iq":
            q = self._q if self._q is not None else self.build_saxs_matrix()[0]
            idx = np.argmin(np.abs(q - qstar))
            print(idx)
            frames = []
            times = []
            I0s = []
            errI0s = []
            for frame, m in enumerate(self.saxs_list):
                frames.append(frame+1)
                times.append(m.time)
                I0s.append(float(m.get_data()['I'][idx]))
                errI0s.append(float(m.get_data()['err_I'][idx]))
            return pd.DataFrame({
                'frame': frames,
                'time': times,
                'I0': I0s,
                'errI0': errI0s})
                # 'I0': self._I[idx, :]})
        # elif kind == 'integrate':
        #     idx_min = np.searchsorted(q, qmin)
        #     idx_max = np.searchsorted(q, qmax)
        #     print(np.trapz(self._I[idx_min:idx_max], q[idx_min:idx_max], axis=0))
        #     return np.trapz(self._I[idx_min:idx_max], q[idx_min:idx_max],
        #                     axis=0)
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
            preprocess_uv_mode="identity",  # later: "custom"
            scale0shift0=scale0shift0,
        )
        self._alignment = alignment
        self._alignment["saxs_kind"] = saxs_kind
        self._alignment["qstar"] = qstar
        self._alignment["qmin"] = qmin
        self._alignment["qmax"] = qmax
        return self._alignment

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
        )

        saxs_time = np.array([m.time for m in self.saxs_list])
        uv_on_saxs = interp(saxs_time)

        self._uv_on_saxs = pd.DataFrame({
            "frame": self._x,
            "time": saxs_time,
            "Abs": uv_on_saxs,
        })

        return self._uv_on_saxs

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
        scale_bounds=(0.5, 5.0),
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
                kind="linear",
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
