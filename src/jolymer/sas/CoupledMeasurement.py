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
                   min_seqi=0, max_seqi=100000, q_beamstop=0.006):
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
        m = SAXS_Measurement(path=SRC_DIR, filename=filename,
                             qmin=q_beamstop)
        m.time = time
        mlist.append(m)
        time += 2.1
        if i>=max_seqi:
            break
        if i<min_seqi:
            continue
    print("got number of ms:", len(mlist))
    out =  Ms(mlist)
    out.name = file_prefix
    return out

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
            df = df.iloc[1:]  # your skip logic
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
                                          bounds=((10, 0), (20, 0.1)))
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
        shift_bounds=(-1000.0, 1000.0),
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
