from dataclasses import dataclass
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from ..Measurement import Measurement
from .. import os_utility as osu
from ..samples.bioMOLECULE import bioMOLECULE, load_ac6
from ..import plot_utility as plu

ac6 = load_ac6()


@dataclass
class onlineUV(Measurement):

    instrument: str = 'no instrument'
    path: str = '/home/johannes/pCloudDrive/jophd/rawdata/onlineUV/20221117/'
    spec_filename: str = '20221117_AC6_0_1_34ul_c01_000001.spc'
    name: str = 'unnamed'
    short_name: str = 'unnamed'
    color: str = None
    marker: str = None
    linestyle: str = ''
    shift: float = 1
    maxtime: float = 1000000
    mintime: float = 0
    refwl: float=280
    outwl: float=280
    alignment_time: float=280

    def get_log_filepath(self, spec_filename=None):
        if spec_filename is None:
            spec_filename = self.spec_filename
        log_filename = f'{spec_filename.split('.spc')[0]}.log'
        return Path(path) / log_filename

    def load_log(self, spec_filename=None):
        log_filepath = self.get_log_filepath(spec_filename=spec_filename)

    def __post_init__(self, sample:bioMOLECULE=ac6):
        # Example: "20221117_AC6_0_1_34ul_c01_000001.spc"
        self.name = self.spec_filename.split("01")[0]
        self.filename = f"{self.name}_onlineUV.dat"
        self.sample = sample
        osu.create_path("ellution")

    def get_spec_filename(self):
        return Path(self.path) / self.spec_filename

    def get_filename(self):
        return Path("ellution") / self.filename

    def get_data(self):
        """Convert .spc → .dat, then read the data."""
        import spc_spectra as spc

        spc_file = spc.File(self.get_spec_filename())
        spc_file.write_file(self.get_filename())
        return self.load_data()

    def load_data(self):
        """Load the whitespace-delimited .dat file as a DataFrame."""

        # Read one line to count columns correctly
        with open(self.get_filename(), "r") as f:
            first_line = f.readline()
            n_cols = len(first_line.split())
        colnames = ["wl"] + [f"Abs{i}" for i in range(n_cols - 1)]
        df = pd.read_csv(
            self.get_filename(),
            sep=r"\s+",
            names=colnames,
            engine="python",
        )
        return df

    def get_scaled_Abs(
        self,
        refwl=None,
        outwl=None,
        alignment_time=None,
        show=False,
        show_wl=None,
    ):
        """
        Return scaled absorbance vs time for one or multiple wavelengths.
        Parameters
        ----------
        refwl : int, optional
            Reference wavelength. If unavailable, the closest measured wavelength
            is used.
        outwl : int or list[int], optional
            Wavelength(s) to scale to the reference wavelength.
        alignment_time : float, optional
            Time point used for scaling.
        show : bool
            Plot reference and scaled absorbance.
        show_wl : list[int], optional
            Additional wavelengths to plot.
        Returns
        -------
        pd.DataFrame
            Scaled absorbance:
                time, Abs, err_Abs
        Notes
        -----
        Scaling assumes:
            Abs_out_scaled = Abs_out * Abs_ref(t_align) / Abs_out(t_align)
        """
        if refwl is None:
            refwl = self.refwl
        if outwl is None:
            outwl = self.outwl
        if alignment_time is None:
            alignment_time = self.alignment_time
        if show_wl is None:
            show_wl = []
        if not isinstance(outwl, (list, tuple, np.ndarray)):
            outwl = [outwl]
        df = self.load_data().T
        uv_channels = {}
        for i in range(20, df.shape[1]):
            wl = int(df.iat[0, i])
            Abs = df.iloc[1:, i].astype(float).values
            time = np.arange(len(Abs)) * 0.946
            ddf = pd.DataFrame({
                "time": time,
                "Abs": Abs,
            })
            ddf = ddf[
                (ddf.time < self.maxtime) &
                (ddf.time > self.mintime)
            ].reset_index(drop=True)
            uv_channels[wl] = ddf
        available_wl = np.array(sorted(uv_channels.keys()))
        refwl_used = available_wl[
            np.argmin(np.abs(available_wl - refwl))
        ]
        ref = uv_channels[refwl_used]
        scaled_channels = []
        scale_factors = []
        for wl in outwl:
            if wl not in uv_channels:
                raise ValueError(
                    f"Wavelength {wl} nm not available. "
                    f"Available range: {available_wl.min()}-{available_wl.max()} nm"
                )
            channel = uv_channels[wl].copy()
            ref_abs = np.interp(
                alignment_time,
                ref.time,
                ref.Abs
            )
            out_abs = np.interp(
                alignment_time,
                channel.time,
                channel.Abs
            )
            if np.isclose(out_abs, 0):
                raise ValueError(
                    f"Cannot scale {wl} nm: absorbance is zero "
                    f"at alignment time {alignment_time}"
                )
            factor = ref_abs / out_abs
            channel["Abs"] *= factor
            scaled_channels.append(channel)
            scale_factors.append(factor)
        Abs_stack = np.vstack([
            x.Abs.values for x in scaled_channels
        ])
        scaled_out = pd.DataFrame({
            "time": scaled_channels[0].time.values,
            "Abs": Abs_stack.mean(axis=0),
            "err_Abs": Abs_stack.std(axis=0),
        })
        scaled_out.attrs.update({
            "refwl_requested": refwl,
            "refwl_used": refwl_used,
            "outwl": outwl,
            "scale_factors": scale_factors,
            "mean_scale_factor": float(np.mean(scale_factors)),
            "alignment_time": alignment_time,
        })
        if show:
            self.plot_scaled_Abs(scaled_out, ref)
        return scaled_out

    def plot_scaled_Abs(self, scaled_out, ref):
        fig, ax = plt.subplots()
        refwl = scaled_out.attrs["refwl_requested"]
        refwl_used = scaled_out.attrs["refwl_used"]
        outwl = scaled_out.attrs["outwl"]
        scale_factors = scaled_out.attrs["scale_factors"]
        mean_scale_factor = scaled_out.attrs["mean_scale_factor"]
        alignment_time = scaled_out.attrs["alignment_time"]
        ax.plot(
            scaled_out.time,
            scaled_out.Abs,
            label=(
                f"scaled {outwl} nm "
                f"(factor={np.mean(scale_factors):.3f})"
            )
        )
        ax.plot(
            ref.time,
            ref.Abs,
            label=f"{refwl_used} nm reference"
        )
        idx = np.argmin(
            np.abs(ref.time.values - alignment_time)
        )
        ax.scatter(
            ref.time.iloc[idx],
            ref.Abs.iloc[idx],
            label="alignment point"
        )
        ax.set_xlabel("time [s]")
        ax.set_ylabel("Abs")
        ax.legend(fontsize='xx-small')
        return ax

    def plot_full_wavelength_map(self,
                              wl_min=210,
                              wl_max=900,
                              refwl=None,
                                 ax=None,
                              alignment_time=None,
                              cmap="plasma",
                              show=True):
        """
        Plot full wavelength-resolved absorbance map (time vs wavelength).
        Uses reference wavelength scaling similar to get_scaled_Abs.
        """
        if refwl is None:
            refwl = self.refwl
        if alignment_time is None:
            alignment_time = self.alignment_time
        df = self.load_data().T
        wavelengths = []
        absorbance_matrix = []
        time_axis = None
        ref_trace = None
        ref_idx = None
        # --- collect data ---
        for i in range(20, 151):
            wl = int(df.iat[0, i])
            if wl < wl_min or wl > wl_max:
                continue
            Abs = df[i][1:].astype(float).values
            time = np.linspace(0, len(Abs), len(Abs)) * 0.946
            # store reference
            if np.isclose(wl, refwl, atol=0.5):
                ref_trace = Abs.copy()
                ref_time = time.copy()
                ref_idx = np.argmin(np.abs(ref_time - alignment_time))
            wavelengths.append(wl)
            absorbance_matrix.append(Abs)
            if time_axis is None:
                time_axis = time
        wavelengths = np.array(wavelengths)
        absorbance_matrix = np.array(absorbance_matrix)
        # --- scaling using reference wavelength ---
        if ref_trace is None:
            raise ValueError("Reference wavelength not found.")
        ref_value = ref_trace[ref_idx]
        scaled_matrix = np.zeros_like(absorbance_matrix)
        for j in range(len(wavelengths)):
            wl_trace = absorbance_matrix[j]
            factor = ref_value / (wl_trace[ref_idx] + 1e-12)
            scaled_matrix[j] = wl_trace * factor
        # --- plot ---
        if show:
            if ax is None:
                fig, ax = plt.subplots(figsize=(4,2))
            else:
                ax=ax
                fig = ax.get_figure()
            extent = [
                time_axis[0],
                time_axis[-1],
                wavelengths.min(),
                wavelengths.max()
            ]
            # sort by wavelength for nicer gradient
            sort_idx = np.argsort(wavelengths)
            wavelengths_sorted = wavelengths[sort_idx]
            matrix_sorted = scaled_matrix[sort_idx]
            # colors = plu.cm_for_l('plasma', wavelengths)
            norm = plt.Normalize(wl_min, wl_max)
            comap = cm.get_cmap(cmap)
            for wl, absorbance in zip(wavelengths_sorted, matrix_sorted):
                color = comap(norm(wl))
                ax.plot(time_axis, absorbance, color=color)
            ax.set_xlabel("Time [s]")
            ax.set_ylabel("Abs [A.U.]")
            # cbar = plt.colorbar(im, ax=ax)
            # cbar.set_label("Scaled Absorbance")
            sm = cm.ScalarMappable(norm=norm, cmap=comap)
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax)
            cbar.set_label("Wavelength (nm)")
        return {
            "time": time_axis,
            "wavelengths": wavelengths,
            "abs_matrix": scaled_matrix
        }

    def to_concentration(self, biomolecule, wl, pathlength_cm=1.0):
        """
        Returns dataframe with:
          time [s]
          c_g_per_L
          c_mol_per_L
          I0_per_cm
        """
        Abs = df.Abs.values
        eps = biomolecule.get_epsilon(wl)   # L / (mol·cm)
        Mw  = biomolecule.get_Mw()           # g / mol
        df['cmol'] = Abs / (eps * pathlength_cm)
        df['c']   = c_mol_L * Mw
        return df

    def to_I0(self, biomolecule, wl, pathlength_cm=1.0):
        df = self.to_concentration(self, biomolecule, wl, pathlength_cm=pathlength_cm)
        df['pseudoI0'] = df.c / biomolecule.get_I0_prefactor()
        return df
