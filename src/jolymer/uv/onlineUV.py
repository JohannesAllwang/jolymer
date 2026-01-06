from dataclasses import dataclass
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ..Measurement import Measurement
from .. import os_utility as osu
from ..samples.bioMOLECULE import bioMOLECULE, ac6


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

    # --------------------------------------------------------
    def __post_init__(self, sample:bioMOLECULE=ac6):
        # Example: "20221117_AC6_0_1_34ul_c01_000001.spc"
        self.name = self.spec_filename.split("01")[0]
        self.filename = f"{self.name}_onlineUV.dat"
        self.sample = sample
        osu.create_path("ellution")

    # --------------------------------------------------------
    def get_spec_filename(self):
        return Path(self.path) / self.spec_filename

    def get_filename(self):
        return Path("ellution") / self.filename

    # --------------------------------------------------------
    def get_data(self):
        """Convert .spc → .dat, then read the data."""
        import spc_spectra as spc

        spc_file = spc.File(self.get_spec_filename())
        spc_file.write_file(self.get_filename())
        return self.load_data()

    # --------------------------------------------------------
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

    # --------------------------------------------------------
    def get_scaled_Abs(self, refwl=None, outwl=None, alignment_time=None,
                       show=False):
        """
        Return scaled Abs vs time for two wavelengths.
        refwl: reference wl (typically 260 nm)
        outwl: wavelength to scale (e.g. 320 nm)
        """
        if refwl is None:
            refwl = self.refwl
        if outwl is None:
            outwl = self.outwl
        if alignment_time is None:
            alignment_time = self.alignment_time

        df = self.load_data().T
        # df.T[i][0] = wavelength
        # df.T[i][1:] = absorption series

        outdict = {}

        for i in range(20, 151):   # if you intentionally meant two indices use [20,150]
            # print(df.iloc[i])
            wl = int(df.iat[0, i])
            # print(wl)
            Abs = df[i][1:].astype(float).values
            time = np.linspace(0, len(Abs), len(Abs)) * 0.946
            ddf = pd.DataFrame({
                "time": time,
                "Abs": Abs,
            })
            ddf = ddf[ddf.time<self.maxtime]
            ddf = ddf[ddf.time>self.mintime]
            if wl == int(refwl):
                outdict["ref"] = ddf
                if refwl == outwl:
                    outdict["scaled_out"] = ddf
                    outdict["scale_factor"] = 1
            elif wl == int(outwl):
                if "ref" not in outdict:
                    raise ValueError("Reference wavelength not processed yet.")
                ddf_ref = outdict["ref"]
                # scale at alignment_time index
                factor = ddf_ref.Abs.iloc[alignment_time] / ddf.Abs.iloc[alignment_time]
                ddf["Abs"] = ddf["Abs"] * factor
                outdict["scaled_out"] = ddf
                outdict["scale_factor"] = factor
        if show:
            fig, ax = plt.subplots()
            ax.plot(outdict["scaled_out"].time, outdict["scaled_out"].Abs, label=f"$Abs_{{{outwl}}} \\times {factor:.2f}$")
            ax.plot(outdict["ref"].time, outdict["ref"].Abs, label=f"$Abs_{{{refwl}}}$ reference")
            ax.plot(outdict["ref"].time.iloc[alignment_time], outdict["ref"].Abs.iloc[alignment_time], label=f"$Abs_{{{refwl}}}$ reference")
            ax.set_xlabel(f"time [s]")
            ax.set_ylabel(f"Abs")
            ax.legend()

        return outdict["scaled_out"]
        # return pd.DataFrame(outdict)

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
