import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Optional
import mdtraj as md
from scipy.constants import Avogadro as NA
from scipy import constants
import periodictable as pt

# based on typical mass density 1.7 g/cm**3 and wl=1.55 A-1

@dataclass
class bioMOLECULE:

    name: str = 'A'
    sequence: str = 'A'
    chemical_formula: str = '' # 'C10H14N5O7P'
    volume: float = 162 # A**3
    SLD: float = 1e-5
    dna: bool = False
    pdb_path: str = 'None'
    traj_path: str = 'None'
    volume: float=None
    DorP: str='aa'

    def get_c2I0_prefactor(self, wl=1.23,
                         solvent_SLD=9.47e-06,
                         **kwargs):
        sld = self.get_SLD(wl=wl, **kwargs)
        Drho = (sld - solvent_SLD) ## A-2
        Drho = Drho * 1e16 ## cm-2
        NA = constants.Avogadro
        V = self.get_volume() ## cm-3
        MW = self.get_MW() ## g
        # df_X['c'] = df_X.I0  / (MW*NA * Drho**2)
        # dfUV['c'] = dfUV.Abs * MW / (epsilon)
        # print('UV max c', dfUV.c.max())
        # print('X-ray max c', df_X.c.max())
        out = (Drho**2 * V**2 * NA) /1000 / MW
        print('volume:', V)
        print('MW:', MW)
        print('sld', sld)
        print('Drho', Drho)
        print('I0 prefactor', out)
        return out
        ## dfUV['I0'] = dfUV.c  / get_I0_prefactor
        ## df_X['c'] = I0 * get_I0_prefactor

    def get_epsilon(self):
        print("get_epsilon not implemented")
        epsilon = 1
        return epsilon

    def get_volume(self, radius=25, radii=None,
                   volume_per_water=30.13):
        if not self.volume is None:
            return self.volume
        import MDAnalysis as mda
        import numpy as np
        u = mda.Universe(self.pdb_path, self.traj_path)
        dna = u.select_atoms("resname DA DT DC DG")
        water = u.select_atoms("resname SOL and name OW")  # select water oxygens
        center = dna.center_of_mass()
        if radii is None:
            radii = [radius-5, radius, radius+5]
        dna_effective_volume = 0
        for radius in radii:
            water_within = water.select_atoms(f"point {center[0]} {center[1]} {center[2]} {radius}")
            num_water = len(water_within)
            water_volume = num_water * 18 * 1e24 / (0.992*NA) # 18 g/mol / 1g/L
            water_volume = num_water * volume_per_water # water volume from column.ipynb: 30.13
            sphere_volume = (4.0/3.0) * np.pi * (radius ** 3)
            dna_effective_volume += sphere_volume - water_volume
        dna_effective_volume = dna_effective_volume / len(radii)
        # print("Number of water molecules in sphere:", num_water)
        # print("water volume (Å³):", water_volume)
        # print("Sphere volume (Å³):", sphere_volume)
        # print("Total water volume (Å³):", water_volume)
        # print("Effective DNA volume (Å³):", dna_effective_volume)
        return dna_effective_volume

    def get_MW(self):
        if self.Mw is not None:
            print("Mw already loaded")
            return self.Mw
        pdb = md.load(self.pdb_path)
        mw = sum(atom.element.mass for atom in pdb.topology.atoms)
        self.Mw = mw
        return mw

    def get_density(self):
        MW = self.get_MW()
        density = MW / (self.volume * constants.Avogadro)
        return density
    # def get_density(self, **kwargs):
    #     """
    #     get volume is in A^3, mass is in g/mol
    #     return is in g/cm^3.0
    #     """
    #     volume = self.get_volume(**kwargs) * 1e-24
    #     volume = self.volume
    #     formula = pt.formula(f"dna: {self.sequence}")
    #     mass = formula.molecular_mass
    #     print('vm', volume, mass)
    #     density = mass / volume
    #     return density

    def get_SLD(self, wl=1.24):
        f = pt.formula(f"{self.DorP}: {self.sequence}")
        f.density = self.get_density()
        sld, _ = pt.xsf.xray_sld(f, wavelength=wl)
        return sld * 1e-6
    def get_delta_rho(self, wl=1.24):
        return (self.get_SLD(wl=wl) - self.solvent_SLD)
    def get_I0_prefactor(self, wl=1.24):
        MW = self.get_MW()
        V = self.volume
        Drho = self.get_delta_rho(wl=wl) * 1e16
        NA = constants.Avogadro
        prefactor = (Drho**2 * V**2 * NA) / (1000 * MW)
        return prefactor
    def c_from_I0(self, I0, wl=1.24):
        return I0 / self.get_I0_prefactor(wl=wl)
    def I0_from_c(self, c, wl=1.24):
        return c * self.get_I0_prefactor(wl=wl)
    def c_from_UV(self, A):
        MW = self.get_MW()
        return A * MW / self.epsilon


class DNA(bioMOLECULE):
    UV_ext = {
                # ϵ260	ϵmax	λmax	Nb	Confc	ϵ260	fd
        'A':	[15.02,	15.04,	259,	8,	0.5,	15.34,	1.0210],
        'dA':	[15.06,	15.08,	259,	8,	0.3,	15.34,	1.0183],
        'C':	[7.07,	8.74,	271,	5,	0.8,	7.60,	1.0753],
        'dC':	[7.10,	8.86,	271,	5,	0.3,	7.60,	1.0705],
        'G':	[12.08,	14.09,	252,	7,	0.3,	12.16,	1.0066],
        'dG':   [12.18,	14.23,	252,	4,	0.4,	12.16,	0.9987],
        'U':	[9.66,	9.78,	262,	5,	0.9,	1.021,	1.0570],
        'T':	[8.56,	9.49,	267,	5,	0.7,	8.70,	1.0165]
        }
    name: str = 'A'
    sequence: str = 'A'
    chemical_formula: str = '' # 'C10H14N5O7P'
    volume: float = 162 # A**3
    SLD: float = 1e-5
    dna: bool = False
    pdb_path: str = 'None'
    traj_path: str = 'None'

    def get_epsilon(self):
        epsilon = 0
        for letter in self.sequence:
            epsilon += self.UV_ext[letter][0]
            # print(letter, epsilon)
        if self.dna:
            h = 0
            for s in self.sequence:
                if s in ['G', 'C']:
                    h += 0.059
                elif s in ['A', 'T']:
                    h += 0.287
            h = h/len(self.sequence)
            print('h', h)
            epsilon = epsilon * (1-h)
        return epsilon

    def get_volume(self, radius=25, radii=None,
                   volume_per_water=30.13):
        return self.volume
        import MDAnalysis as mda
        import numpy as np
        u = mda.Universe(self.pdb_path, self.traj_path)
        dna = u.select_atoms("resname DA DT DC DG")
        water = u.select_atoms("resname SOL and name OW")  # select water oxygens
        center = dna.center_of_mass()
        if radii is None:
            radii = [radius-5, radius, radius+5]
        dna_effective_volume = 0
        for radius in radii:
            water_within = water.select_atoms(f"point {center[0]} {center[1]} {center[2]} {radius}")
            num_water = len(water_within)
            water_volume = num_water * 18 * 1e24 / (0.992*NA) # 18 g/mol / 1g/L
            water_volume = num_water * volume_per_water # water volume from column.ipynb: 30.13
            sphere_volume = (4.0/3.0) * np.pi * (radius ** 3)
            dna_effective_volume += sphere_volume - water_volume
        dna_effective_volume = dna_effective_volume / len(radii)
        # print("Number of water molecules in sphere:", num_water)
        # print("water volume (Å³):", water_volume)
        # print("Sphere volume (Å³):", sphere_volume)
        # print("Total water volume (Å³):", water_volume)
        # print("Effective DNA volume (Å³):", dna_effective_volume)
        return dna_effective_volume

    def get_MW(self):
        pdb = md.load_pdb(f"pymol/{self.name}.pdb")
        mw = sum(atom.element.mass for atom in pdb.topology.atoms)

        print(csample.name, f"Molecular weight: {mw:.3f} Da")
        MW = mw
        return MW
        formula = pt.formula(f"dna: {self.sequence}")
        return formula.molecular_mass

    def get_density(self, **kwargs):
        """
        get volume is in A^3, mass is in g/mol
        return is in g/cm^3.0
        """
        volume = self.get_volume(**kwargs) * 1e-24
        volume = self.volume
        formula = pt.formula(f"dna: {self.sequence}")
        mass = formula.molecular_mass
        print('vm', volume, mass)
        density = mass / volume
        return density

    def get_SLD(self, wl=1.23, **kwargs):
        f = pt.formula(f'dna: {self.sequence}')
        f.density = self.get_density(**kwargs)
        sld, _ = pt.xsf.xray_sld(f, wavelength=wl)
        return sld*1e-6

class PROTEIN(bioMOLECULE):
    name: str = 'A'
    sequence: str = 'A'
    chemical_formula: str = '' # 'C10H14N5O7P'
    volume: float = 162 # A**3
    SLD: float = 1e-5
    dna: bool = False
    pdb_path: str = 'None'
    traj_path: str = 'None'
    vbar: float = 0.733 # cm^3/g

    def estimate_volume_from_vbar(self):
        MW = self.get_MW()
        self.volume = MW * self.vbar / constants.Avogadro
        return self.volume


from importlib.resources import files
DATA = files("jolymer.data")

def load_bsa():
    bsa = PROTEIN(name="BSA",
                  pdb_path=str(DATA / "bsa.gro"),
                  traj_path=str(DATA / "bsa.xtc"),
                  epsilon=43824,
                  Mw=66430,
                  vbar=0.733, )
    bsa.estimate_volume_from_vbar()
    return bsa


# https://molbiotools.com/dnacalculator.php

def load_ac6():
    ac6 = DNA(name='AC6',
              pdb_path=str(DATA / "ac6.gro"),
              traj_path=str(DATA / "ac6.gro"),
            sequence='ACACACACACAC',
            SLD=1.31e-05)
    ac6.SLD17 = 1.5e-05
    ac6.epsilon_molbio = 119200
    ac6.volume = 2.81345e-21
    ac6.volume = 2.41345e-21
    ac6.volume = 2.9e-21 * 0.8
    ac6.epsilon = 13254 * 1.0
    return ac6
def load_gt6():
    gt6 = DNA(name='GT6', SLD=1.26e-05,
              pdb_path=str(DATA / "gt6.gro"),
              traj_path=str(DATA / "gt6.gro"),
            sequence='GTGTGTGTGTGT')
    gt6.epsilon_molbio = 114000
    gt6.SLD17 = 1.5e-05
    gt6.volume = 2.99177e-21
    gt6.volume = 3.99177e-21
    gt6.volume = 3.0e-21 * 1.2
    gt6.epsilon = 12384 * 1.4
    return gt6
def load_ac6gt6():
    ac6gt6 = DNA(name='AC6GT6', SLD=1.66e-05,
              pdb_path=str(DATA / "ac6gt6.gro"),
              traj_path=str(DATA / "ac6gt6.gro"),
            sequence='ACACACACACACGTGTGTGTGTGT',
            dna=True)
    ac6gt6.epsilon_molbio = 199048 * 0.6
    ac6gt6.SLD17 = 1.5e-05
    ac6gt6.volume = 5.83022e-21
    ac6gt6.volume = 7.33022e-21
    ac6gt6.volume = 6.0e-21 * 1.2
    ac6gt6.epsilon = 21202 * 1.0
    return ac6gt6
h2o = {'sld': 9.47e-06}
# csamples = [ac6, gt6, ac6gt6]

# for csample in csamples:
#     pdb = md.load_pdb(f"/home/johannes/jophd/reports/odna/pymol/{csample.name}.pdb")
#     mw = sum(atom.element.mass for atom in pdb.topology.atoms)

#     print(csample.name, f"Molecular weight: {mw:.3f} Da")
#     csample.Mw = mw

# 	This work	l/mmol/cm			Previous https://academic.oup.com/nar/article/32/1/e13/1195443


