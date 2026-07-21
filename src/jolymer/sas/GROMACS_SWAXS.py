from scipy import optimize
import numpy as np
import re
from os.path import join
from matplotlib.colors import LogNorm
from pylab import cm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd

import shutil
import os
from tqdm.auto import tqdm

import sasmodels
from sasmodels import data as sasmodels_data
import pyFAI
import fabio
from pathlib import Path

import MDAnalysis as mda

from .. import os_utility as osu
from .SAXS_Measurement import SAXS_Measurement
from .. import jocolors

from dataclasses import dataclass

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

@dataclass
class GROMACS_SWAXS(SAXS_Measurement):

    instrument = 'gromacs swaxs'

    rawpath: str= ''
    path: str = '~'
    filename: str='waxs_final.xvg'
    name: str='unnamed'
    angular_unit: str='nm'
    gromacs_path: str='not defined'
    maxq: float=np.inf
    color: str=None
    marker: str='o'
    linestyle: str='-'
    mdpath: str=None
    md_basename: str=None
    spectra_filename: str=None
    proj_filename: str=None
    rmsf_filename: str=None
    rms_filename: str=None
    rg_filename: str=None
    eigenval_filename: str=None
    swaxspot_filename: str=None
    npt_filename: str='npt.gro'
    # npt_filename: str='empty.tpr'
    xtc_filename: str=None
    u: mda.Universe=None

    plot_standardwidth: float=3.25
    plot_standardheight: float=3.25
    plot_shortheight: float=2
    plot_longheight: float=4

    min_index: int=0

    def __post_init__(self):
        if self.spectra_filename is None:
            self.spectra_filename = self.md_basename+'_spectra.xvg'
        if self.proj_filename is None:
            self.proj_filename = 'proj_'+self.md_basename+'.xvg'
        if self.rmsf_filename is None:
            self.rmsf_filename = 'rmsf_'+self.md_basename+'.xvg'
        if self.rms_filename is None:
            self.rms_filename = 'rms_'+self.md_basename+'.xvg'
        if self.rg_filename is None:
            self.rg_filename = 'rg_'+self.md_basename+'.xvg'
        if self.swaxspot_filename is None:
            self.swaxspot_filename = self.md_basename+'_pot.xvg'
        if self.xtc_filename is None:
            self.xtc_filename = self.md_basename+'_center.xtc'
            self.xtc_filename = self.md_basename+'_fit.xtc'
        # if self.eigenval_filename is None:
        #     self.eigenval_filename = self.md_basename+'_spectra.xvg'

    def get_eigenval(self, mdpath=None, **kwargs):
        if mdpath is None:
            mdpath = self.mdpath
        # if eigenval_filename is None:
        eigenval_filename = self.eigenval_filename
        if eigenval_filename is None:
            return None
        # print(mdpath, eigenval_filename)
        filename = Path(mdpath, eigenval_filename)
        with open(filename) as file:
            out = {}
            i = 1
            data = []
            for line in file:
                if not line.startswith(('@', '#', '&')):
                    data.append([float(x) for x in line.split()])
        df = pd.DataFrame({'pca': [int(row[0]) for row in data],
                           'RMSF': [row[1] for row in data]}).set_index('pca')
        return df

    def get_swaxspot(self, path=None, filename=None):
        if path is None:
            path = self.mdpath
        if filename is None:
            filename = self.swaxspot_filename
        with open(Path(path,filename)) as f:
            outdict = {'time': []}
            qs = []
            for line in f:
                if len(line.split('legend')) > 1:
                    qeqq = line.split('legend')[1]
                    q = float(line.split('=')[1].replace('"', '')) / 10
                    qs.append(q)
                    outdict[q] = []
                elif len(line.split()) > 3 and not line[0] in ['#', '@']:
                    numbers = [float(n) for n in line.split()]
                    outdict['time'].append(numbers.pop(0))
                    for number, q in zip(numbers, qs):
                        outdict[q].append(number)
        df = pd.DataFrame(outdict)
        # df.q = df.q/10
        df.time = df.time / 1000
        return df

    def get_filename(self):
        return join(self.path, self.filename)

    def get_spectra_filename(self):
        return join(self.path, self.spectra_name)

    def get_rg(self, path=None, filename=None):
        if path is None:
            path = self.mdpath
        if filename is None:
            filename = self.rg_filename
        full_path = Path(path, filename)
        print('rg_filename', full_path)
        gdict = self.read_xvg(full_path, include1=False, every_n=1, max_out=np.inf,
                 names=['step', 'Rg', 'Rx', 'Ry', 'Rz'],
                 onlyfirst=False)
        df = gdict['dfs'][0]
        df['time'] = df.step/1000
        return df

    # def get_data(self, engine='pandas', **kwargs):
    #     """
    #     Reads the data file at self.path / self.filename
    #     Alternatively, a path can be provided as a keyword argument.
    #     """
    #     out =  self.read_xvg(self.get_filename(), names=['q', 'I', 'err_I'],
    #                          include1=True, onlyfirst=True)
    #     return out['dfs'][0]

    def get_gromacs(self, path=None, filename='waxs_final.xvg'):
        gromacs_path = path
        if path is None:
            gromacs_path = self.gromacs_path
        filename = join(gromacs_path, filename)
        with open(filename) as f:
            i = 0
            for line in f:
                if line == '@type xydy\n':
                    skiprows = i+1
                if line == '&\n':
                    print(line)
                    nrows = i - skiprows
                    break
                i += 1
        print('skiprows', skiprows)
        print('nrows', nrows)
        df = pd.read_csv(filename,
                 sep=r'\s+', names=['q', 'I', 'errI'], skiprows=skiprows,
                 nrows = nrows,
                 dtype=np.float64)
        return df

    def get_proj(self, mdpath=None, filename='proj.xvg', pc='pc', **kwargs):
        eigenval_df = self.get_eigenval(mdpath=mdpath)
        # N_eigenval = len(eigenval_df)
        if mdpath is None:
            mdpath = self.mdpath
        filename = Path(mdpath, filename)
        with open(filename) as file:
            out = {}
            i = 1
            data = []
            for line in file:
                if line.startswith(('&')):
                    df = pd.DataFrame({'x': [row[0] for row in data],
                                       'y': [row[1] for row in data]})
                    df['time'] = df.x/1000
                    out['time'] = df.time
                    # RMSFi = float(eigenval_df.loc[i].iloc[0])
                    out[f'{pc}{i}'] = df.y
                    data = []
                    i += 1
                if not line.startswith(('@', '#', '&')):
                    data.append([float(x) for x in line.split()])
        outdf = pd.DataFrame(out)
        return outdf

    def calc_fbs(self,
                selection=''):
        pass

    def calc_rhoR(self,
                  selection=''):
        pass

    def calc_Rg(self,
                selection=''):
        pass

    def get_NinR(
            self,
            selection="resname SOL",
            reference="name P",
            reference_index=5,
            R=30.0,
            step=1,
            center="atom",
        ):
            """
            Parameters
            ----------
            selection : str
            reference : str
                Atom selection defining the reference.
            reference_index : int
                Which atom to use if multiple atoms match.
            R : float
                Radius in Å.
            step : int
                Trajectory stride.
            center : {"atom", "com"}
                Use either a single atom or the center of mass as reference.
            Returns
            -------
            DataFrame
                time (ns), NinR
            """
            u = self.get_u()
            particles = u.select_atoms(selection)
            ref_atoms = u.select_atoms(reference)
            if center == "atom":
                ref = ref_atoms[reference_index]
            elif center == "com":
                ref = ref_atoms
            else:
                raise ValueError(...)
            times = []
            counts = []
            for ts in tqdm(u.trajectory[::step], desc="Counting waters"):
                if center == "atom":
                    ref_pos = ref.position
                else:
                    ref_pos = ref.center_of_mass()
                r = np.linalg.norm(
                    particles.positions - ref_pos,
                    axis=1,
                )
                counts.append(np.count_nonzero(r < R))
                times.append(ts.time / 1000)
            return pd.DataFrame(
                {
                    "time": times,
                    "N": counts,
                }
            )

    @staticmethod
    def pick_chi2(df, n=4, min_dt=5.0, min_time=0, ascending=True):
        # sort by chi2 ascending
        df_sorted = df.sort_values("chi2", ascending=ascending)
        chosen = []
        for _, row in df_sorted.iterrows():
            t = row["time"]
            if t < min_time:
                continue
            # only accept if time difference >= min_dt from all already chosen
            if all(abs(t - c["time"]) >= min_dt for c in chosen):
                chosen.append(row)
            if len(chosen) == n:
                break
        return pd.DataFrame(chosen).sort_values('time')

    def get_u(self, rebuild=False):
        if not self.u is None:
            if not rebuild:
                print('u was already loaded')
                return self.u
        from MDAnalysis.transformations import unwrap, center_in_box, wrap
        from MDAnalysis.analysis import align
        npt_path = Path(self.mdpath, self.npt_filename)
        xtc_path = Path(self.mdpath, self.xtc_filename)
        print('loading', xtc_path)
        u = mda.Universe(npt_path, xtc_path)
        nucleic = u.select_atoms("nucleic")
        workflow = [wrap(u.atoms, compound='residues')]
        u.trajectory.add_transformations(*workflow)
        # align.AlignTraj(u, u, select='nucleic', in_memory=True).run()
        self.u = u
        return u

    def get_medoid_path(self, time0=0, timef=1000, path=None):
        if path is None:
            outpath = self.mdpath
        out = Path(outpath, f"medoid_{self.md_basename}_{time0:.0f}ps{timef:.0f}ps.pdb")
        return out

    def get_Dmax(self, time0=0, timef=1000):
        from scipy.spatial.distance import pdist
        u = self.get_u()
        sel = u.select_atoms("nucleic and name P")
        outd = {'time': [],
                'Dmax': []}
        start_frame = int(time0 / u.trajectory.dt)
        end_frame   = int(timef / u.trajectory.dt)
        for ts in u.trajectory[start_frame:end_frame]:
            coords = sel.atoms.positions
            Dmax = pdist(coords).max()  # ✅ memory-efficient
            outd['time'].append(ts.time/1000)
            outd['Dmax'].append(Dmax)
        out = pd.DataFrame(outd)
        return out

    def save_medoid(self, time0, timef, selection="name P", outpath=None):
        """
        Compute the medoid frame (closest to all others) using a full distance matrix,
        optionally restricted to a selection (e.g., backbone or P atoms).
        Parameters
        ----------
        time0, timef : float
            Start and end times (ps) to consider for medoid selection.
        selection : str
            Atom selection string, default to "name P" (DNA phosphates).
        """
        import numpy as np
        from MDAnalysis.analysis import rms
        from pathlib import Path
        out = {}
        u = self.get_u()
        atoms = u.select_atoms(selection)
        # select frames in the given time window
        # frames = [ts.frame for ts in u.trajectory if time0 <= ts.time <= timef]
        first_frame_time = u.trajectory[0].time
        dt = u.trajectory.dt  # ps per frame
        start_idx = int((time0 - first_frame_time) / dt)
        end_idx   = int((timef - first_frame_time) / dt)
        frames = list(range(start_idx, end_idx + 1))
        M = len(frames)
        print(M)
        N = atoms.n_atoms
        # extract coordinates for selected frames
        coords = np.zeros((M, N, 3))
        for i, f in enumerate(frames):
            u.trajectory[f]
            coords[i] = atoms.positions.copy()
        # flatten coordinates for RMSD calculations
        coords_flat = coords.reshape(M, -1)  # shape: [M, N*3]
        # compute full pairwise RMSD matrix
        dist_matrix = np.zeros((M, M))
        for i in range(M):
            diff = coords_flat[i] - coords_flat  # broadcasted difference
            rmsd_i = np.sqrt(np.mean(diff**2, axis=1))
            dist_matrix[i] = rmsd_i
        # compute average distance of each frame to all others
        avg_dist = dist_matrix.mean(axis=1)
        # medoid = frame with minimal average distance
        best_idx_local = int(np.argmin(avg_dist))
        best_frame = frames[best_idx_local]
        print(f"Medoid (closest to all others): frame {best_frame}, RMSD={avg_dist[best_idx_local]:.3f} Å")
        # save PDB of medoid
        u.trajectory[best_frame]
        out_path = self.get_medoid_path(time0=time0, timef=timef, path=outpath)
        u.atoms.write(out_path)
        rmsd_to_medoid = dist_matrix[best_idx_local]  # shape (M,)
        rmsd_df = pd.DataFrame({'time': frames,
                                'rmsd': rmsd_to_medoid})
        rmsd_df.time = rmsd_df.time * dt / 1000
        out['rmsd_df'] = rmsd_df
        out['medoid_path'] = out_path
        return out

    def rmsd_full_traj_to_medoid(self, medoid_path, selection="nucleic"):
        """
        Compute RMSD of a full trajectory to a given medoid structure.
        Parameters
        ----------
        traj_universe : MDAnalysis.Universe
            Universe with trajectory loaded.
        medoid_path : str
            Path to saved medoid PDB.
        selection : str
            Atom selection (default: 'nucleic').
        dt : float
            Time step between frames in ps.
        """
        from MDAnalysis.analysis import rms
        traj_universe = self.get_u()
        dt = traj_universe.trajectory.dt
        print('mpath', medoid_path)
        ref = mda.Universe(medoid_path)  # medoid structure as reference
        R = rms.RMSD(
            traj_universe,
            ref,
            select=selection,
            groupselections=[
                "nucleic and name P"                                     # NMP
            ],
            center=True,
            superposition=True,
            verbose=True
        )
        step = 1
        R.run(step=step) #, stop=7500)
        rmsd_df = pd.DataFrame({
            "frame": R.rmsd[:,0].astype(int),
            "time": R.rmsd[:,1] / 1000,  # ns
            "rmsd": R.rmsd[:,2],               # Å
            "rmsdPP": R.rmsd[:,3]               # Å
        })
        return rmsd_df

    def _plot_bestfit_rg_segments(self, df, dfrg, medoid_times,
                                  window=2.0, ax=None):
        """
        Plot Rg traces from the best-fit intervals on a continuous
        artificial x-axis.
        Parameters
        ----------
        window : float
            Number of ns before each medoid to include.
        """
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        xoffset = 0.0
        for medoid_time, color in zip(
                medoid_times,
                [jocolors.tstum1,
                 jocolors.tstum5,
                 jocolors.tstum3,
                 jocolors.tstum7]):
            mask_swaxs = (
                (df.time >= medoid_time - window) &
                (df.time <= medoid_time)
            )
            seg_swaxs = df.loc[mask_swaxs].copy()
            if len(seg_swaxs) == 0:
                continue
            x_swaxs = (
                seg_swaxs.time -
                seg_swaxs.time.iloc[0] +
                xoffset
            )
            ax.plot(
                x_swaxs,
                seg_swaxs.Rg,
                color=color,
                lw=2,
                label=f'{medoid_time:.1f} ns'
            )
            mask_gyrate = (
                (dfrg.time >= medoid_time - window) &
                (dfrg.time <= medoid_time)
            )
            seg_gyrate = dfrg.loc[mask_gyrate].copy()
            if len(seg_gyrate):
                x_gyrate = (
                    seg_gyrate.time -
                    seg_gyrate.time.iloc[0] +
                    xoffset
                )
                ax.plot(
                    x_gyrate,
                    seg_gyrate.Rg,
                    color=color,
                    ls='--',
                    alpha=0.7
                )
            xoffset += window
        ax.set_ylabel(r'$R_G$ [nm]')
        ax.set_xticks([])
        ax.set_xlabel('')
        ax.legend(title='Best-fit intervals')
        return ax

    def plot_rg_chi2(self, ax=None, medoid_times=None, gyrate=True):
        color_chi2 = jocolors.tstum1
        color_rg = jocolors.tstum7
        rgref = SAXS_Measurement.get_rg(self, df=self.get_data(), qmin=0.3, qmax=0.8)['Rg']
        print('rgref', rgref)
        dfrg = self.get_rg()
        print('rg_filename', self.rg_filename)
        df = self.get_all()
        if medoid_times is None:
            medoid_times = self.pick_chi2(df)['time']
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        axR = ax.twinx()
        ax.plot(df.time, df.chi2, color=color_chi2)
        axR.plot(df.time, df.Rg, color=jocolors.tstum7, label='GROMACS-SWAXS')
        axR.plot([df.time[0], df.time.max()], [rgref, rgref], color=color_rg,
                linestyle='--', label='Guinier Experiment')
        if gyrate:
            axR.plot(dfrg.time, dfrg.Rg, color=color_rg,
                    linestyle='-', marker='', label='gmx gyrate',
                    alpha=0.5)
        ax.set_xlabel('time [ns]')
        ax.set_ylabel('$\\chi^2$', color=color_chi2)
        axR.set_ylabel('$R_G$ [nm]', color=color_rg)
        # axR.legend()
        m1color = jocolors.tstum1
        m2color = jocolors.tstum5
        m3color = jocolors.tstum3
        m4color = jocolors.tstum7
        medoid_colors = [m1color,
                         m2color,
                         m3color,
                         m4color]
        for medoid_time, color in zip(medoid_times, medoid_colors):
            ax.axvspan(medoid_time-0.5, medoid_time, facecolor=color, alpha=0.7)  # adjust alpha for transparency
            mask_swaxs = (df.time >= medoid_time - 2) & (df.time <= medoid_time)
            mask_gyrate = (dfrg.time >= medoid_time - 2) & (dfrg.time <= medoid_time)
            mean_swaxs = df.loc[mask_swaxs, 'Rg'].mean()
            mean_chi2 = df.loc[mask_swaxs, 'chi2'].mean()
            std_swaxs  = df.loc[mask_swaxs, 'Rg'].std()
            mean_gyrate = dfrg.loc[mask_gyrate, 'Rg'].mean()
            std_gyrate  = dfrg.loc[mask_gyrate, 'Rg'].std()
            print(
                f"Medoid @ {medoid_time:.2f} ns:\n"
                f"Chi2 @ {mean_chi2:.2f}:\n"
                f"  GROMACS-SWAXS: {mean_swaxs:.3f} ± {std_swaxs:.3f} nm\n"
                f"  gmx gyrate:    {mean_gyrate:.3f} ± {std_gyrate:.3f} nm"
            )
        mask_swaxs = (df.time >= 50) & (df.time <= 63)
        mask_gyrate = (dfrg.time >= 50) & (dfrg.time <= 63)
        mean_swaxs = df.loc[mask_swaxs, 'Rg'].mean()
        std_swaxs  = df.loc[mask_swaxs, 'Rg'].std()
        mean_gyrate = dfrg.loc[mask_gyrate, 'Rg'].mean()
        std_gyrate  = dfrg.loc[mask_gyrate, 'Rg'].std()
        print(
            "\n25–50 ns average:\n"
            f"  GROMACS-SWAXS: {mean_swaxs:.3f} ± {std_swaxs:.3f} nm\n"
            f"  gmx gyrate:    {mean_gyrate:.3f} ± {std_gyrate:.3f} nm"
        )
        fig_best, ax_best = plt.subplots(figsize=(6, 3))
        self._plot_bestfit_rg_segments(
            df=df,
            dfrg=dfrg,
            medoid_times=medoid_times,
            window=3.0,  # change to 3.0 if desired
            ax=ax_best
        )
        outdict = {'axchi2': ax,
                   'axrg': axR,
                   'fig': ax.get_figure(),
                   'fig_best': fig_best,
                   'ax_best': ax_best}
        return outdict


    def plot_waxspot3D(self, cmap='viridis', **kwargs):
        figsize = (self.plot_standardwidth, self.plot_standardheight)
        if 'figsize' in kwargs:
            figsize = kwargs.pop('figsize')
        df = self.get_swaxspot(**kwargs)
        from mpl_toolkits.mplot3d import Axes3D  # not strictly necessary, but keeps IDEs happy
        q_vals = df.columns[1:].astype(float)
        time_vals = df.time
        Z = df.iloc[:, 1:].values
        T, Q = np.meshgrid(time_vals, q_vals, indexing="ij")
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(T, Q, Z, cmap=cmap)
        ax.invert_yaxis()
        ax.set_xlabel("Time [ns]")
        ax.set_ylabel("$q$ [Å$^{-1}$]")
        ax.set_zlabel("$E_{saxs}$ [kJ/mol]")
        return ax

    def get_average_spectra(self, filename=None, path=None,
                     plot=True, get_Rg=False, maxqfit=np.inf,
                     rerun_filename=None,
                     time_interval=(0, np.inf), angular_unit='A',
                      every_n=1, max_out=1000, index_list=None,
                            **kwargs):
        if path is None:
            path = self.mdpath
        if filename is None:
            filename = self.spectra_filename
        outdict = {'time': [],
                   'df': [],
                   'chi2': [],
                   'Rg': [],
                   'err_Rg': [],
                   'Rchi2': [],
                   'I0': [],
                   'err_I0': []}
        file = Path(path, filename)
        if filename is None:
            file = self.get_spectra_filename()
        df_data = self.get_data(angular_unit=angular_unit, scale=self.shift)
        if 'ax' in kwargs:
            ax = kwargs.pop('ax')
        if 'ax_res' in kwargs:
            ax = kwargs.pop('ax_res')
        out = self.get_waxs_spectra(file, every_n=every_n, max_out=max_out)
        df_data['I'] = df_data.I
        df_data['err_I'] = df_data.err_I
        dfs = out['dfs'][1::]
        times = out['times'][1::]
        if not rerun_filename is None:
            rerun_out = self.get_waxs_spectra(Path(path, rerun_filename), every_n=every_n, max_out=max_out)
            rerun_dfs = rerun_out['dfs'][1::]
        if not index_list is None:
            dfs = [out['dfs'][i] for i in index_list]
            times = [out['times'][i] for i in index_list]
        dfs_for_avg = []
        for df, time in zip(dfs, times):
            time_ns = time * 2e-6
            if time_ns<time_interval[0]:
                continue
            if time_ns>time_interval[1]:
                continue
            df['fit'] = df.I
            df = df[df.q>0.15]
            if angular_unit == "A":  # Ångström
                df.loc[:, "q"] /= 10
            dfs_for_avg.append(df)
        df = self.get_average_df(dfs_for_avg)

        scale = None
        offset = None
             # p0=None, scale=None, offset=None, bounds=(-np.inf, np.inf))
        print('maxqfit', maxqfit, 'dfqmax', df.q.max())
        if maxqfit < df.q.max():
            if rerun_filename is None:
                df_maxqfit = df.copy()
            else:
                df_maxqfit = ''
            df_maxqfit = df_maxqfit[df_maxqfit.q < maxqfit]
            odict = self.scale_and_offset_fit(df_data, df_maxqfit)
            scale = odict['scale']
            offset = odict['offset']
            print('maxq', scale, offset)
        odict = self.scale_and_offset_fit(df_data, df, scale=scale, offset=offset)
        df = odict['df']
        chi2 = odict['chi2']
        if get_Rg:
            qmin, qmax = 0.1, 1
            dfr = df.copy()
            if angular_unit == 'A':
                dfr.q = dfr.q * 10
            Rgdict = SAXS_Measurement.get_rg(self, df=dfr, plot=plot, ax=ax,
                                             qmin=0.1, qmax=1)
            Rg = Rgdict['Rg']
            # print('rgdict', Rgdict)
            for key in Rgdict.keys():
                # print(Rgdict[key])
                if key == 'chi2':
                    outdict['Rchi2'].append(Rgdict[key])
                else:
                    outdict[key].append(Rgdict[key])
        outdict['df'].append(df)
        outdict['chi2'].append(chi2)
        return outdict

    def plot_spectra(self, filename=None, path=None, plot=True,
                     get_Rg=False, maxqfit=np.inf, rerun_filename=None,
                     Rg_filter=(0, np.inf), chi2_filter=(0, np.inf),
                     angular_unit='A', every_n=1, max_out=1000, index_list=None,
                     spectra_colors=None, spectra_legend=True,
                     inset=False, inset_xlim=None, inset_ylim=None,
                     spectra_timeis=[0,200], q_cutoff=0.15, **kwargs):
        path = self.mdpath if path is None else path
        filename = self.spectra_filename if filename is None else filename
        file = Path(path, filename) if filename is not None else self.get_spectra_filename()
        outdict = {'time': [], 'df': [], 'chi2': [], 'Rg': [],
                   'err_Rg': [], 'Rchi2': [], 'I0': [], 'err_I0': []}
        df_data = self.get_data(angular_unit=angular_unit, scale=self.shift)
        ax = kwargs.pop('ax') if 'ax' in kwargs else None
        ax_res = kwargs.pop('ax_res') if 'ax_res' in kwargs else None
        if plot:
            if ax is None:
                fig = plt.figure(figsize=(self.plot_standardwidth, self.plot_longheight))
                gs = fig.add_gridspec(2, 1, height_ratios=[3, 1])
                ax = fig.add_subplot(gs[0])
                ax_res = fig.add_subplot(gs[1], sharex=ax)
            else:
                fig = ax.get_figure()
            ax_inset = None
            if inset:
                ax_inset = inset_axes(ax, width="50%", height="50%", loc="upper right")
                ax_inset.set_xscale('log')
                ax_inset.set_yscale('log')
                ax_inset.tick_params(axis='both', which='both', direction='in',
                                     pad=2.0)
                ax_inset = self.plot_data(ax=ax_inset, label=f'{self.name}', marker='o', linestyle='', scale=self.shift, unit=angular_unit, **kwargs)
            if ax_res is None:
                _, ax_res = plt.subplots()
            ax_res.plot([df_data.q[df_data.q<maxqfit].min(), df_data[df_data.q<maxqfit].q.max()], [0, 0])
            ax = self.plot_data(ax=ax, label=f'{self.name}', marker=self.marker, linestyle='', scale=self.shift, unit=angular_unit, **kwargs)
        out = self.get_waxs_spectra(file, every_n=every_n, max_out=max_out)
        dfs, times = out['dfs'][1::], out['times'][1::]
        print('il', index_list)
        if index_list is not None:
            dfs = [dfs[i] for i in index_list]
            times = [times[i] for i in index_list]
        for i, (df, time) in enumerate(zip(dfs, times)):
            df = self._prepare_spectrum(df, angular_unit=angular_unit, q_cutoff=q_cutoff)
            df, chi2 = self._fit_spectrum(df_data, df, maxqfit=maxqfit)
            if not (chi2_filter[0] <= chi2 <= chi2_filter[1]): continue
            time_ns = time * 2e-6
            if get_Rg:
                Rgdict = self._compute_rg(df, angular_unit=angular_unit, plot=plot, ax=ax)
                Rg = Rgdict['Rg']
                if not (Rg_filter[0] <= Rg <= Rg_filter[1]): continue
                for k, v in Rgdict.items():
                    if k == 'chi2': outdict['Rchi2'].append(v)
                    else: outdict[k].append(v)
            outdict['df'].append(df)
            outdict['chi2'].append(chi2)
            outdict['time'].append(time_ns)
            if plot:
                # self._plot_single(ax, ax_res, df, chi2, time_ns, i, spectra_colors, spectra_legend)
                self._plot_single(ax, ax_res, df, chi2, time_ns, i, spectra_colors, spectra_legend,
                                  maxqfit=maxqfit, ax_inset=ax_inset)
        if plot:
            ax.legend(fontsize='x-small')
            ax_res.set_xlabel(ax.get_xlabel())
            ax_res.set_ylabel('Residuals')
            # if inset:
                # ax_inset.set_xticks([])
                # ax_inset.set_yticks([])
            outdict['ax_inset'] = ax_inset
            outdict['fig'], outdict['ax'], outdict['ax_res'] = fig, ax, ax_res
        return outdict

    def _prepare_spectrum(self, df, angular_unit='A', q_cutoff=0.15):
        df = df.copy()
        df['fit'] = df.I
        df = df[df.q > q_cutoff]
        if angular_unit == "A": df.loc[:, "q"] /= 10
        return df

    def _fit_spectrum(self, df_data, df, maxqfit=np.inf):
        scale = offset = chi2 = None
        if maxqfit < df.q.max():
            df_fit = df[df.q < maxqfit]
            odict = self.scale_and_offset_fit(df_data, df_fit)
            scale, offset, chi2 = odict['scale'], odict['offset'], odict['chi2']
        odict = self.scale_and_offset_fit(df_data, df, scale=scale, offset=offset)
        return odict['df'], chi2 if chi2 is not None else odict['chi2']

    def _compute_rg(self, df, angular_unit='A', plot=False, ax=None):
        qmin, qmax = max(0.01, self.qmin), max(0.04, self.qmin + 0.6)
        dfr = df.copy()
        if angular_unit == 'A': dfr.q *= 10; qmin, qmax = 10*qmin, 10*qmax
        return SAXS_Measurement.get_rg(self, df=dfr, plot=plot, ax=ax, qmin=qmin, qmax=qmax)

    def _plot_single(self, ax, ax_res, df, chi2, time_ns, i, spectra_colors, spectra_legend, maxqfit=np.inf, ax_inset=None):
        color = spectra_colors[i] if spectra_colors is not None else None
        label = spectra_legend[i] if isinstance(spectra_legend, list) else (f'$t = {time_ns:.1f}$ ns; $\\chi^2={chi2:.1f}$' if spectra_legend else None)
        ax.errorbar(df.q, df.I, fmt='-', color=color, label=label)
        ax_res.errorbar(df[df.q<maxqfit].q, df[df.q<maxqfit].res, fmt='-', color=color)
        if ax_inset is not None: ax_inset.errorbar(df.q, df.I, fmt='-', color=color, label=label)

    def get_Ree(self, selection="nucleic and name P",
                time_interval=(0, np.inf), step=1):
        u = self.get_u()
        sel = u.select_atoms(selection)
        out = {"time": [], "Ree": []}
        for ts in u.trajectory[::step]:
            t_ns = ts.time / 1000
            if not (time_interval[0] <= t_ns <= time_interval[1]):
                continue
            r0 = sel.positions[0]
            rN = sel.positions[-1]
            Ree = np.linalg.norm(rN - r0)
            out["time"].append(t_ns)
            out["Ree"].append(Ree)
        return pd.DataFrame(out)

    def get_Ree_stats(self, **kwargs):
        df = self.get_Ree(**kwargs)
        Ree2 = df.Ree.values**2
        return {
            "Ree_mean": df.Ree.mean(),
            "Ree_std": df.Ree.std(),
            "Ree2_mean": Ree2.mean(),
            "Ree2_std": Ree2.std()
        }

    def get_all(self, mdpath=None,
                plot_spectra=False,
                spectra_filename=None,
                proj_filename=None,
                rmsf_filename=None,
                waxspot_name=None):
        if mdpath is None:
            mdpath = self.mdpath
        if spectra_filename is None:
            spectra_filename = self.spectra_filename
        if proj_filename is None:
            proj_filename = self.proj_filename
        if rmsf_filename is None:
            rmsf_filename = self.rmsf_filename
        print(spectra_filename)
        try:
            spectra_dict = self.plot_spectra(filename=spectra_filename,
                                         path=mdpath,
                                         every_n=1,
                                         max_out=100000,
                                         plot=plot_spectra,
                                         angular_unit='nm',
                                         get_Rg=True)
            for key in spectra_dict.keys():
                print(key, len(spectra_dict[key]))
            spectra_df = pd.DataFrame({'time': spectra_dict['time'],
                                   'chi2': spectra_dict['chi2'],
                                   'I0': spectra_dict['I0'],
                                   'err_I0': spectra_dict['err_I0'],
                                   'Rg': spectra_dict['Rg'],
                                   'err_Rg': spectra_dict['err_Rg']})
            # print(spectra_dict)
        except Exception as e:
            print('failed to load spectra dict:')
            print(e)
            spectra_filename = None
        try:
            proj_df = self.get_proj(mdpath=mdpath, filename=proj_filename)
            # print(proj_df)
        except Exception as e:
            proj_filename = None
            print('no proj')
        if not rmsf_filename is None:
            try:
                rmsf_df = self.get_proj(mdpath=mdpath, filename=rmsf_filename, pc='rmsf')
            except:
                print(f'no file called {rmsf_filename}')
        # print('proj_df', proj_df)
        # print('spectra_df', spectra_df)
        if proj_filename is None:
            return spectra_df
        elif spectra_filename is None:
            return proj_df
        else:
            merged_df = pd.merge_asof(proj_df, spectra_df[['time', 'chi2', 'Rg']], on='time', direction='forward')
            return merged_df

    def plot_rg(self, ax=None, rgkwargs={}, **kwargs):
        out = {}
        if ax is None:
            fig, ax = plt.subplots(figsize=(3,1))
        df = self.get_rg(**rgkwargs)
        out['df'] = df
        annotation = '{}, $R_g = {:.2f} \\pm {:.2f}$ nm'.format(
                        self.name, df.Rg.mean(), df.Rg.std())
        ax.errorbar(df.time, df.Rg,
                    marker='',
                    linestyle='-',
                    **kwargs)
        ax.annotate(annotation, xy=(0.01, 0.8), xycoords='axes fraction',
                    fontsize='small', color='#e5e7fa', bbox=dict(facecolor='#6e79ad'))

        ax.set_ylabel('$R_G$ [nm]')
        ax.set_xlabel('time [ns]')
        out['ax'] = ax
        return out

    def plot_pick_chi2(self, n=5, min_dt=10):
        df = self.get_all()
        pick_chi2df = self.pick_chi2(df, n=5, min_dt=10)
        pick_chi2 = [-int((-10*t)//6) for t in pick_chi2df.time]
        print(pick_chi2df)
        print(pick_chi2)
        outdict = self.plot_spectra(maxqfit=19.0,
                          every_n=1, index_list=pick_chi2)
        ax = outdict['ax']
        ax_res = outdict['ax_res']
        ax.legend()
        ax_res.legend('')
        plt.tight_layout()
        return outdict

    def plot_best_worst(self, best_kwargs={}, worst_kwargs={},
                        average_kwargs={},
                        worst=False, average=True, best=True,
                        spectra_legend=False, **kwargs):
        df = self.get_all()
        best_chi2df = self.pick_chi2(df, n=1)
        if worst:
            worst_chi2df = self.pick_chi2(df, n=1, min_time=20, ascending=False)
        pick_chi2 = [-int((-10*t)//6) for t in best_chi2df.time]
        if not best:
            pick_chi2 = []
        if worst:
            pick_chi2 +=\
                    [-int((-10*t)//6) for t in worst_chi2df.time]
        print(pick_chi2)
        ax, ax_res = None, None
        if 'ax' in kwargs:
            ax = kwargs.pop('ax')
        if 'ax_res' in kwargs:
            ax_res = kwargs.pop('ax_res')
        outdict = self.plot_spectra(maxqfit=1.7,
                                    every_n=1, index_list=pick_chi2,
                                    ax=ax, ax_res=ax_res,
                                    spectra_legend=spectra_legend,
                                    spectra_colors=[jocolors.tstum22,
                                                    jocolors.tstum8])
        if not average:
            return outdict
        ax = outdict['ax']
        ax_res = outdict['ax_res']
        average_dict = self.get_average_spectra(maxqfit=19.0, every_n=1,
                                           time_interval=(100, 150))
        average_df = average_dict['df'][0]
        average_chi2 = average_dict['chi2'][0]
        print(average_kwargs)
        if not 'label' in average_kwargs:
            avg_label = f'average; $\\chi^2={average_chi2:.1f}$'
            average_kwargs['label'] = avg_label
        if not 'color' in average_kwargs:
            average_kwargs['color'] = jocolors.tstum23
        if not 'marker' in average_kwargs:
            average_kwargs['marker'] = ''
        if not 'linestyle' in average_kwargs:
            average_kwargs['linestyle'] = '-'
        ax.errorbar(average_df.q, average_df.I, **average_kwargs)
        ax.legend()
        ax_res.legend('')
        plt.tight_layout()
        return outdict

    @staticmethod
    def read_xvg(filename, include1=False, every_n=1, max_out=np.inf,
                 names=['q', 'I', 'errI', 'I2', 'I3', 'I4'],
                 onlyfirst=False):
        out = {'dfs': [],
               'times': []}
        with open(filename) as f:
            start_row = 0
            time = 1
            plot = 0 if include1 else 1
            for i, line in enumerate(f):
                if line[0] == '@':
                    start_row = i+1
                elif line[0] == '#':
                    start_row = i+1
                elif line == '\n':
                    start_row = i+1
                elif line == '&\n':
                    end_row = i-2
                    if start_row>0 and (plot%every_n)==0 and len(out['dfs'])<max_out:
                        df = pd.read_csv(filename,
                                         skiprows=start_row,
                                         nrows=end_row-start_row, sep=r'\s+',
                                         names=names)
                        if onlyfirst:
                            return df
                        out['dfs'].append(df)
                        out['times'].append(time)
                    time += 1
                    plot += 1
            df = pd.read_csv(filename, skiprows=start_row, sep=r'\s+',
                            names=names)
            out['dfs'].append(df)
            out['times'].append(time)
        return out

    @staticmethod
    def eta_filter(x, dt, tau):
        alpha = dt / tau
        y = np.zeros_like(x)
        y[0] = x[0]
        for i in range(1, len(x)):
            y[i] = alpha * x[i] + (1 - alpha) * y[i-1]
        return y

    @staticmethod
    def get_crysol(path):
        df = pd.read_csv(path, skiprows=2, sep=r'\s+',
                         names=['q', 'I', 'errI', 'fit'])
        return df

    @staticmethod
    def get_average_df(list_of_dfs):
        # Concatenate all DataFrames
        df_all = pd.concat(list_of_dfs, ignore_index=True)

        # Weighted average per q
        def weighted_avg(group):
            I = group["I"].values
            err = group["err_I"].values
            w = 1 / (err**2)  # inverse variance weights
            I_avg = np.sum(I * w) / np.sum(w)
            err_avg = np.sqrt(1 / np.sum(w))
            return pd.Series({"I": I_avg, "err_I": err_avg})

        df_avg = df_all.groupby("q").apply(weighted_avg).reset_index()
        return df_avg

    def save_time_time_rmsd(self,
                            selection="nucleic",
                            timestep=1,
                            out_file="rmsd_matrix.txt",
                            tmp_coords="coords.dat"):
        """
        Compute a time-time RMSD matrix for the entire trajectory (sampled),
        using no-fit coordinates (trajectory assumed pre-fitted externally),
        and store results line by line to avoid RAM blow-up.

        Parameters
        ----------
        selection : str
            MDAnalysis atom selection.
        timestep : int
            Frame stride (e.g., 10 means every 10th frame).
        out_file : str
            Output text file where matrix rows are appended.
        tmp_coords : str
            Temporary file for numpy.memmap coordinate storage.
        """
        import numpy as np
        import MDAnalysis as mda
        from pathlib import Path
        u = self.get_u()  # assumes already nojump+fit
        atoms = u.select_atoms(selection)
        # --------------------------------------------------------------------------
        # PASS 1: COLLECT COORDS INTO MEMMAP
        # --------------------------------------------------------------------------
        print("[1/2] Sampling trajectory...")
        # First count frames
        n_frames = 0
        for i, ts in enumerate(u.trajectory):
            if i % timestep == 0:
                n_frames += 1

        n_atoms = atoms.n_atoms
        ndim = n_atoms * 3
        # Prepare memory map for coordinates
        tmp_coords = Path(tmp_coords)
        if tmp_coords.exists():
            tmp_coords.unlink()  # ensure clean
        coords_mm = np.memmap(tmp_coords, dtype='float32', mode='w+',
                              shape=(n_frames, ndim))
        # Fill memmap
        idx = 0
        for i, ts in enumerate(u.trajectory):
            print('wtf long')
            if i % timestep != 0:
                continue
            pos = atoms.positions.astype('float32').reshape(-1)
            coords_mm[idx, :] = pos
            idx += 1
        coords_mm.flush()
        print(f"  Stored {n_frames} sampled frames in memmap at '{tmp_coords}'")
        # --------------------------------------------------------------------------
        # PASS 2: COMPUTE RMSD MATRIX ROW-BY-ROW
        # --------------------------------------------------------------------------
        print("[2/2] Computing pairwise RMSD matrix...")
        out_file = Path(out_file)
        if out_file.exists():
            out_file.unlink()  # overwrite if exists

        # Re-open as readonly for safety
        coords_mm = np.memmap(tmp_coords, dtype='float32', mode='r',
                              shape=(n_frames, ndim))
        with out_file.open("w") as fout:
            for i in range(n_frames):
                # Load ith frame once
                ci = coords_mm[i]

                # Compute RMSD(i,j) for j>=0
                diffs = coords_mm - ci  # broadcast (n_frames, ndim)
                rmsd_vec = np.sqrt(np.mean(diffs**2, axis=1))

                # Write one line: space-separated RMSD values
                fout.write(" ".join(f"{v:.4f}" for v in rmsd_vec) + "\n")

                if (i+1) % 10 == 0:
                    print(f"  Row {i+1}/{n_frames} done...")

        print(f"Done! RMSD matrix written to '{out_file}'")
        print("Pro tip: use `np.loadtxt` or pandas to visualize later.")


from dataclasses import dataclass, field
import numpy as np
import pandas as pd


@dataclass
class Analysis:

    pars: list[str] = field(default_factory=list)
    df_pars: list[str] = field(default_factory=list)

    name: str = "Base Class"
    selection_solute: str = (
        "not resname SOL NA CL and not name H*"
    )
    selection_bb: str = "name P"
    selection_solvent: str = "resname SOL and not name H*"
    selection_ions: str = "resname NA CL"

    def run(self, u, *, step=1, **kwargs):
        solute = u.select_atoms(self.selection_solute)
        bb_atoms = u.select_atoms(self.selection_bb)
        solvent = u.select_atoms(self.selection_solvent)
        ions = u.select_atoms(self.selection_ions)
        results = []
        aux_results = {
            par: []
            for par in self.df_pars
        }
        # for ts in u.trajectory[::step]:
        for ts in tqdm(u.trajectory[::step], desc="run"):
            frame_results, frame_aux = self.calc_function(
                solute=solute,
                bb_atoms=bb_atoms,
                solvent=solvent,
                ions=ions,
                ts=ts,
                **kwargs,
            )
            results.append(frame_results)
            for par in self.df_pars:
                aux_results[par].append(frame_aux[par])
        mean_aux_results = {
            par: np.mean(aux_results[par], axis=0)
            for par in self.df_pars
        }
        return pd.DataFrame(results), mean_aux_results

    def calc_function(self, *args, **kwargs):
        raise NotImplementedError

    def save_json(
        self,
        filename,
        results,
        mean_aux_results=None,
    ):
        output_dir = Path("hydration_output")
        output_dir.mkdir(exist_ok=True)
        results.to_json(
            output_dir / f"{filename}_results.json",
            orient="table",
        )
        if mean_aux_results is not None:
            for dfpar in self.df_pars:
                aux_df = pd.DataFrame(mean_aux_results[dfpar])
                aux_df.to_json(
                    output_dir / f"{filename}_{dfpar}.json",
                    orient="table",
                )

    def load_json(self, filename):
        output_dir = Path("hydration_output")
        results = pd.read_json(
            output_dir / f"{filename}_results.json",
            orient="table",
        )
        mean_aux_results = pd.read_json(
            output_dir / f"{filename}_aux.json",
            orient="table",
        )
        return results, mean_aux_results

@dataclass
class HydrationAnalysis(Analysis):

    pars: list[str] = field(default_factory=lambda: [
        "time",
        "V_eff",
        "rho_bulk",
        "rho_solute",
    ])
    df_pars: list[str] = field(default_factory=lambda: [
        "rho_hist"
        ])
    name: str = "hydration"
    selection_solvent: str = "resname SOL and name OW"

    def calc_function(
        self,
        solute,
        bb_atoms,
        solvent,
        ions,
        ts,
        r0=30.0,
        r_min_bulk=30.0,
        r_max_bulk=35.0,
        dr=0.5,
    ):
        V_sphere = (4 / 3) * np.pi * r0**3
        ref_pos = solute.center_of_mass()
        diff = solvent.positions - ref_pos
        r = np.linalg.norm(diff, axis=1)
        inside = r < r0
        N_in = np.sum(inside)
        bins = np.arange(0, r_max_bulk + dr, dr)
        hist, edges = np.histogram(r, bins=bins)
        r_centers = 0.5 * (edges[:-1] + edges[1:])
        shell_volumes = (
            4 / 3
            * np.pi
            * (edges[1:]**3 - edges[:-1]**3)
        )
        rho_r = hist / shell_volumes
        solute_mask = r_centers <= r0
        rho_solute = (
            np.sum(hist[solute_mask])
            / np.sum(shell_volumes[solute_mask])
        )
        rho_hist = hist / shell_volumes
        rho_hist_df = pd.DataFrame({
            'r': rho_r,
            'rho': rho_hist})
        bulk_mask = (
            (r_centers >= r_min_bulk)
            & (r_centers <= r_max_bulk)
        )
        rho_bulk = (
            np.sum(hist[bulk_mask])
            / np.sum(shell_volumes[bulk_mask])
        )
        V_eff_rho = (
            (rho_bulk - rho_solute)
            * V_sphere
            / rho_bulk
        )
        results = {
            "time": ts.time,
            "V_eff": V_eff_rho,
            "rho_bulk": rho_bulk,
            "rho_solute": rho_solute,
        }
        frame_aux = {'rho_hist': rho_hist_df}
        return results, frame_aux

@dataclass
class FbsAnalysis(Analysis):

    pars: list[str] = field(default_factory=lambda: [
        "time",
        "fbs",
        "distance_matrices",
    ])
    df_pars: list[str] = field(default_factory=lambda: [
        "gbs"
        ])
    name: str = "BS"
    select_group: str = f"not name P O1P OP1 O2P OP2 O5' C5' C4' O4' C3' O3' C2' C1'"

    def calc_function(
        self,
        solute,
        bb_atoms,
        solvent,
        ions,
        ts,
        r_max=10.0,
        dr=0.1,
        cutoff=5.0,
    ):
        n_res = len(solute.residues)
        base_coms = np.array([
            residue.atoms.select_atoms(self.select_group).center_of_mass()
            for residue in solute.residues
        ])
        dist_mat = np.full((n_res, n_res), np.nan)
        for i in range(n_res):
            for j in range(i + 1, n_res):
                if i == 12 and j == 13:
                    continue
                dist = distance_array(
                    base_coms[i][None, :],
                    base_coms[j][None, :],
                    box=ts.dimensions,
                )[0, 0]
                dist_mat[i, j] = dist
                dist_mat[j, i] = dist
        neighbor_dists = np.array([
            dist_mat[i, i + 1]
            for i in range(n_res - 1)
            if i != 12
        ])
        r_edges = np.arange(0, r_max + dr, dr)
        r = 0.5 * (r_edges[:-1] + r_edges[1:])
        gbs = np.histogram(neighbor_dists, bins=r_edges, density=True)[0]
        gbs_df = pd.DataFrame({'r': r, 'gbs': gbs})
        fbs = np.sum(neighbor_dists < cutoff) / len(neighbor_dists)
        results = {
            "time": ts.time,
            "fbs": fbs,
            "distance_matrices": dist_mat.tolist(),
        }
        frame_aux = {'gbs': gbs_df}
        return results, frame_aux

@dataclass
class RadiusAnalysis(Analysis):

    pars: list[str] = field(default_factory=lambda: [
        "time",
        "Rg",
        "Ree",
    ])
    df_pars: list[str] = field(default_factory=list)
    name: str = "Radius"
    selection_solute: str = (
        "not resname SOL NA CL"
    )

    def calc_function(self, solute, bb_atoms, solvent, ions, ts):
        Rg = solute.radius_of_gyration()
        r0 = bb_atoms.positions[0]
        rN = bb_atoms.positions[-1]
        Ree = np.linalg.norm(rN - r0)
        results = {
            "time": ts.time,
            "Rg": Rg,
            "Ree": Ree,
        }
        frame_aux = {}
        return results, frame_aux
