import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib
from .beaucage import beaucage
from dataclasses import dataclass, field

from .SAXS_Measurement import SAXS_Measurement
from .. import os_utility as osu


def binArray(data, axis, binstep, binsize, func=np.nanmean):
    data = np.array(data)
    dims = np.array(data.shape)
    argdims = np.arange(data.ndim)
    argdims[0], argdims[axis] = argdims[axis], argdims[0]
    data = data.transpose(argdims)
    data = [func(np.take(data, np.arange(int(i*binstep),
            int(i*binstep+binsize)), 0), 0)
            for i in np.arange(dims[axis]//binstep)]
    data = np.array(data).transpose(argdims)
    return data


@dataclass
class Ms:

    ms: list[SAXS_Measurement] = field(default_factory=list)
    model: str = ""
    name: str = ""


    def __index__(self, index):
        return self.ms[index]

    def __iter__(self):
        return iter(self.ms)

    def __len__(self):
        return len(self.ms)

    def build_saxs_matrix(self):
        I_list = []
        sigma_list = []
        q_values = None
        for m in self.ms:
            df = m.get_data(angular_unit='nm')
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


    def get_workdir(self):
        dirname = f'./workdirs/ms_{self.name}'
        return dirname

    def make_workdir(self, dirname=None):
        if dirname is None:
            dirname = self.get_workdir()
        osu.create_path('workdirs')
        osu.create_path(dirname)
        # shutil.copyfile(self.get_filename(), join(dirname, 'data.dat'))

    def make_plot(self, **kwargs):
        if 'ax' in kwargs:
            ax = kwargs.pop('ax')
        else:
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4),
                                   sharex=True, sharey=True, squeeze=True)
        if 'ylim' in kwargs:
            ylim = kwargs.pop('ylim')
            ax.set_ylim(*ylim)
        if 'xlim' in kwargs:
            xlim = kwargs.pop('xlim')
            ax.set_xlim(*xlim)
        ax.tick_params(axis="x", bottom=True, top=True,
                       labelbottom=True, labeltop=False)
        ax.tick_params(axis="y", left=True, right=True,
                       labelleft=True, labelright=False)
        return ax, kwargs

    def fits(self, fits=True, shiftby=1, dataqmin=0, dataqmax=2300,
             error=True, **kwargs):
        ax, kwargs = self.make_plot(**kwargs)
        shift = shiftby**len(self.ms)
        if 'shifts' in kwargs:
            shifts = kwargs.pop('shifts')
        else:
            shifts = [shift / (shiftby ** i) for i in range(len(self.ms))]
        for i, m in enumerate(self.ms):
            shift = shifts[i]
            print(m.label, shift)
            dfall = m.get_data(cout=False)
            df = dfall[m.dataqmin: m.dataqmax]
            if error:
                ax.errorbar(df.q, df.I*shift, df.err_I*shift,
                            marker=m.marker, color=m.color,
                            linestyle='', label=m.label,
                            elinewidth=0.2, **kwargs)
            else:
                ax.plot(df.q, df.I*shift,
                        marker=m.marker, color=m.color,
                        linestyle='', label=m.label, **kwargs)
            df = dfall[m.iqmin: m.iqmax]
            if fits:
                m.fit_dict, m.fit_df = m.model.fit(m, bounds=m.bounds,
                                                   iqmax=m.iqmax, p0=m.p0,
                                                   iqmin=m.iqmin,
                                                   fixed_parameters=m.fixed_pars)
                ax.errorbar(m.fit_df.q, m.fit_df.fit*shift, marker='',
                            color='black')
                m.partext = m.model.get_text(m.fit_dict)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend(fontsize='xx-small')
        ax.set_xlabel('$q \\,\\mathrm{[nm^{-1}]}$')
        ax.set_ylabel('$I \\,\\mathrm{[cm^{-1}]}$')

    def data(self, shiftby=1, **kwargs):
        return self.fits(fits=False, shiftby=shiftby, **kwargs)

    def res(self, **kwargs):
        if 'ax' in kwargs:
            ax = kwargs.pop('ax')
        else:
            fig, ax = plt.subplots(nrows=1, ncols=1,
                                   figsize=(5, 4),
                                   sharex=True, sharey=True,
                                   squeeze=True)
        if 'ylim' in kwargs:
            ylim=kwargs.pop('ylim')
        else:
            ylim = [None, None]
        for m in self.ms:
            label = '{}; $\\chi^2 = {:.1f}$'.format(
                    m.label, m.fit_dict['chi2']
                    )
            df = m.fit_df
            ydata = df.res/df.err_I
            ax.plot(df.q, ydata, marker=m.marker, color=m.color,
                    linestyle='', label=label, **kwargs)
        ax.legend(fontsize='xx-small')
        ax.set_ylabel('Normalized Residuals')
        ax.set_xlabel('$q\\mathrm{\\,[nm^{-1}]}$')

    def kratky(self, **kwargs):
        ax, kwargs = self.make_plot(**kwargs)
        for m in self.ms:
            df = m.get_data(cout=False)[m.iqmin:m.iqmax]
            xdata = df.q
            ydata = 1000 * df.I * df.q**2
            if 'binnum' in kwargs:
                binnum = kwargs.pop('binnum')
                xdata = binArray(df.q, 0, binnum, binnum)
                ydata = 1000 * binArray(df.I, 0, binnum, binnum) * xdata**2
            label = m.label
            ax.errorbar(xdata, ydata, marker=m.marker, color=m.color,
                        linestyle='', label=label, elinewidth=0.2, **kwargs)
            # ax.legend()
            ax.legend(fontsize='xx-small')
            #             ax.annotate(m.partext, xy=(0.0, 0.0), xycoords='axes fraction')
            # ax.grid()
            ax.set_ylim(0,5)
            ax.set_xlim(0, 2.5)
            ax.set_xlabel('$q \\,\\mathrm{[nm^{-1}]}$')
            #         axes[0][1].set_xlabel('$q$ [1/nm]')
            ax.set_ylabel('$I\\cdot q^2 \\mathrm{\\,[0.001nm^{-2}cm^{-1}]}$')

    def q3I(self, **kwargs):
        ax, kwargs = self.make_plot(**kwargs)
        for m in self.ms:
            df = m.get_data(cout=False)[m.iqmin:m.iqmax]
            label = m.label
            ax.errorbar(df.q, df.q**3 * df.I * 10_000,  marker=m.marker, color=m.color,
                        linestyle='', label=label, elinewidth=0.2, **kwargs)
            #             ax.annotate(m.partext, xy=(0.0, 0.0), xycoords='axes fraction')
            # ax.grid()
            ax.set_ylim(0, 0.0025)
            ax.set_xlim(0, 0.6)
            ax.set_xlabel('$q \\,\\mathrm{[nm^{-1}]}$')
            #         axes[0][1].set_xlabel('$q$ [1/nm]')
            ax.set_ylabel('$I\\cdot q^3 \\mathrm{\\,[nm^{-3}cm^{-1}]}$')
            ax.set_ylabel('$I\\cdot q^3 \\mathrm{\\,[A.U.]}$')
        ax.legend(fontsize='xx-small')

    def debye(self, **kwargs):
        pass

    def get_results(self):
        par_dict = {}
        for par in self.model.parameters:
            par_dict[par] = []
            par_dict[f'std_{par}'] = []
            for m in self.ms:
                if par in m.fit_dict:
                    par_dict[par].append(m.fit_dict[par])
                    par_dict[f'std_{par}'].append(m.fit_dict['std_'+par])
                else:
                    par_dict[par].append(None)
                    par_dict[f'std_{par}'].append(None)
                    print(f'{par} not in fit_dict')
        self.df = pd.DataFrame(par_dict)
        return self.df

    def save_fit_results(self, path, **kwargs):
        self.df = self.get_results()
        self.df.to_csv(path)
        return self.df

    def get_rgs_filename(self, name='rgs.dat'):
        out = f'{self.get_workdir()}/{name}'
        return out

    def _load_rgs(self):
        filename = self.get_rgs_filename()
        outdf = pd.read_csv(filename, sep=r'\s+')
        return outdf

    def get_rgs(self, load=True, **kwargs):
        rgs_filename = self.get_rgs_filename()
        if load and os.path.exists(rgs_filename):
            print('rgs loaded from existing file', rgs_filename)
            return self._load_rgs()
        outdict = {'Rg': [],
                   'err_Rg': [],
                   'I0': [],
                   'err_I0': [],
                   'chi2': []}
        for m in self.ms:
            try:
                rdict = m.get_rg(**kwargs)
                for key in rdict.keys():
                    outdict[key].append(rdict[key])
            except Exception as e:
                for key in outdict.keys():
                    outdict[key].append(0)
        outdf = pd.DataFrame(outdict)
        osu.create_path(self.get_workdir())
        outdf.to_csv(rgs_filename, sep='\t', float_format='{:.7e}'.format, index=False)
        return outdf

    def plot_par(self, par, **kwargs):
        ax, kwargs = self.make_plot(**kwargs)
        df = self.get_results()
        xdata = [m.label for m in self.ms]
        ydata = df[par]
        err_ydata = df[f'std_{par}']
        ax.errorbar(xdata, ydata, err_ydata, **kwargs)
        return ax

    def markdown_table(self, fixed_pars=[], e_pars=[], **kwargs):
        out = '| parameter |'
        for m in self.ms:
            out += f'{m.label} |'
        out += '\n | --- |'
        for m in self.ms:
            out += f'--- |'
        out += '\n'
        for par in self.model.parameters:
            out += f'| {par} |'
            for m in self.ms:
                if not f'std_{par}' in m.fit_dict:
                    pass
                elif m.fit_dict[f'std_{par}'] == 'fixed':
                    out += '{:.2e} fix |'.format(m.fit_dict[par])
                elif par in e_pars:
                    print(m.fit_dict[f'std_{par}'])
                    out += '{:.2e} |'.format(m.fit_dict[par], m.fit_dict['std_'+par])
                else:
                    out += '{:.3g} ± {:.3g} |'.format(m.fit_dict[par], m.fit_dict['std_'+par])
            out += '\n'
        if 'beaucage_scale' in self.model.parameters:
            out += f'| beaucage_C |'
            for m in self.ms:
                out += '{:.2e} |'.format(beaucage.get_C(m.fit_dict['beaucage_scale'],
                                                        m.fit_dict['beaucage_rg'],
                                                        m.fit_dict['beaucage_exp']))
            out += '\n'
        return out

    def markdown_sasview(self):
        pass

    def markdown_latex_table(self, **kwargs):
        out = '| parameter |'
        for m in self.ms:
            out += f'{m.label} |'
        out += '\n | --- |'
        for m in self.ms:
            out += f' --- |'
        out += '\n'
        for par in self.model.parameters:
            out += f'| ${self.model.pardict["par"]["latex"]}$ |'
            for m in self.ms:
                out += '{:.2f} \\pm {:.2f} |'.format(m.fit_dict[par], m.fit_dict['std_'+par])
            out += '\n'
        return out

    def to_excel(self, fits=True, onlyraw=False,
                 outpath='~/test_saxs_to_excel.xlsx',
                 shiftby=10, **kwargs):
        datadic = {}
        fitdic = {}
        colordic = {'labels': [m.label for m in self.ms],
                    'color': [matplotlib.colors.to_hex(
                        m.color) for m in self.ms]}
        shift = shiftby ** len(self.ms)
        if 'shifts' in kwargs:
            shifts = kwargs.pop('shifts')
        else:
            shifts = [shift / (shiftby ** i) for i in range(len(self.ms))]
        for i, m in enumerate(self.ms):
            shift = shifts[i]
            dfall = m.get_data(cout=False)
            df = dfall[m.dataqmin: m.dataqmax]
            datadic[f'{m.label}_q'] = df.q
            datadic[f'{m.label}'] = df.I
            datadic[f'{m.label}_err_I'] = df.err_I
            if fits:
                m.fit_dict, m.fit_df = m.model.fit(
                        m, bounds=m.bounds,
                        iqmax=m.iqmax, p0=m.p0,
                        iqmin=m.iqmin,
                        fixed_parameters=m.fixed_pars)
                fitdic[f'{m.label}_q'] = df.q
                fitdic[f'{m.label}'] = df.I * shift
                fitdic[f'{m.label}_err_I'] = df.err_I * shift
                fitdic[f'{m.label}_q_fit'] = m.fit_df.q
                fitdic[f'{m.label}_I_fit'] = m.fit_df.fit * shift
        df_data = pd.DataFrame(datadic)
        df_color = pd.DataFrame(colordic)
        if fits:
            df_fit = pd.DataFrame(fitdic)
            df_par = self.get_results()
        with pd.ExcelWriter(outpath) as writer:
            df_data.to_excel(writer, sheet_name='data')
            df_color.to_excel(writer, sheet_name='color')
            if fits:
                df_fit.to_excel(writer, sheet_name='fit')
                df_par.to_excel(writer, sheet_name='par')
