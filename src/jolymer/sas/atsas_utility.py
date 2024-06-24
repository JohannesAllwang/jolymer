import subprocess
import getpass
import io
import os
from os.path import join
import pandas as pd
import re

from dataclasses import dataclass

atsas_bin = f"C:\\Users\\{getpass.getuser()}\\atsas\\bin"

class ATSAS:

    path: str = '~'
    filename: str='data.dat'
    angular_unit: str='$\\mathrm{nm^{-1}}$'
    color: str=None
    marker: str=''
    linestyle: str='-'
    cwd: str='.'

    def run_program(self, program, *args, cwd=None, **kwargs):
        runcwd = self.cwd
        if cwd is None:
            runcwd = self.cwd
        subprocess.run([f'/home/johannes/repos/ATSAS-3.2.0-1/bin/{program}',
                        *args], cwd=runcwd)

    def run_crysol(self, file=None, datafile=None, prefix='test',
                   smax=2, units=1, cwd='crysol', ns=500,
                   other_arguments=['--explicit-hydrogens']):
        if file is None:
            file = join(self.cwd, self.filename)
        if datafile is None:
            datafile = self.filepath

        self.run_program('crysol',
                       file, datafile, f'--ns={ns}', f'--smax={smax}',
                       f'--units={units}',
                       *other_arguments, '--lm=60',
                       f'--prefix={prefix}')

    def run_gnom(self, path=None, filename=None, prefix='test',
                 outfile='gnom.out', pofrfile='gnom.pofr', fitfile='gnom.fit',
                 radius=500, other_arguments=[]):
        if filename is None:
            filename = self.filename
        if path is None:
            path = self.cwd

        filepath = filename
        # filepath = join(path, filename)
        # outfile = join(path, outfile)
        # pofrfile = join(path, pofrfile)
        # fitfile = join(path, fitfile)
        # print('cwd2', self.cwd)

        self.run_program('datgnom', filepath, '-o', outfile,
                         '-r', f'{radius}', *other_arguments)
        self.run_program('out2pofr', outfile, '-o', pofrfile, *other_arguments)
        self.run_program('out2fit', outfile, '-o', fitfile, *other_arguments)

    def get_fit(self, path, engine='pandas', **kwargs):
        """
        Reads the data file at self.path / self.filename
        Alternatively, a path can be provided as a keyword argument.
        """
        df = pd.read_csv(path, sep='\s+', header=None, skipfooter=20,
                         names=['q', 'I', 'errI', 'fit'])
        return df

    def get_pofr(self, path, engine='pandas', **kwargs):
        """
        Reads the data file at self.path / self.filename
        Alternatively, a path can be provided as a keyword argument.
        """
        df = pd.read_csv(path, sep='\s+', header=None,
                         names=['r', 'pr', 'errpr'])
        return df

    def get_rg(self, path):
        with open(path, 'r') as file:
            content = file.read()
            # Using regular expression to find the Real space Rg value
            match = re.search(r'Real space Rg:\s+([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*\+\-\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', content)
            if match:
                rg_value = float(match.group(1))
                uncertainty = float(match.group(2))
                return rg_value, uncertainty
            else:
                return None, None



def run_crysol(file='', datafile='', prefix='test',
                   smax=2, units=1, cwd='crysol', ns=500,
               other_arguments=['--explicit-hydrogens']):
    subprocess.run(['/home/johannes/repos/ATSAS-3.2.0-1/bin/crysol',
                    file, datafile, f'--ns={ns}', f'--smax={smax}',
                    f'--units={units}',
                    *other_arguments, '--lm=60',
                    f'--prefix={prefix}'],
                   cwd=cwd)
    return file

def get_outpath(m, name):
    outfile = r"{}".replace('\\\\', '\\').replace('/', '\\').format(join(m.path, name)).replace('/', '\\')
    return outfile

def run_program(name, *args):
    cmd = ['powershell', f'./{name}.exe', *args]
    return subprocess.run(cmd, cwd=atsas_bin)

def primus(m, toplot=[ 'processed_averaged' ]):
    allnames = ""
    if 'processed_averaged' in toplot:
        fullpath = r"{}".replace('\\\\', '\\').replace('/', '\\').format(join(m.get_filename())).replace('/', '\\')
        allnames += f" '{fullpath}'"
    if 'averaged_sample' in toplot:
        fil = m.get_averaged_fullpath(buf=False)
        fullpath = r"{}".replace('\\\\', '\\').replace('/', '\\').format(fil).replace('/', '\\')
        allnames += f" '{fullpath}'"
    if 'averaged_buffer' in toplot:
        fil = m.get_averaged_fullpath(buf=True)
        fullpath = r"{}".replace('\\\\', '\\').replace('/', '\\').format(fil).replace('/', '\\')
        allnames += f" '{fullpath}'"
    if 'absolutes_sample' in toplot:
        for fil in m.get_absolute_fullpaths(buf=False):
            fullpath = r"{}".replace('\\\\', '\\').replace('/', '\\').format(fil).replace('/', '\\')
            allnames += f" '{fullpath}'"
    if 'absolutes_buffer' in toplot:
        for fil in m.get_absolute_fullpaths(buf=True):
            fullpath = r"{}".replace('\\\\', '\\').replace('/', '\\').format(fil).replace('/', '\\')
            allnames += f" '{fullpath}'"
    # if 'parents_sample' in toplot:
    #     for fil in m.get_paget_absolute_fullpaths(buf=True):
    #         fullpath = r"{}".replace('\\\\', '\\').replace('/', '\\').format(fil).replace('/', '\\')
    #         allnames += f" '{fullpath}'"
    out = run_program('primusqt', allnames)
    return out



def datgnom(m, r=20):
    infile = r"{}".replace('\\\\', '\\').format(m.get_filename()).replace('/', '\\')
    outfile = r"{}".replace('\\\\', '\\').replace('/', '\\').format(join(m.path, 'datgnom.out')).replace('/', '\\')
    pofrfile = r"{}".replace('\\\\', '\\').replace('/', '\\').format(join(m.path, 'ou2pofr.dat')).replace('/', '\\')
    fitfile = r"{}".replace('\\\\', '\\').replace('/', '\\').format(join(m.path, 'out2fit.dat')).replace('/', '\\')
    out = run_program('datgnom', f"'{infile}'", '-o', f"'{outfile}'", '-r', f'{r}')
    out = run_program('out2pofr', f"'{outfile}'", '-o', f"'{pofrfile}'")
    out = run_program('out2fit', f"'{outfile}'", '-o', f"'{fitfile}'")
    print(out)
    return out

def datmif(m, r=20):
    infile = r"{}".replace('\\\\', '\\').format(m.get_filename()).replace('/', '\\')
    outfile = r"{}".replace('\\\\', '\\').replace('/', '\\').format(join(m.path, 'datmif.out')).replace('/', '\\')
    # pofrfile = r"{}".replace('\\\\', '\\').replace('/', '\\').format(join(m.path, 'ou2pofr.dat')).replace('/', '\\')
    # fitfile = r"{}".replace('\\\\', '\\').replace('/', '\\').format(join(m.path, 'out2fit.dat')).replace('/', '\\')
    # out = run_program('datmif', f"'{infile}'", '-o', f"'{outfile}'", '-r', f'{r}')
    out = run_program('datmif', f"'{infile}'", '-o', f"'{outfile}'", '-r', f'{r}', '--first=1', '--last=2300', '--dmax=140')
    # out = run_program('out2pofr', f"'{outfile}'", '-o', f"'{pofrfile}'")
    # out = run_program('out2fit', f"'{outfile}'", '-o', f"'{fitfile}'")
    print(out)
    return out

def get_out2fit(m):
    fitfile = r"{}".replace('\\\\', '\\').replace('/', '\\').format(join(m.path, 'out2fit.dat')).replace('/', '\\')
    # df = pd.DataFrame(get_outpath(m, 'out2fit.dat'))
    df = pd.DataFrame(fitfile)
    return df

def gnom(m, r=20):
    infile = r"{}".replace('\\\\', '\\').format(m.get_filename()).replace('/', '\\')
    outfile = r"{}".replace('\\\\', '\\').replace('/', '\\').format(join(m.path, 'gnom.out')).replace('/', '\\')
    out = run_program('gnom', f"'{infile}'", '-o', f"'{outfile}'", '-r', f'{r}')
    print(out)
    return out

def plot_datgnom(m):
    outfile = r"{}".replace('\\\\', '\\').replace('/', '\\').format(join(m.path, 'datgnom.out')).replace('/', '\\')
    out = run_program('gnomplotqt', f"'{outfile}'")
    # for line in open(outfile):
    #     print(line)
    return out

def print_datgnom(m):
    outfile = r"{}".replace('\\\\', '\\').replace('/', '\\').format(join(m.path, 'datgnom.out')).replace('/', '\\')
    for line in open(outfile):
        print(line)

def autorg(m):
    infile = r"{}".replace('\\\\', '\\').format(m.get_filename()).replace('/', '\\')
    outfile = r"{}".replace('\\\\', '\\').replace('/', '\\').format(join(m.path, 'autorg.csv')).replace('/', '\\')
    out = run_program('autorg', f"'{infile}'", '-o', f"'{outfile}'")
    print(out)
    return out

def bodies(m, rg):
    infile = r"{}".replace('\\\\', '\\').format(m.get_filename()).replace('/', '\\')
    outfile = r"{}".replace('\\\\', '\\').replace('/', '\\').format(join(m.path, 'bodies.csv')).replace('/', '\\')
    # outfile = r"{}".replace('\\\\', '\\').replace('/', '\\').format(join(m.path, 'datgnom.out')).replace('/', '\\')
    out = run_program('bodies', f"'{infile}'", '-o', f"'{outfile}'", '--rg=', f"'{rg}'")
    print(out)
    return out

def dammin(m, rg):
    """This one would need a GNOMFILE"""
    out = run_program('dammin')
    return out

def datft(m, rg):
    """--dmax is in A
    --particle-bead-radius is in A
    it takes forever and easyly gives some integer overflow error."""
    out = run_program('datft')
    return out
