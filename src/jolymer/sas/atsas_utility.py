import subprocess
import getpass
import io
import os
from os.path import join
import pandas as pd

atsas_bin = f"C:\\Users\\{getpass.getuser()}\\atsas\\bin"

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
