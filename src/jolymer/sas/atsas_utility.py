import subprocess
import getpass
import io
import os
from os.path import join

atsas_bin = f"C:\\Users\\{getpass.getuser()}\\atsas\\bin"

def run_program(name, *args):
    cmd = ['powershell', f'./{name}.exe', *args]
    return subprocess.run(cmd, cwd=atsas_bin)

def primus(m, fname=''):
    fullpath = r"{}".replace('\\\\', '\\').replace('/', '\\').format(join(m.path, fname)).replace('/', '\\')
    if fname=='':
        fullpath = r"{}".replace('\\\\', '\\').replace('/', '\\').format(join(m.get_filename())).replace('/', '\\')
    out = run_program('primusqt', f"'{fullpath}'")
    # for line in open(outfile):
    #     print(line)
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
