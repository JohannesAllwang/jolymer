import subprocess
import getpass
import io
import os
from os.path import join

atsas_bin = f"C:\\Users\\{getpass.getuser()}\\atsas\\bin"

def run_program(name, *args):
    cmd = ['powershell', f'./{name}.exe', *args]
    return subprocess.run(cmd, cwd=atsas_bin)

def datgnom(m, r=20):
    infile = r"{}".replace('\\\\', '\\').format(m.get_filename()).replace('/', '\\')
    outfile = r"{}".replace('\\\\', '\\').replace('/', '\\').format(join(m.path, 'datgnom.out')).replace('/', '\\')
    out = run_program('datgnom', f"'{infile}'", '-o', f"'{outfile}'", '-r', f'{r}')
    print(out)
    return out

def autorg(m):
    infile = r"{}".replace('\\\\', '\\').format(m.get_filename()).replace('/', '\\')
    outfile = r"{}".replace('\\\\', '\\').replace('/', '\\').format(join(m.path, 'autorg.csv')).replace('/', '\\')
    out = run_program('autorg', f"'{infile}'", '-o', f"'{outfile}'")
    print(out)
    return out
