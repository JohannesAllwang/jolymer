import subprocess
import getpass
import io
import os
from os.path import join
import pandas as pd

gendist_bin = f"C:\\Users\\{getpass.getuser()}\\repes\\gen11\\Genr11\\Basic0.0"

def get_outpath(m, name):
    outfile = r"{}".replace('\\\\', '\\').replace('/', '\\').format(join(m.path, name)).replace('/', '\\')
    return outfile

def run_program(name, *args):
    cmd = ['powershell', f'./{name}.exe', *args]
    return subprocess.run(cmd, cwd=gendist_bin)

def create_repes_dir(m):
    pass

def make_average_ALV(m):
    pass

def get_mo(m, phi, A='A'):
    pass

def write_ALV(m, dfg2, dfI=None):
    with open(filename) as f:
        write('ALV-5000/E-WIN Data\n')
        write(f'Date :	"02.03.20"\n')
        write(f'Time :	"10:24 AM"\n')
        write(f'Samplename : 	""\n')
        write('SampMemo(0) : 	""\n')
        write('SampMemo(1) : 	""\n')
        write('SampMemo(2) : 	""\n')
        write('SampMemo(3) : 	""\n')
        write('SampMemo(4) : 	""\n')
        write('SampMemo(5) : 	""\n')
        write('SampMemo(6) : 	""\n')
        write('SampMemo(7) : 	""\n')
        write('SampMemo(8) : 	""\n')
        write('SampMemo(9) : 	""\n')
        write(f'Temperature [K] :	     293.20000\n')
        write(f'Viscosity [cp]  :	       1.15900\n')
        write('\n')
        write('"Correlation"\n')
        for row in dfg2:
            write('  {}\t  {}\t  {}'.format(
                q, g2, err_g2
                ))
        write('\n')
        write('"Count Rate"\n')
        if not dfI==None:
            for row in dfI:
                write('  {}\t  {}\t  {}'.format(
                t, CRA, CRB
                ))

def change_pamet(path, errors='exper.', par1=0, par2=5, par3=1000):
    with open(join(gendist_bin, 'pamet.pam'), 'w') as f:
        f.write('c: [OSDisk]' + '\n')
        f.write(path + '\n')
        f.write(path + '\n')
        f.write('ALV-5000/E' + '\n')
        f.write(str(errors) + '\n')
        f.write(f' {par1} \n')
        f.write(str(par2) + '\n')
        f.write(str(par3) + '\n')


def run_gendist(path):
    change_pamet(path)
    out = run_program('gen_M')
    return out
