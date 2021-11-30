import subprocess
import getpass
import io
import os
from os.path import join
import pandas as pd

gendist_bin = f"C:\\Users\\{getpass.getuser()}\\repes\\gen11\\Genr11\\Basic0.0"
testpath = f"C:\\Users\\{getpass.getuser()}\\repes\\testdata\\"
testpath1 = join(testpath, 'Jan20_P123_01')

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
    filename = join(m.phidls_path, 'test.ASC')
    with open(filename, 'w') as f:
        f.write('ALV-5000/E-WIN Data\n')
        f.write(f'Date :	"00.00.00"\n')
        f.write(f'Time :	"00:00 AM"\n')
        f.write(f'Samplename : 	""\n')
        f.write('SampMemo(0) : 	""\n')
        f.write('SampMemo(1) : 	""\n')
        f.write('SampMemo(2) : 	""\n')
        f.write('SampMemo(3) : 	""\n')
        f.write('SampMemo(4) : 	""\n')
        f.write('SampMemo(5) : 	""\n')
        f.write('SampMemo(6) : 	""\n')
        f.write('SampMemo(7) : 	""\n')
        f.write('SampMemo(8) : 	""\n')
        f.write('SampMemo(9) : 	""\n')
        f.write('Temperature [K] :	     {:.5f}\n'.format(m.get_TK()))
        f.write('Viscosity [cp]  :	       {:.5f}\n'.format(1.15900))
        f.write('Refractive Index:	       1.36100\n')
        f.write('Wavelength [nm] :	     632.80000\n')
        f.write('Angle [°]       :	      30.00000\n')
        f.write('Duration [s]    :	        30.000000\n')
        f.write('Runs            :	         1\n')
        f.write('Mode            :	"Pseudo Cross Correlation"\n')
        f.write('MeanCR0 [kHz]   :	      141.50000\n')
        f.write('MeanCR1 [kHz]   :	      133.20000\n')

        print('repes_utility.write_ALV gives the viscosity of water to the file...')
        print('Also time and date are arbitrary')
        f.write('\n')
        f.write('"Correlation"\n')
        for idx, row in dfg2.iterrows():
            f.write('  {:.5f}\t  {:.5f}\t  {:.5f}'.format(
                row.t, row.g2, row.err_g2
                ))
            f.write('\n')
        f.write('\n')
        # f.write('"Count Rate"\n')
        # if not dfI==None:
            # for row in dfI:
            #     f.write('  {}\t  {}\t  {}'.format(
                # t, CRA, CRB
                # ))

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
