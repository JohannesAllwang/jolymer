import subprocess
import getpass
import io
import os
from os.path import join
import pandas as pd

gendist_bin = f"C:\\Users\\{getpass.getuser()}\\repes\\gen11\\Genr11\\Basic0.0"
testpath = f"C:\\Users\\{getpass.getuser()}\\repes\\testdata\\"
testpath1 = join(testpath, 'Jan20_P123_01')
if getpass.getuser() == 'johannes':
    gendist_bin = '/home/johannes/repes/gen11/Genr11/Basic0.0/'
    testpath = "~/repes/testdata/"
    testpath1 = join(testpath, 'Jan20_P123_01')


def get_outpath(m, name):
    outfile = r"{}".replace('\\\\', '\\').replace('/', '\\').format(join(m.path, name)).replace('/', '\\')
    return outfile


def create_repes_dir(m):
    pass


def make_average_ALV(m):
    pass


def get_mo(m, phi, A='A'):
    pass

def write_ALV5000(m, dfg2, dfI=None, filename=None):
    if filename is None:
        filename = join(m.phidls_path, 'test.ASC')
    print(filename)
    with open(filename, 'w') as f:
        f.write('ALV-5000/E-WIN Data\n')
        f.write('Date :	"5.12.2017"\n')
        f.write('Time :  "00:00:00"\n')
        f.write('Samplename : 	"lol"\n')
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
        f.write('''Temperature [K] :	     297.99188
Viscosity [cp]  :	       0.89343
Refractive Index:	       1.37200
Wavelength [nm] :	     632.79999
Angle [°]       :	      39.99900
Duration [s]    :	       120
Runs            :	         1
Mode            :	"SINGLE CROSS CH0 "
MeanCR0 [kHz]   :	       2.92202
MeanCR1 [kHz]   :	       3.58007
''')
        # f.write('Temperature [K] :	     {:.5f}\n'.format(m.get_TK()))
        # f.write('Viscosity [cp]  :	       {:.5f}\n'.format(1.15900))
        # f.write('Refractive Index:	       1.36100\n')
        # f.write('Wavelength [nm] :	     632.80000\n')
        # f.write('Angle [°]       :	      30.00000\n')
        # f.write('Duration [s]    :	        30\n')
        # f.write('Runs            :	         1\n')
        # f.write('Mode            :	"mod3D"\n')
        # f.write('MeanCR0 [kHz]   :	      141.50000\n')
        # f.write('MeanCR1 [kHz]   :	      133.20000\n')


def write_ALV6000(m, dfg2, dfI=None, filename=None):
    if filename is None:
        filename = join(m.phidls_path, 'test.ASC')
    print(filename)
    with open(filename, 'w') as f:
        with open('/home/johannes/jodls/sample6000.asc', 'r') as sample:
            for line in sample:
                print(line)
                # f.write(line)
    return 0
    with open(filename, 'w') as f:
        f.write('ALV-6000/E-WIN Data\n')
        f.write('Date :	"5.12.2017"\n')
        f.write('Time :  "00:00:00"\n')
        f.write('Samplename : 	"lol"\n')
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
        f.write('''Temperature [K] :	     297.99188
Viscosity [cp]  :	       0.89343
Refractive Index:	       1.37200
Wavelength [nm] :	     632.79999
Angle [°]       :	      39.99900
Duration [s]    :	       120
Runs            :	         1
Mode            :	"SINGLE CROSS CH0 "
MeanCR0 [kHz]   :	       2.92202
MeanCR1 [kHz]   :	       3.58007
''')
        # f.write('Temperature [K] :	     {:.5f}\n'.format(m.get_TK()))
        # f.write('Viscosity [cp]  :	       {:.5f}\n'.format(1.15900))
        # f.write('Refractive Index:	       1.36100\n')
        # f.write('Wavelength [nm] :	     632.80000\n')
        # f.write('Angle [°]       :	      30.00000\n')
        # f.write('Duration [s]    :	        30\n')
        # f.write('Runs            :	         1\n')
        # f.write('Mode            :	"mod3D"\n')
        # f.write('MeanCR0 [kHz]   :	      141.50000\n')
        # f.write('MeanCR1 [kHz]   :	      133.20000\n')

        print('repes_utility.write_ALV gives the viscosity of water to the file...')
        print('Also time and date are arbitrary')
        f.write('\n')
        f.write('"Correlation"\n')
        for idx, row in dfg2.iterrows():
            # f.write('  {:.5E}\t  {:.5E}\n'.format(
            #     row.t, row.g2
            #     ))
            f.write('  1.25000E-004	  6.21406E-001\n')
        f.write('\n')
        return 0
        f.write('"Count Rate"\n')
        for idx, row in dfg2.iterrows():
            f.write(' {:.5f}\t {:.5f}\t {:.5f}'.format(
                row.t, 1, 1
                ))
            f.write('\n')
        f.write('\n')
        f.write('Monitor Diode	 414198.98\n')
        # f.write('\n')
        f.write('''
"Cumulant 1.Order"
FluctuationFreq. [1/ms]	 8.3350E-002
DiffCoefficient [µm²/s]	 9.5990E-001
Hydrodyn. Radius [nm]	 2.5439E+002

"Cumulant 2.Order"
FluctuationFreq. [1/ms]	 1.4396E-001
DiffCoefficient [µm²/s]	 1.6579E+000
Hydrodyn. Radius [nm]	 1.4728E+002
Expansion Parameter µ2	 1.2359E-002

"Cumulant 3.Order"
FluctuationFreq. [1/ms]	 2.0369E-001
DiffCoefficient [µm²/s]	 2.3458E+000
Hydrodyn. Radius [nm]	 1.0410E+002
Expansion Parameter µ2	 4.5287E-002
Expansion Parameter µ3	 5.4684E-003

''')

        f.write('"StandardDeviation"\n')
        for idx, row in dfg2.iterrows():
            f.write('  {:.5E}\t  {:.5E}'.format(
                row.t, row.err_g2
                ))
            f.write('\n')
        f.write('''
"Special"
  Temperature Server [K] :	     298.1600
  Temperature Server [C] :	     25.0100''')


def change_pamet(path, errors='exper.', par1=0, par2=5, par3=1000):
    out = []
    with open(join(gendist_bin, 'pamet.pam'), 'r') as f:
        for line in f:
            out.append(line)
    if getpass.getuser() == 'johannes':
        with open(join(gendist_bin, 'pamet.pam'), 'w') as f:
            f.write('z:\n')
            f.write(path + '\n')
            f.write(path + '\n')
            f.write('ALV-6000/E\n')
            f.write(str(errors) + '\n')
            f.write('0\n')
            f.write('5\n')
            f.write('1000')
            # f.write(str(par2) + '\n')
            # f.write(str(par3) + '\n')
    else:
        with open(join(gendist_bin, 'pamet.pam'), 'w') as f:
            f.write('c: [OSDisk]' + '\n')
            f.write(path + '\n')
            f.write(path + '\n')
            f.write('ALV-6000/E' + '\n')
            f.write(str(errors) + '\n')
            f.write(f' {par1} \n')
            f.write(str(par2) + '\n')
            f.write(str(par3) + '\n')
    return out


def run_gendist(path):
    change_pamet(path)
    cmd = ['powershell', './gen_M.exe']
    if getpass.getuser() == 'johannes':
        cmd = ['wine', 'start', './gen_M.exe']
    return subprocess.run(cmd, cwd=gendist_bin)
