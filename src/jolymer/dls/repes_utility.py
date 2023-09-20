import subprocess
import getpass
import io
import os
from os.path import join
import pandas as pd
import numpy as np

from .cun import Cun

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


def write_ALV5000(m, dfg2, dfI=None, filename=None,
                  encoding='utf_16_be'):
    if filename is None:
        filename = join(m.phidls_path, 'test.ASC')
    print(filename)
    # encoding = 'utf_16_be'
    # encoding = 'utf_16_le'
    with open('/home/johannes/jodls/sample6000.asc', 'r',
              encoding=encoding) as sample:
        with open(filename, 'w',
                  encoding=encoding) as f:
            for line in sample:
                f.write(line)


def write_ALV6000(m, dfg2, dfI=None, filename=None):
    if filename is None:
        filename = join(m.phidls_path, 'test.ASC')
    print(filename)
    print('TODO: Make nice metadata!')
    with io.open(filename, 'w', encoding='utf_16_be') as f:
    # with io.open(filename, 'w', encoding='utf-8') as f:
        # with io.open('/home/johannes/jodls/sample6000.asc', 'r',
        #           encoding='utf_16_be') as sample_file:
            # for line in sample_file:
            #     f.write(line)
                # break
        # f.write(u'test\r\n')
        f.write('ALV-6000/E-WIN Data')
        f.write(u'\r\n')
        f.write('Date :	"5.13.2017"')
        f.write(u'\r\n')
        f.write('Time :	"20:31:39"')
        f.write(u'\r\n')
        f.write(f'Samplename : 	"{m.samplestring}"')
        f.write(u'\r\n')
        f.write('SampMemo(0) : 	""')
        f.write(u'\r\n')
        f.write('SampMemo(1) : 	""')
        f.write(u'\r\n')
        f.write('SampMemo(2) : 	""')
        f.write(u'\r\n')
        f.write('SampMemo(3) : 	""')
        f.write(u'\r\n')
        f.write('SampMemo(4) : 	""')
        f.write(u'\r\n')
        f.write('SampMemo(5) : 	""')
        f.write(u'\r\n')
        f.write('SampMemo(6) : 	""')
        f.write(u'\r\n')
        f.write('SampMemo(7) : 	""')
        f.write(u'\r\n')
        f.write('SampMemo(8) : 	""')
        f.write(u'\r\n')
        f.write('SampMemo(9) :   ""')
        f.write(u'\r\n')
        f.write('Temperature [K] :	     297.99188')
        f.write(u'\r\n')
        f.write('Viscosity [cp]  :	       0.89343')
        f.write(u'\r\n')
        f.write('Refractive Index:	       1.37200')
        f.write(u'\r\n')
        f.write('Wavelength [nm] :	     632.79999')
        f.write(u'\r\n')
        f.write('Angle [°]       :	      39.99900')
        f.write(u'\r\n')
        f.write('Duration [s]    :	       120')
        f.write(u'\r\n')
        f.write('Runs            :	         1')
        f.write(u'\r\n')
        f.write('Mode            :	"SINGLE CROSS CH0 "')
        f.write(u'\r\n')
        f.write('MeanCR0 [kHz]   :	       2.92202')
        f.write(u'\r\n')
        f.write('MeanCR1 [kHz]   :	       3.58007')
        f.write(u'\r\n')
        f.write(u'\r\n')

        print('repes_utility.write_ALV gives the viscosity of water to the file...')
        print('Also time and date are arbitrary')
        f.write('"Correlation"')
        f.write(u'\r\n')
        for idx, row in dfg2.iterrows():
            f.write('  {:.5E}\t  {:.5E}'.format(
                row.t, row.g2
                ))
            # f.write('  1.25000E-004	  6.21406E-001')
            f.write('\r\n')
        f.write('\r\n')
        f.write('"Count Rate"')
        f.write('\r\n')
        for idx, row in dfg2.iterrows():
            CRA = np.random.random()
            CRB = np.random.random()
            f.write(' {:.5f}\t {:.5f}\t {:.5f}'.format(
                row.t, CRA, CRB
                ))
            f.write('\r\n')
        f.write('\r\n')
        f.write('Monitor Diode	 414198.98\n')
        f.write('\r\n')
        f.write('"Cumulant 1.Order"')
        f.write('\r\n')
        f.write('FluctuationFreq. [1/ms]	 8.3350E-002')
        f.write('\r\n')
        f.write('DiffCoefficient [µm²/s]	 9.5990E-001')
        f.write('\r\n')
        f.write('Hydrodyn. Radius [nm]	 2.5439E+002')
        f.write('\r\n')
        f.write('\r\n')
        f.write('"Cumulant 2.Order"')
        f.write('\r\n')
        f.write('FluctuationFreq. [1/ms]	 1.4396E-001')
        f.write('\r\n')
        f.write('DiffCoefficient [µm²/s]	 1.6579E+000')
        f.write('\r\n')
        f.write('Hydrodyn. Radius [nm]	 1.4728E+002')
        f.write('\r\n')
        f.write('Expansion Parameter µ2	 1.2359E-002')
        f.write('\r\n')
        f.write('\r\n')
        f.write('"Cumulant 3.Order"')
        f.write('\r\n')
        f.write('FluctuationFreq. [1/ms]	 2.0369E-001')
        f.write('\r\n')
        f.write('DiffCoefficient [µm²/s]	 2.3458E+000')
        f.write('\r\n')
        f.write('Hydrodyn. Radius [nm]	 1.0410E+002')
        f.write('\r\n')
        f.write('Expansion Parameter µ2	 4.5287E-002')
        f.write('\r\n')
        f.write('Expansion Parameter µ3	 5.4684E-003')
        f.write('\r\n')
        f.write('\r\n')

        f.write('"StandardDeviation"\n')
        for idx, row in dfg2.iterrows():
            f.write('  {:.5E}\t  {:.5E}'.format(
                row.t, row.err_g2
                ))
            f.write('\r\n')
        f.write('\r\n')
        f.write('"Special"')
        f.write('\r\n')
        f.write('Temperature Server [K] :	     298.1600')
        f.write('\r\n')
        f.write('Temperature Server [C] :	     25.0100')
        f.write('\r\n')


def change_pamet(path, errors='exper.', par1=0, par2=5, par3=1000):
    if getpass.getuser() == 'johannes':
        with open(join(gendist_bin, 'pamet.pam'), 'w',
                  encoding='utf_16_be') as f:
            f.write('z:\n')
            f.write(path + '\n')
            f.write(path + '\n')
            f.write('ALV-6000/E\n')
            f.write(str(errors) + '\n')
            f.write(' 0 \n')
            f.write('5\n')
            f.write('1000\n')
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


def run_gendist(path):
    # change_pamet(path)
    cmd = ['powershell', './gen_M.exe']
    if getpass.getuser() == 'johannes':
        cmd = ['wine', 'start', './gen_M.exe']
    return subprocess.run(cmd, cwd=gendist_bin)


class REPES(Cun):

    def __init__(self):
        self.name = 'repes'

    def get_moA(self, m, phi, **kwargs):
        return m.get_moA(phi, **kwargs)

    def get_average_tau(self, m, phi, mean='log', rmin=0, rmax=np.inf, A='A'):
        df = m.get_Arl(phi, rmin=rmin, rmax=rmax, A=A)
        if mean == 'log':
            out = np.sum(df.dist * np.log(df.t)) / np.sum(df.dist)
            out = np.exp(out)
        elif mean == 'arithmetic':
            out = np.sum(df.dist * df.t)
        return out

    def get_rh(self, m, phi, mean='log', rmin=0, rmax=np.inf, A='A'):
        df = m.get_Arl(phi, rmin=rmin, rmax=rmax, A=A)
        if mean == 'log':
            out = np.sum(df.dist * np.log(df.rapp)) / np.sum(df.dist)
            out = np.exp(out)
        elif mean == 'arithmetic':
            out = np.sum(df.dist * df.rapp)
        return out

    def get_phidlstable(self, m, rmin=0, rmax=np.inf):
        tempdict = {'phi': [],
                    'qq': [],
                    'Dapp': [],
                    'Rapp': [],
                    'tau': [],
                    'beta': [],
                    'SqBeta': [],
                    'Baseline': [],
                    'Gamma': []}
        for phi in m.angles:
            qq = m.qq(phi)
            tau = self.get_average_tau(m, phi, rmin=rmin, rmax=rmax)
            resultdict = self.get_moA(m, phi)
            SqBeta = resultdict['SqBeta']
            Baseline = resultdict['Baseline']
            tempdict['phi'].append(phi)
            tempdict['qq'].append(qq)
            tempdict['Dapp'].append(1 / (qq * tau))
            tempdict['Rapp'].append(m.Rfromt(tau, qq))
            tempdict['tau'].append(tau)
            tempdict['Gamma'].append(1/tau)
            tempdict['SqBeta'].append(SqBeta)
            tempdict['beta'].append(SqBeta ** 2)
            tempdict['Baseline'].append(Baseline)
        df = pd.DataFrame(data=tempdict)
        return df



repes = REPES()
