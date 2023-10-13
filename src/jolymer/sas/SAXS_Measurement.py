"""
"""

from scipy import optimize
import numpy as np
import re
from os.path import join
from matplotlib.colors import LogNorm
from pylab import cm
import matplotlib.pyplot as plt
import pandas as pd

import sasmodels
from sasmodels import data as sasmodels_data
import pyFAI
import jscatter as js
import fabio

from .. import database_operations as dbo
from ..Measurement import Measurement

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
class SAXS_Measurement(Measurement):

    instrument = 'no instrument'

    rawpath: str= ''
    path: str = '~'
    filename: str='merge_001.dat'

    def get_filename(self):
        return join(self.path, self.filename)

    def get_data(self, engine='pandas', **kwargs):
        path = self.get_filename()
        self.data1d = sasmodels_data.load_data(path, **kwargs)
        if engine == 'sasview':
            return self.data1d
        elif engine == 'pandas':
            outdict = {'q': self.data1d.x,
                       'I': self.data1d.y,
                       'err_I': self.data1d.dy}
            return pd.DataFrame(outdict)

    def get_filename_sasImage(self):
        with open(self.get_filename(), 'r') as f:
            lines = f.readlines()
        matched_lines = [line for line in lines if re.search('Sample filename', line)]
        filename = matched_lines[0].split()[0]
        return join(self.rawpath, filename)

    def get_sasImage(self):
        return fabio.open(self.get_filename_sasImage()).data


    def pyfai_integrate1d(self):
        masked_image = self.get_masked()
        nbins=200.
        sdd = self.sasImage.detector_distance[0]
        centerx, centery = masked_image.center
        pixelsizex, pixelsizey = masked_image.pixel_size
        ai = pyFAI.azimuthalIntegrator.AzimuthalIntegrator(dist=sdd,
                                                           poni1=centerx*pixelsizex,
                                                           poni2=centery*pixelsizey,
                                                           detector='pilatus300k')
        # ai.setFit2D(sdd, centerx, centery)
        ai.wavelength = masked_image.wavelength[0] * 10**-10

        q, I, err_I = ai.integrate1d(data=masked_image.data, npt=nbins,
                                     unit="q_nm^-1", mask=masked_image.mask,
                                     error_model='poisson', correctSolidAngle=True)
        dict={'q':q, 'I':I/self.exposure_time, 'err_I':err_I/self.exposure_time}
        df = pd.DataFrame(dict)
        return df

    def show_sasImage(self, **kwargs):
        fig, ax = plt.subplots()
        img = self.get_sasImage()
        im = ax.matshow(img, cmap=cm.viridis, origin='lower',
                        norm=LogNorm(vmin=0.01, vmax=10000), **kwargs)
        _colorbar(im)

def gen_guinier_fitfunc(alpha):
    def inner(q, Rg, A):
        if alpha == 0:
            pre = 1
        elif alpha == 1 or alpha == 2:
            pre = alpha * np.pi * q ** -alpha
        else:
            raise TypeError('alpha needs to be in 0,1,2')

        I = pre * A * np.exp(-Rg ** 2 * q ** 2 / (3 - alpha))
        return I
    return inner

def guinier_porod_3D(q, Rg1, s1, Rg2, s2, G2, dd):
    q = np.atleast_1d(q)

    # define parameters for smooth transitions
    Q1 = (1 / Rg1) * ((dd - s1) * (3 - s1) / 2) ** 0.5
    Q2 = ((s1 - s2) / (2 / (3 - s2) * Rg2 ** 2 - 2 / (3 - s1) * Rg1 ** 2)) ** 0.5
    G1 = G2 / (np.exp(-Q2 ** 2 * (Rg1 ** 2 / (3 - s1) -
                                  Rg2 ** 2 / (3 - s2))) * Q2 ** (s2 - s1))
    D = G1 * np.exp(-Q1 ** 2 * Rg1 ** 2 / (3 - s1)) * Q1 ** (dd - s1)

    # define functions in different regions
    def _I1_3regions(q):
        res = G2 / q ** s2 * np.exp(-q ** 2 * Rg2 ** 2 / (3 - s2))
        return res

    def _I2_3regions(q):
        res = G1 / q ** s1 * np.exp(-q ** 2 * Rg1 ** 2 / (3 - s1))
        return res

    def _I3_3regions(q):
        res = D / q ** dd
        return res

    I = np.piecewise(q, [q < Q2, (Q2 <= q) & (q < Q1), q >= Q1],
                     [_I1_3regions, _I2_3regions, _I3_3regions])
    return I

