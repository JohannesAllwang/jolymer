# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 15:12:01 2020

@author: xcill
"""

import pyFAI
import matplotlib.pyplot as plt
from .. import database_operations as dbo
import jscatter as js
import pandas as pd
import os
from .. import Sample
from ..Measurement import Measurement



class _goc(Measurement):
    """ganesha measurement at only configuration 2 or 3"""
    def __init__(self, path, filename):
        self.path = path
        self.filename_tiff = path + filename + '.tiff'
        self.filename_masked = path + filename + '_masked.tiff'
        self.sasImage = js.sas.sasImage(self.filename_tiff)
        self.exposure_time = self.sasImage.exposure_time[0]            
        self.transmission_factor = self.sasImage.transmission_factor[0]
        self.config = filename[-1]
        self.filename_pyfai = f'{path}{filename}_pyfai.csv'
        self.filename_fit2d = f'{path}{filename}_fit2d'
    
    def get_pyfai(self):
        df = pd.read_csv(self.filename_pyfai)
        return df

    def get_fit2d(self, iqmin=0, iqmax=-1):
        skiprows = 5
        df =  pd.read_csv(self.filename_fit2d, skiprows = skiprows,
                            delimiter = '  ', header=None, names = ['q', 'I'])[iqmin:iqmax]
        df['Ibyt'] = df.I / self.exposure_time
        return df

    def get_masked(self):
        return js.sas.sasImage(self.filename_masked)

    def plot_masked(self, **kwargs):
        masked = self.get_masked()
        figure = masked.show(scale = 'symlog', colorMap = 'ocean', badcolor = 'red', axis = 'pixel', **kwargs)

    def pyfai_integrate1d(self):
        masked_image = self.get_masked()
        nbins=200.
        sdd = self.sasImage.detector_distance[0]
        centerx, centery = masked_image.center
        pixelsizex, pixelsizey = masked_image.pixel_size
        ai = pyFAI.azimuthalIntegrator.AzimuthalIntegrator(dist=sdd, poni1=centerx*pixelsizex, poni2=centery*pixelsizey, detector='pilatus300k')
        # ai.setFit2D(sdd, centerx, centery)
        ai.wavelength = masked_image.wavelength[0] * 10**-10

        q, I, err_I = ai.integrate1d(data=masked_image.data, npt=nbins, unit="q_nm^-1", mask=masked_image.mask, error_model='poisson', correctSolidAngle=True)
        dict={'q':q, 'I':I/self.exposure_time, 'err_I':err_I/self.exposure_time}
        df = pd.DataFrame(dict)
        return df
    
    def get_Iuncorrected(self, iqmin=0, iqmax=-1):
        filename = os.path.join(self.path, f'uncorrected{self.config}.csv')
        df = pd.read_csv(filename)[iqmin:iqmax]
        df['Ibyt'] = df.I / self.exposure_time
        return df
    
    def get_Icorrected(self, bg, dark, iqmin = 8):
        Ts = self.transmission_factor
        Tbg = bg.transmission_factor
        dfsample = self.get_Iuncorrected()
        dfdark = dark.get_Iuncorrected()
        dfbg = bg.get_Iuncorrected()
        
        A = dfsample.Ibyt - dfdark.Ibyt
        B = dfbg.Ibyt - dfdark.Ibyt
        B = B * Ts / Tbg
        Icorr = A - B
        frame = {'q': dfsample.q, 'I': Icorr}
        out = pd.DataFrame(frame)[0:-1]
        return out
    
    def plot_Iuncorrected(self, fit2d = False, label = None, scale = 1, iqmin=0, iqmax=-1, **kwargs):
        df = self.get_Iuncorrected(iqmin=iqmin, iqmax=iqmax)
        if fit2d:
            df = self.get_fit2d(iqmin=iqmin, iqmax=iqmax)
        if 'figure' in kwargs:
            fig, ax = kwargs['figure']
            kwargs.pop('figure')

        else:
            fig, ax = plt.subplots()
            ax.set_xlabel('$q$ [1/nm]')
            ax.set_ylabel('$I$ [1/cm]')
        if 'shift' in kwargs:
            scale = kwargs['shift']
            kwargs.pop('shift')
        if 'marker' in kwargs:
            marker=kwargs['marker']
            kwargs.pop('marker')
        else:
            marker='.'
        if 'linestyle' in kwargs:
            linestyle=kwargs['linestyle']
            kwargs.pop('linestyle')
        else:
            linestyle=''

        markers, caps, bars = ax.errorbar(df.q, df.Ibyt * scale, marker = marker, 
                    linestyle=linestyle, label = label, elinewidth=0.2,  **kwargs)
        # [bar.set_alpha(0.2) for bar in bars]
        
        ax.set_xscale('log')
        ax.set_yscale('log')
        return fig, ax

class Ganesha(Measurement):

    def __init__(self, ganesha_id, **kwargs):
        with dbo.dbopen() as c:
            c.execute("SELECT * FROM ganesha_measurements WHERE id=?;", (ganesha_id,))
            self.id, self.datestring, self.samplestring, self.comment = list(c.fetchone())
        print(2)
        print(self.samplestring)
        if self.samplestring == None:
            self.sample = None
        else:
            sample_type, sample_id = self.samplestring.split('_')
            self.sample = Sample.get_sample(sample_type, sample_id)
        self.path = os.path.join(self.rawdatapath, 'ganesha{0:03}/'.format(self.id))
        # get the filenames for configuration 2 and 3:
        self.filename2 = 'conf2'
        self.filename3 = 'conf3'
        self.config2 = _goc(self.path, self.filename2)
        self.config3 = _goc(self.path, self.filename3)
        
    def get_config(self, n):
        if n==2:
            return self.config2
        elif n==3:
            return self.config3
    
    def get_bg(self):
        buffer = self.sample.buffer
        bg = Ganesha(buffer.ganesha_id)
        return bg
    
    def get_dark(self):
        return Ganesha(0)
    
    def get_data(self):
        df = pd.read_csv(self.path + 'Ifinal.csv')
        return df
    
    def plot_data(self, label = None, scale = 1, **kwargs):
        df = self.get_data()
        if 'figure' in kwargs:
            fig, ax = kwargs['figure']
            kwargs.pop('figure')

        else:
            fig, ax = plt.subplots()
            ax.set_xlabel('$q$ [1/nm]')
            ax.set_ylabel('$I$ [1/cm]')
        if 'shift' in kwargs:
            scale = kwargs['shift']
            kwargs.pop('shift')
        if 'marker' in kwargs:
            marker=kwargs['marker']
            kwargs.pop('marker')
        else:
            marker='.'
        if 'linestyle' in kwargs:
            linestyle=kwargs['linestyle']
            kwargs.pop('linestyle')
        else:
            linestyle=''

        markers, caps, bars = ax.errorbar(df.q, df.I * scale, marker = marker, 
                    linestyle=linestyle, label = label, elinewidth=0.2,  **kwargs)
        # [bar.set_alpha(0.2) for bar in bars]
        
        ax.set_xscale('log')
        ax.set_yscale('log')
        return fig, ax
        

        
        

        
        
