#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 12:32:47 2020

@author: johannes
"""

from . import database_operations as dbo
# from .tresy_samples import *

import pandas as pd
import numpy as np
import datetime as dt
from scipy import optimize, constants
from pathlib import Path
from dataclasses import dataclass

def n_water(T, wl, rho=1000):
    """
    Thormählen, I. / Straub, J. / Grigull, U.
    Refractive Index of Water and Its Dependence on Wavelength, Temperature, and Density
    1985-10

    Journal of Physical and Chemical Reference Data , Vol. 14, No. 4
    AIP Publishing
    p. 933-945

    Args:
    rho : density in kg/m^3
    wl  : wavelength in m
    T   : Temperature in K

    Range of validity is 182 nm < wl < 2770 nm
    """
    # print('refractive index of water calculated from Thormählen')
    # numerical coefficients:
    a1 = 3.036167e-3
    a2 = 0.052421
    a3 = 2.117579e-1
    a4 = -5.195756e-2
    a5 = 5.922248e-2
    a6 = -1.918429e-2
    a7 = 2.582351e-3
    a8 = -2.352054e-4
    a9 = 3.964628e-5
    a10 = 3.336153e-2
    a11 = -4.008264e-2
    a12 = 8.339681e-3
    a13 = - 1.054741e-2
    a14 = 9.491575e-3

    wl = wl*1e9

    _rho = rho / 1000
    _wl = wl / 589
    _T = T/273.15
    # x is (n^2 - 1) / (n^2 + 1) / 1/_rho
    R = a1 / (_wl**2 - a2) + a3 +\
            (a4 + a5*_wl + a6 * _wl**2 + a7 * _wl**3 + a8 * _wl**4) * (_wl**2) + \
            a9/_rho + (a10 + a11*_wl + a12* _wl**2) * (_wl**2) * _T +\
            (a13 + a14*_wl) * (_wl) * (_T**2)

    n = np.sqrt((1 + 2 * R * _rho) / (1 - _rho * R))
    return n

def visc_water(TK):
    """
    Article (Kestin1978Viscosity)

    Kestin, Joseph / Sokolov, Mordechai / Wakeham, William A.
    Viscosity of liquid water in the range -8 textdegreeC to 150 textdegreeC
    1978-07

    Journal of Physical and Chemical Reference Data , Vol. 7, No. 3
    AIP Publishing
    p. 941-948
    """

    # print('using visc of water from Kestin')
    t = TK - 273.15
    mu20 = 1002.1e-6

    R = (20 - t) / (t+96) * \
        (1.2378 - 1.303e-3 * (20-t) \
        + 3.06e-6 * (20-t)**2 \
        + 2.55e-8 * (20-t)**3)
    visc = np.exp(R) * mu20
    return visc


class Polymer:

    def __init__(self, pid):
        self.type = 'polymer'
        with dbo.dbopen() as c:
            c.execute("SELECT * FROM polymers WHERE id=?", (pid,))
            values = c.fetchone()
        self.id, self.name, self.short_name, self.comment = values

class Protein:

    def __init__(self, id):
        self.type = 'protein'
        with dbo.dbopen() as c:
            c.execute("SELECT * FROM proteins WHERE id=?", (id,))
            values = c.fetchone()
        self.id, self.CAS_number, self.name, \
            self.short_name, self.EC_number, self.molecular_weight,\
            self.pI, self.reference_size, self.comment = values

class Polysaccharide:

    def __init__(self, id):
        self.type = 'polysaccharide'
        with dbo.dbopen() as c:
            c.execute("SELECT * FROM polysaccharides WHERE id=?", (id,))
            values = c.fetchone()
        self.id, self.CAS_number, self.name, self.short_name, self.EC_number,\
            self.molecular_weight, self.reference_size, self.comment = values

@dataclass
class Sample:

    id: str = 'X'
    type: str = 'generic'
    datestring: str = 'xxxx-xx-xx'
    timestring: str = 'xx:xx'
    # buffer: Sample = None

    def get_samplestring(self):
        return f'{self.type}_{self.id}'

    def get_pdt(self):
        "returns date in datetime format and boolean for if time is includet or not"
        year = int(self.datestring[0:4])
        month = int(self.datestring[4:6])
        day = int(self.datestring[6:9])

        if self.timestring!= None:
            hour, minute = [int(x) for x in self.timestring.split(':')]
            date = dt.datetime(year, month, day, hour=hour, minute=minute)
            return date, True
        else:
            # date = df.date(year, month, day)
            # return date, False
            print('Error')

    def get_pHmeasurements(self):
        query = f"""
        SELECT * FROM pH_measurements WHERE sample='{self.get_samplestring()}';
        """
        print(query)
        with dbo.dbopen() as c:
            c.execute(query)
            phlist = c.fetchall()
        return phlist

    def get_pH(self):
        phlist = self.get_pHmeasurements()
        if len(phlist)==0:
            return None
        elif len(phlist)==1:
            ph = phlist[0][4]
        else:
            print('multiple pH measurements for this sample')
            print('Using the last one')
            ph = phlist[-1][4]
        return ph

    def get_n(self, TK, wl):
        return self.buffer.get_n(TK, wl)

    def get_viscosity(self, TK):
        return self.buffer.get_viscosity(TK)


class Buffer(Sample):

    def __init__(self, bid, pH, comment):
        self.type = 'ideal_buffer'
        self.id = bid
        self.name = bid
        self.pH = pH
        self.comment = comment

    @classmethod
    def from_db(cls, bid):
        with dbo.dbopen() as c:
            c.execute("SELECT id, pH, comment FROM buffers WHERE id=?",
                      (bid,))
            values = c.fetchone()
        return cls(*values)

    def get_viscosity(self, TK):
        if self.id == 0:
            return visc_water(TK)
        elif self.name == 'ETOH':
            if TK == 293.15:
                return 1.159
        else:
            print('using the viscosity of pure H2O (In the buffer class)', visc_water(TK))
            return visc_water(TK)

    def get_n(self, TK, wl):
        if self.id == 0:
            return n_water(TK, wl)
        else:
            print('Using the refractive index of water. (In the Buffer Class)', n_water(TK, wl))
            return n_water(TK, wl)


class BufferSample(Buffer):

    def __init__(self, bid, datestring, ideal_id, comment):
        self.id = bid
        self.name = bid
        self.datestring = datestring
        self.comment = comment
        self.ideal = Buffer.from_db(ideal_id)

    @classmethod
    def from_db(cls, bid):
        with dbo.dbopen() as c:
            query = """SELECT id, datestring, ideal_id, comment
            FROM buffer_samples WHERE id=?"""
            c.execute(query, (bid,))
            values = c.fetchone()
        return cls(*values)


class mix_Sample(Sample):

    def __init__(self, id):
        with dbo.dbopen() as c:
            c.execute("SELECT * FROM mix_samples WHERE id=?", (id,))
            values = c.fetchone()
            self.id, self.datestring, self.timestring, self.comment = values

            c.execute("SELECT * FROM mix_components WHERE mix_id=?", (id,))
            components = c.fetchall()
            self.components = []
            self.volume = 0.
            for comp in components:
                _, samplestring, mycroliter = comp
                stype, sid = samplestring.split('_')
                sample = get_sample(stype, sid)
                self.componets.append([sample, mycroliter])
                self.volume += mycroliter

    def get_concentration(self, type, id):
        pass

    def get_c(self, componentstring):
        for comp in self.components:
            sample, mycroliter  = comp
            string = sample.get_samplestring()
            if samplestring == componentstring:
                return mycroliter/self.totalvolume


def get_polymer(type, id):
    if type=='protein':
        return Protein(id)
    elif type == 'polysaccharide':
        return Polysaccharide(id)
    elif type == 'polymer':
        return Polymer(id)
    else:
        raise Exception('No polymer type of this sort found')


class GenericSample(Sample):

    def __init__(self, gid, datestring, timestring, polymer, gpl, buffer_id,
                 comment):
        self.id = gid
        self.datestring = datestring
        self.gpl = gpl
        self.comment = comment
        self.timestring = timestring
        polymer_type, polymer_id = polymer.split('_')
        self.polymer = get_polymer(polymer_type, polymer_id)
        try:
            self.buffer = Buffer(buffer_id)
            print('Buffer is ', buffer_id)
        except:
            print('No buffer found')
            self.buffer = None
        self.type = 'generic'

    @classmethod
    def from_db(cls, gid):
        with dbo.dbopen() as c:
            query = """SELECT id, creation_date, ptime, polymer, gpl,
            buffer_id, comment
            FROM generic_samples
            WHERE id=?"""
            c.execute(query, (gid,))
            values = c.fetchone()
            print(values)
        return cls(*values)

    def get_molpg(self):
        gpl = self.gpl
        # molecular mass is saved in Dalton, which is g/mol (not SI)
        gpmol = self.polymer.molecular_weight
        molpl = gpl / gpmol
        # Actually we have to devide this by the gpl of water:
        molpg = molpl / 1000
        return molpg

    def get_molpl(self):
        gpl = self.gpl
        Dalton = self.polymer.molecular_weight # g/mol
        molpl = gpl/Dalton
        return molpl


typedict = {
    'pps': pps_Sample,
    # 'buffer': Buffer,
    'buffer': BufferSample,
    'trytry': TryTrySample,
    'polymer': polymer_Sample,
    'generic': GenericSample,
}

nodbdict = {
    'tresy': TresySample
    }



def get_sample(stype, sid):
    if stype in nodbdict:
        return nodbdict[stype](sid)
    if stype in typedict:
        return typedict[stype].from_db(sid)
    else:
        print('No samples of that type are implemented...')
