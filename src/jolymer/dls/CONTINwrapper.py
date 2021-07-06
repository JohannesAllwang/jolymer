#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 14:04:51 2020

@author: johannes
"""


import io
import codecs
import collections
import locale
import numbers
import os
import sys
import time
import getpass

import subprocess
from math import pi

import numpy as np

# As CONTIN uses cgs units we use these here too

kb = 1.3806505e-16  # cm2 g s-2 K-1

path = os.path.realpath(os.path.dirname(__file__))

# http://bugs.pythoorg/issue1731717
# this workaround at http://bugs.pythoorg/issue1236
# subprocess._cleanup = lambda: None
subprocess._cleanup = lambda: None

# try:
#     p = subprocess.Popen(['powershell', 'command -v ./jolib/dls/contin' ], shell=True,
#     # p = subprocess.Popen(['powershell', 'ls ./jolib/dls' ], shell=True,
#                          bufsize=0,
#                          stdin=subprocess.PIPE,
#                          stdout=subprocess.PIPE,
#                          stderr=subprocess.PIPE)
#     continexe = p.communicate()[0].strip()
# except:
#     print('Connection with contin failed!')
#     continexe = ''
# print(continexe)

def w2f(word):
    """converts if possible to float"""
    try:
        return float(word)
    except ValueError:
        return word


def continfit(measurement, seq_number, Ngrid=256, tmin=-2, bgr=0, distribution='x', RDG=-1, timescale=1e1, **kw):
    r"""
    https://jscatter.readthedocs.io/en/latest/dls.html
    """
    # contin is called in a shell like:  contin <input.txt >output.txt
    # we mimic this by subprocess.Popen
    # create a pipe to CONTIN with the input to stdin and stdout to
    # outputpath = os.path.realpath(os.path.dirname(__file__))

    # if continexe == '':
    #     raise Exception('There is no contin executable found. ' +
    #                     'Please compile and place executable at callable path. ' +
    #                     'See documentation of DLS module.')

    typ = -1 # g2-1
    l = 632.  # in nm
    T = measurement.get_TK()
    n = measurement.sample.buffer.get_n(T, measurement.get_wl())
    if n == None:
        n = 1.333   
        print('Using n = 1.333 from water at 25 degC')
    a = measurement.phifromseq(seq_number)
    v = measurement.sample.buffer.get_viscosity(T)
    if v == None:
        v = 0.0010016
        print('Using viscosity of water at 25 cdeg')
    if bgr != 0:
        bgr = 1  # no background

    # wavevector in 1/cm
    qq = lambda n, ll, theta: 4. * pi * n / (ll * 1e-7) * np.sin(np.deg2rad(theta) / 2.)
    # hydrodynamic radius  Rh(relaxationtime) in cm, gamma is relaxation time visc in cPoise
    Rh = lambda gamma, q, T, visc: kb * T / (pi * 0.06 * visc) * q ** 2 * gamma
    # mass_weight(t)
    massweight = lambda ci, qq, Ri: ci * Ri ** 3 * qq ** 6 / (np.sin(qq * Ri) - qq * Ri * np.cos(qq * Ri)) ** 2
    # number weight
    numberweight = lambda ci, qq, Ri: ci * qq ** 6 / (np.sin(qq * Ri) - qq * Ri * np.cos(qq * Ri)) ** 2
    ##########################################################

    # DEFINES THE kind of the distribution kernel
    distr = {'m': 1,
             'D': 2, 'L': 2,
             'r': 3,
             'u': 4, 'x': 4, 'd': 4, 'T': 4}

    edist = ('x', 'T', 'd')  # here we set R21-R23 explicitly in the later header part
    if distribution[0] == 'x':  # l^R23*exp(-R21*t*l^R22) with  R21=1,R22=-1,R23=0  ==> l=1/T -> exp(-t/T)
        R21 = 1
        R22 = -1
        R23 = 0

    # write header für CONTIN als ASCII inputfile           #######################################
    last = 1  # for multi dataset evaluation this is -1 (false) except for the last one
    # we process here ONLY single files
    elements = 40  # just to have an array for it; last 2 lines are "end" and NY
    header = np.array([''] * elements, dtype='|S70')  # a single fortran line
    header[0] = 'filename'.ljust(70)  # the loaded file
    header[1] = 'LAST'.ljust(6) + str().rjust(5) + ('%15.4E' % last)  # see above
    # header[2]='GMNMX'.ljust(6)+str(1).rjust(5)+('%15.4E' % gmin)     # first point of the distribution to fit
    # header[3]='GMNMX'.ljust(6)+str(2).rjust(5)+('%15.4E' % gmax)     # last  point of the distribution to fit
    header[4] = 'IWT'.ljust(6) + str().rjust(5) + (
            '%15.4E' % 5)  # fit strategy how to determine errors -> 5 from a prefit, results in 2 fits but good errors
    # unweighted fit IWT=1 ->errors equal; IWT=4 direct input of errors not implemented
    header[5] = 'NERFIT'.ljust(6) + str().rjust(5) + (
            '%15.4E' % 0)  # number of points around a point to determine error; safety margin default 10; we use 0
    header[6] = 'NINTT'.ljust(6) + str().rjust(5) + (
            '%15.4E' % -1)  # number of equally spaced sets in tk; <0 means direct input as used here
    header[7] = 'IFORMT'.ljust(6) + str().rjust(20)  # format of time variable for direct input
    header[8] = '(1E12.5)'.ljust(26)  # 1 in a row
    header[9] = 'IFORMY'.ljust(6) + str().rjust(20)  # format of y variable for direct input correlation
    header[10] = '(1E12.5)'.ljust(26)  # 1 in a row
    if 'IGRID' in kw.keys():
        # Grid=2 is log grid ; 1 is equally spaced grid; default 2= log grid
        header[11] = 'IGRID'.ljust(6) + str().rjust(5) + ('%15.4E' % float(kw['IGRID']))
    header[12] = 'NLINF'.ljust(6) + str().rjust(5) + ('%15.4E' % bgr)  # allows a single const background , 0 no bgr
    header[13] = 'NG'.ljust(6) + str().rjust(5) + ('%15.4E' % Ngrid)  # Ngrid  points between gmin,gmax
    header[14] = 'DOUSIN'.ljust(6) + str().rjust(5) + (
            '%15.4E' % 1)  # Do User INPUT ; to use the below given values anyway this is the default
    header[15] = 'IUSER'.ljust(6) + str(10).rjust(5) + (
            '%15.4E' % distr[distribution[0]])  # selects the kernel see help above
    header[16] = 'RUSER'.ljust(6) + str(15).rjust(5) + ('%15.4E' % n)  # refractive index
    header[17] = 'RUSER'.ljust(6) + str(16).rjust(5) + ('%15.4E' % l)  # wavelength  in nm
    header[18] = 'RUSER'.ljust(6) + str(17).rjust(5) + ('%15.4E' % a)  # scattering angle in degrees
    header[19] = 'RUSER'.ljust(6) + str(18).rjust(5) + (
            '%15.4E' % T)  # absolute Temperature in K or proportionality constant
    header[20] = 'RUSER'.ljust(6) + str(19).rjust(5) + ('%15.4E' % v)  # viscosity in centipoise
    header[25] = 'RUSER'.ljust(6) + str(10).rjust(5) + ('%15.4E' % typ)  # (0) means dont change; input is g1;
    # (1) input is intensity correlation g2; =>  calculate (g2/R21-1)^0.5
    # (-1) input is g2-1,  takes only the square root
    # ALV and Zetasizer Data are g2-1 -> -1
    header[26] = 'LUSER'.ljust(6) + str(3).rjust(5) + (
            '%15.4E' % RDG)  # use rayleighDebyeGans formfactor (set to true => 1) or const 1 (set to false => -1 )
    if 'R16' in kw.keys(): header[17] = 'RUSER'.ljust(6) + str(16).rjust(5) + ('%15.4E' % float(kw['R16']))
    if 'WALL' in kw.keys(): header[24] = 'RUSER'.ljust(6) + str(24).rjust(5) + ('%15.4E' % float(
        kw['WALL']))  # wall thickness  in cm in RDG for simulating hollow spheres =0 normal sphere
    if 'ALPS1' in kw.keys(): header[27] = 'ALPST'.ljust(6) + str(1).rjust(5) + (
            '%15.4E' % float(kw['ALPS1']))  # take this alpha in preliminary error analysis
    if 'ALPS2' in kw.keys(): header[28] = 'ALPST'.ljust(6) + str(2).rjust(5) + (
            '%15.4E' % float(kw['ALPS2']))  # take this alpha in final analysis and as choosen solution
    if 'I18' in kw.keys(): header[29] = 'IUSER'.ljust(6) + str(18).rjust(5) + (
            '%15.4E' % float(kw['I18']))  # formfactor average over 2 I18+1 points
    if distribution[0] in edist or 'R21' in kw.keys(): header[22] = 'RUSER'.ljust(6) + str(21).rjust(5) + (
            '%15.4E' % R21)
    if distribution[0] in edist or 'R23' in kw.keys(): header[23] = 'RUSER'.ljust(6) + str(23).rjust(5) + (
            '%15.4E' % R23)
    if distribution[0] in edist or 'R22' in kw.keys(): header[21] = 'RUSER'.ljust(6) + str(22).rjust(5) + (
            '%15.4E' % R22)
    if 'IQUAD' in kw.keys(): header[30] = 'IQUAD'.ljust(6) + str().rjust(5) + (
            '%15.4E' % float(kw['IQUAD']))  # quadrature default=3 Simpson       1 is direct; 2 is trapezoidal
    if 'NORDER' in kw.keys(): header[30] = 'NORDER'.ljust(6) + str().rjust(5) + (
            '%15.4E' % float(kw['NORDER']))  # order regularization; default 2
    if 'PLEVEL' in kw.keys():
        word = ('%5.2f' % float(kw['PLEVEL']))
        header[31] = 'PLEVEL'.ljust(6)
        header[32] = (word * 4)  # 0.1<PLEVEL<0.9  best is not to use it, default =0.5
    if 'NONNEG' in kw.keys():
        header[33] = 'NONNEG'.ljust(6) + str().rjust(5) + (
                '%15.4E' % float(kw['NONNEG']))  # no negative values in the solution default is 1;
    else:
        #   default is 1 as a distribution has no negative values
        header[33] = 'NONNEG'.ljust(6) + str().rjust(5) + ('%15.4E' % 1)
    header[-2] = 'END'.ljust(26)
    header[-1] = 'NY'.ljust(6) + str(0).rjust(5)  # Number of datapoints is set per file later
    # ende header für CONTIN als ASCII inputfile         ##################################

    # a default grid min max is needed if it is not given explicitly
    # to transform tmin and tmax  according to this function  from kernel exp   1=t*R21*l**R22  -> l=(t*R21)**-1/R22
    if distribution[0] == 'L':
        transk = lambda t, n, l, a, T, v, R22, R21: 1 / t
    elif distribution[0] == 'm':
        transk = lambda t, n, l, a, T, v, R22, R21: (t * T * qq(n, l, a) ** 2) ** (-1 / R22)
    elif distribution[0] == 'D':
        transk = lambda t, n, l, a, T, v, R22, R21: (t * qq(n, l, a) ** 2) ** (-1)
    elif distribution[0] == 'r':
        transk = lambda t, n, l, a, T, v, R22, R21: (t * kb * T * qq(n, l, a) ** 2 / (0.06 * pi * v))
    elif distribution[0] == 'd':
        transk = lambda t, n, l, a, T, v, R22, R21: (t * kb * T * qq(n, l, a) ** 2 / (0.06 * pi * v))
    elif distribution[0] == 'u' or distribution[0] == 'x':
        transk = lambda t, n, l, a, T, v, R22, R21: (t * R21) ** (-1 / R22)
    elif distribution[0] == 'T':
        transk = lambda t, n, l, a, T, v, R22, R21: (t * R21) ** (-1 / R22)

    ###################################
    # now look at the data
    idata = 0
    for i in [1]:
        data = measurement.get_data(seq_number)
        idata += 1
        print('evaluate Nr.:', idata)
        try:
            file = measurement.name
        except AttributeError:
            file = time.strftime("%y%m%d%H%M%S", time.localtime())
        header[0] = file
        # data.extract_comm(deletechars=':[]()"')  # find extra parameters in comments

        # extract datetime from comments from ALV instrument
        try:
            timestr = list(filter(lambda ll: ll.startswith('Time') or ll.startswith('Date'), measurement.comment))[0]
            timestr = timestr.translate(None, '"').split()[2]
            data.datetime = time.mktime(time.strptime(timestr, "%d.%m.%Y%H:%M:%S"))
        except:
            data.datetime=0
        # take values from data
        # T = measurement.get_TK()
    
        # n = measurement.sample.buffer.get_n(T)
        # a = measurement.phifromseq(seq_number)
        # v = measurement.sample.buffer.get_viscosity(T)
        print('T:', T, 'v:', v)
        l = 632.
        tmin = 0
        tmax = len(data)-1
        itmin = tmin
        itmax = tmax

        # or override it
        if l != 0:
            contin_wavevector = qq(n, l, a)  # in 1/cm
        else:
            contin_wavevector = 0
        # create header for contin with parameters from datafile
        try:
            if 'qtmin' in kw:
                tmin = float(kw['qtmin']) / contin_wavevector ** 2
            if 'qtmax' in kw:
                tmax = float(kw['qtmax']) / contin_wavevector ** 2
        except AttributeError:
            raise Exception('dont use qtmin / qtmax with Laplace option. wavevector is zero')
        # searches where tmin,tmax fit into list and outputs the place
        # if tmin <= 0:
            # itmin = -tmin  # if negative purge first points
        # if tmax <= 0:
            # itmax = tmax  # same but if negative count from behind
        try:
            if 'gmin' not in kw and 'qgmin' not in kw:
                gmin = transk(data.loc[itmin, 't'], n, l, a, T, v, R22,
                              R21)  # calculate min max of the grid if not given
            elif 'qgmin' in kw:
                gmin = float(kw['qgmin']) / contin_wavevector ** 2
            else:
                gmin = float(kw['gmin'])
            if 'gmax' not in kw and 'qgmax' not in kw:
                gmax = transk(data.loc[itmax, 't'], n, l, a, T, v, R22, R21)
            elif 'qgmax' in kw:
                gmax = float(kw['qgmax']) / contin_wavevector ** 2
            else:
                gmax = float(kw['gmax'])
        except ZeroDivisionError:
            print('wavevector is zero; use qgmax/qgmin only with non laplace option')
        header[2] = 'GMNMX'.ljust(6) + str(1).rjust(5) + ('%15.4E' % min(gmin, gmax))  # fit interval min
        header[3] = 'GMNMX'.ljust(6) + str(2).rjust(5) + ('%15.4E' % max(gmin, gmax))  # fit interval max

        # header[19] = 'RUSER'.ljust(6) + str(18).rjust(5) + ('%15.4E' % T)
        header[16] = 'RUSER'.ljust(6) + str(15).rjust(5) + ('%15.4E' % n)
        header[17] = 'RUSER'.ljust(6) + str(16).rjust(5) + ('%15.4E' % l)
        header[18] = 'RUSER'.ljust(6) + str(17).rjust(5) + ('%15.4E' % a)
        header[20] = 'RUSER'.ljust(6) + str(19).rjust(5) + ('%15.4E' % v)

        contin_cwd = f"C:\\Users\\{getpass.getuser()}\\LRZ Sync+Share\\master-thesis\\packaging\\src\\jolymer\\dls"
        p = subprocess.Popen(['powershell', './contin-windows.exe'],  shell=True,
                              cwd=contin_cWd,
                              bufsize=0,
                              stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)


        # set the number of values
        lenX = len(data)
        header[-1] = 'NY'.ljust(6) + str(lenX).rjust(5)
        # now write the header and the data to a input buffer
        input = io.BytesIO()
        input.writelines([b' ' + line + b'\n' for line in header if line != b''])
        input.writelines([b' ' + (b'%8.5E' % line) + b'\n' for line in np.array(data.t)])
        input.writelines([b' ' + (b'%8.5E' % line) + b'\n' for line in np.array(data.g2)])
        # to check input
        if 'write' in kw or 'w' in kw:
            with open('./' + 'input.con', 'w') as f:
                f.writelines(input.getvalue().decode('utf-8'))

        # now run contin in a shell like environment with pipes
        # with input in stdin and collect the stout and error in output and error
        (output, error) = p.communicate(input.getvalue())
        output = output.decode('utf-8')
        error = error.decode('utf-8')
        input.close()
        if error != '':
            print(file, '      :', error)
        if len(output) == 0:
            print('there was nothing in output yet')
#             return

        # CONTIN finished and we look at the output ------------------------------------

        # sort the output to the data
        rawoutputfile = kw['rawoutputfile']
        with open(rawoutputfile, 'w') as f:
            for line in output:
                f.write(line)
        
#         outblocks = output.split('UNREGULARIZED VARIABLES')[-1].split(file)
#         if len(outblocks) < 3:
#             print('last lines of CONTIN output')
#             print(outblocks[-1])
#             raise Exception('CONTIN ended with no result; use w=1 to get output for analysis')
#         # blocks of different alpha in result (without prefit results);
#         # first repeat input ; then preliminary res; last two choosen solution
#         # second last is fit  ; last is result distribution

#         # take the fit block after Abscissa;splitlines and take only lenX;split line to words 
#         # and convert to float
#         for line in outblocks[-2].split('ABSCISSA\n')[-1].splitlines()[:lenX]:
#             print(line)
#         temp = np.r_[[[float(vv) for vv in line[:22].split()] for line in
#                       outblocks[-2].split('ABSCISSA\n')[-1].splitlines()[:lenX]]].T
#         contin_result_fit = temp[[1, 0]]  # resort to have xtime in 0 and fit_y in 1
#         contin_fits = []
#         fitqualities = []
#         peakss = []
#         momentEntireSolutions = []
#         baselines = []

#         for k in np.arange(len(outblocks))[1:-2]:  # all solutions with different alpha
#             chosen = outblocks[k].splitlines()
            
            
#             # take the chosen block line 6 to 6+Ngrid;split into words and convert to float;
#             # sometimes D is used instead of E for float 1E5
#             # order of temp -> [y, error_y, t]
#             temp = np.r_[[[float(vv) for vv in line[:31].replace('D', 'E').split()] 
#                           for line in chosen[6:Ngrid + 6]]].T
#             # fit quality in 3rd line of chosen last chosen block
#             try:
#                 RR = Rh(temp[2] * timescale, contin_wavevector, T,
#                         v)  # hydrodynamic radius temp[2] is correlation time]
#                 ci_massw = massweight(temp[0], contin_wavevector,
#                                       RR)  # mass weighted contribution ci  temp[0] is fraction
#                 ci_numw = numberweight(temp[0], contin_wavevector, RR)  # number weighted contribution
#                 # tohave [t,y,error_y] and hydrodyn Radius,  mass weight and number weighting the result output
#                 contin_fit = [np.c_[temp[[2, 0, 1]].T, RR, ci_massw, ci_numw].T]
#                 contin_fits.append(contin_fit)
#             except:
#                 contin_fit = [np.c_[temp[[2, 0, 1]].T].T]
#                 contin_fits.append(contin_fit)

#             # fitquality from contin output
#             fitquality = {}
#             name = [aa.lstrip() for aa in chosen[2].lstrip().split('    ')]
#             line = chosen[3].split()
#             if line[0] == '*':
#                 value = [float(li) for li in line[1:]]
#             else:
#                 value = [float(li) for li in line[:]]
#             for i in range(len(name)):
#                 fitquality[name[i]] = value[i]
#             fitqualities.append(fitquality)
#             # look at the peaks found at the end of solution block
                    
            
#         # now the chosen solution from CONTIN, which is the last given (repeated)
#         chosen = outblocks[-1].splitlines()
#         contin_alpha = float(chosen[3].split()[0])  # this is the choosen solution alpha
#         contin_alphalist = [f['ALPHA'] for f in fitqualities]  #
#         contin_bestFit = contin_fits[contin_alphalist.index(contin_alpha)]  # this is the choosen solution
#         data['fit'] = np.array(contin_result_fit[1]**2)
#         data['res'] = np.array(data.g2 - data.fit)

#     return data, contin_bestFit[0], contin_alpha

