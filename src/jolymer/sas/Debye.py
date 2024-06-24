"""
Module to calculate the scattering function corresponding to a structure file.
The calculation is based on the debye equation only.
"""

import numpy as np
import xraydb
from matplotlib import pyplot as plt
import pandas as pd

def get_coords_atom(file_name):
    coords = []
    atom_type = []
    with open(file_name) as f:
        for line in f:
            if line[0:6].strip() == 'ATOM': # need to spit by number of characters as columns join at times
                coords.append([line[30:38].strip(), line[38:46].strip(), line[46:54].strip()]) # get x, y, z coodinates
                atom_type.append(line[76:78].strip()) # read element symbol column
                #atom_type.append(line[12])          # might have to change depending on if you have element symbol i.e. 'C'
    return np.array(coords).astype('float64'), atom_type

def _debye_calc(coords, atom_type, q):
    f = []
    print('starting scattering factor calc (f), might take a while')
    for item in atom_type: f.append(xraydb.f0(item, q/4/np.pi)) # scattering vecor as sin(theta)/lambda
    print('finished scattering factor calc')

    #calculate Debye equation part porting from matlab

    X1, X2 = np.meshgrid(coords[:,0], coords[:,0])
    Y1, Y2 = np.meshgrid(coords[:,1], coords[:,1])
    Z1, Z2 = np.meshgrid(coords[:,2], coords[:,2])

    euc_dist_all = np.sqrt((X1- X2)**2 + (Y1 - Y2)**2+(Z1 - Z2)**2) # get all Euclidian distances
    #print(euc_dist_all)
    S = []
    for index, item in enumerate(q): # loop though q
        fq = []
        for f_item in f: fq.append(f_item[index]) # get f for all atoms at that q
        fq1, fq2 = np.meshgrid(fq, fq)
        part1 = np.float64(fq1*fq2*np.sin(item*euc_dist_all)) / np.float64(item*euc_dist_all) # f1*f2*sin(qr)/qr
        S.append(np.sum(np.array(fq)**2) + np.nansum(np.nansum(part1))) #below
        # first part corrects for divide by 0. sin(qr)/qr should be 1 when r=0, second is sum
        print("Current q : " + str(item))
    return S

def debye(filename, ref_df=None):
    if ref_df is None:
        q = np.linspace(start_q,end_q,num_q)
    else:
        q = np.array(ref_df.q)
    coords, atom_type = get_coords_atom(filename)
    S = _debye_calc(coords, atom_type, q)
    outdf = pd.DataFrame({'q': q, 'I': S})
    return outdf


