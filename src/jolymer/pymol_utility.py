from __future__ import print_function

import pymol
from pymol import cmd
import numpy as np
import pandas as pd
import math
import subprocess
import getpass

import matplotlib.pyplot as plt
import matplotlib.colors as mc

from . import plot_utility as plu


def run_program(name, *args):
    cmd = ['powershell', f'{name}', *args]
    return subprocess.run(cmd, cwd=pdb_directory)


# pdb2pqr wrapper stuff:
pdb_directory = f"C:\\Users\\{getpass.getuser()}\\LRZ Sync+Share\\master-thesis\\pdb_files"

def pdb2pqr_help():
    out = run_program('pdb2pqr30', '--help')
    return out

def pdb2pqr(infile, outfile, *args):
    out = run_program('pdb2pqr30', infile, outfile, *args)
    return out

cmd.set('ray_opaque_background', 0)

def rgyrate(selection='(all)', quiet=1):
    '''
DESCRIPTION

    Radius of gyration

USAGE

    rgyrate [ selection ]
    '''
    try:
        from itertools import izip
    except ImportError:
        izip = zip
    quiet = int(quiet)
    model = pymol.cmd.get_model(selection).atom
    x = [i.coord for i in model]
    mass = [i.get_mass() for i in model]
    xm = [(m*i,m*j,m*k) for (i,j,k),m in izip(x,mass)]
    tmass = sum(mass)
    rr = sum(mi*i+mj*j+mk*k for (i,j,k),(mi,mj,mk) in izip(x,xm))
    mm = sum((sum(i)/tmass)**2 for i in izip(*xm))
    rg = math.sqrt(rr/tmass - mm)
    if not quiet:
        print("Radius of gyration: %.2f" % (rg))
    return rg

# import center_of_mass
# fetch 1c3y, finish=1, multiplex=0

# com 1c3y, state=1
#Create a pseudoatom representing the 1c3y COG and store it as "1c3y_COM"
#The "1c3y_COM" object will contain 1 state only

# com 1c3y, state=1, object=COG
#Create a pseudoatom representing the 1c3y COG and store it as "COG"
#The "COG" object will contain 1 state only

# com 1c3y, state=1, object=COM, mass=1
#Create a pseudoatom representing the 1c3y COM and store it as "COM"
#The "COM" object will contain 1 state only

# com 1c3y, object=COM, mass=1
#Create a single pseudoatom containing the COM for each state found in 1c3y and store it as "COM"
#The "COM" object will contain MULTIPLE states!



def com(selection, state=None, mass=None, object=None, quiet=1, **kwargs):
    quiet = int(quiet)
    if (object == None):
        try:
            object = cmd.get_legal_name(selection)
            object = cmd.get_unused_name(object + "_COM", 0)
        except AttributeError:
            object = 'COM'
    cmd.delete(object)

    if (state != None):
        x, y, z = get_com(selection, mass=mass, quiet=quiet)
        if not quiet:
            print("%f %f %f" % (x, y, z))
        cmd.pseudoatom(object, pos=[x, y, z], **kwargs)
        cmd.show("spheres", object)
    else:
        for i in range(cmd.count_states()):
            x, y, z = get_com(selection, mass=mass, state=i + 1, quiet=quiet)
            if not quiet:
                print("State %d:%f %f %f" % (i + 1, x, y, z))
            cmd.pseudoatom(object, pos=[x, y, z], state=i + 1, **kwargs)
            cmd.show("spheres", 'last ' + object)

cmd.extend("com", com)


def get_com(selection, state=1, mass=None, quiet=1):
    """
 DESCRIPTION

    Calculates the center of mass

    Author: Sean Law
    Michigan State University
    slaw (at) msu . edu
    """
    quiet = int(quiet)

    totmass = 0.0
    if mass != None and not quiet:
        print("Calculating mass-weighted COM")

    state = int(state)
    model = cmd.get_model(selection, state)
    x, y, z = 0, 0, 0
    for a in model.atom:
        if (mass != None):
            m = a.get_mass()
            x += a.coord[0] * m
            y += a.coord[1] * m
            z += a.coord[2] * m
            totmass += m
        else:
            x += a.coord[0]
            y += a.coord[1]
            z += a.coord[2]

    if (mass != None):
        return x / totmass, y / totmass, z / totmass
    else:
        return x / len(model.atom), y / len(model.atom), z / len(model.atom)

cmd.extend("get_com", get_com)

def plot_rg():
    pymol.cmd.extend("rgyrate", rgyrate)
    pymol.cmd.set('sphere_transparency', 0.5)

    # pymol.cmd.fetch('2xwu')
    pymol.cmd.hide('everything')
    pymol.cmd.show('cartoon', 'chain A')

    r = rgyrate('chain A and polymer')
    com('chain A and polymer', object='com', vdw=r)

    pymol.util.cbc()
    print('test2')

cmd.set_color('red', mc.to_rgb(plu.tum_dred))
cmd.set_color('blue', mc.to_rgb(plu.tum_s1))
