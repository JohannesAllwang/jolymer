from __future__ import print_function

# import pymol
# from pymol import cmd
import numpy as np
import pandas as pd
import math
import subprocess
import getpass
from os.path import join

import matplotlib.pyplot as plt
import matplotlib.colors as mc

from . import plot_utility as plu
from . import os_utility as osu


pdbdir = r'/home/johannes/LRZ Sync+Share/master-thesis/pdb_files/'
pdbdir = None

def run_pml(name, cwd='.'):
    cmd = ['pymol', f'{name}']
    return subprocess.run(cmd, cwd=cwd)


# pdb2pqr wrapper stuff:
pdb_directory = f"C:\\Users\\{getpass.getuser()}\\LRZ Sync+Share\\master-thesis\\pdb_files"


def pdb2pqr(name, pH=7, cwd='.'):
    cmd = ['pdb2pqr30',
           '--titration-state-method=propka',
           f'--with-ph={pH}',
           '--ff=PARSE',
           f'--apbs-input={name}_pH{pH}.in',
           '--drop-water',
           f'{name}.pdb',
           f'{name}_pH{pH}.pqr']
    return subprocess.run(cmd, cwd=cwd)


def apbs(name, pH=7, cwd='.'):
    cmd = ['apbs', f'{name}_pH{pH}.in']
    return subprocess.run(cmd, cwd=cwd)


def surfacecharge(name, pH=7, cwd='.'):
    pdb2pqr(name, pH=pH, cwd=cwd)
    apbs(name, pH=pH, cwd=cwd)
    pml_file = f'{name}_pH{pH}.pml'
    with open(join(cwd, pml_file), 'w') as f:
        f.write('# Drag this script into an open PyMOL window\n')

        f.write('# Load the files\n')
        f.write('from jolymer.pymol_utility import *\n')
        f.write("cmd.set('ray_opaque_background', 0)\n")
        f.write(f'load {name}_pH{pH}.pqr, molecule\n')
        f.write(f'load {name}_pH{pH}.pqr-PE0.dx, electrostaticmap\n')

        f.write('# Set scale for coloring protein surface\n')
        f.write('ramp_new espramp, electrostaticmap, [ -3, 0, 3]\n')

        f.write('# Show the surface\n')
        f.write('show surface\n')

        f.write('# Set surface colors from dx\n')
        f.write('set surface_color, espramp\n')
        f.write('set surface_ramp_above_mode\n')
        f.write('ray')

        # f.write('# Setup export\n')
        # f.write('set pse_export_version, 1.7\n')

        # f.write('# Save file as .pse')
        # f.write('save hz90dpvisq_APBS.pse')
    # run_pml(pml_file)


def rgyrate(selection='(all)', quiet=1):
    '''
    DESCRIPTION

    Radius of gyration

    USAGE

    rgyrate [ selection ]
    '''
    import pymol
    try:
        from itertools import izip
    except ImportError:
        izip = zip
    # cmd.extend("com", com)
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
    if mass is not None and not quiet:
        print("Calculating mass-weighted COM")

    state = int(state)
    model = cmd.get_model(selection, state)
    x, y, z = 0, 0, 0
    for a in model.atom:
        if (mass is not None):
            m = a.get_mass()
            x += a.coord[0] * m
            y += a.coord[1] * m
            z += a.coord[2] * m
            totmass += m
        else:
            x += a.coord[0]
            y += a.coord[1]
            z += a.coord[2]

    if (mass is not None):
        return x / totmass, y / totmass, z / totmass
    else:
        return x / len(model.atom), y / len(model.atom), z / len(model.atom)


def plot_rg():
    cmd.extend("get_com", get_com)
    pymol.cmd.extend("rgyrate", rgyrate)
    pymol.cmd.set('sphere_transparency', 0.5)

    # pymol.cmd.fetch('2xwu')
    pymol.cmd.hide('everything')
    pymol.cmd.show('cartoon', 'chain A')

    r = rgyrate('chain A and polymer')
    com('chain A and polymer', object='com', vdw=r)

    pymol.util.cbc()
    print(r, 'nm')

from pymol.cgo import *
from math import *
from pymol import cmd
from re import *

# def spectrumbar (*args, **kwargs):

#     """
#     Author Sean M. Law
#     University of Michigan
#     seanlaw_(at)_umich_dot_edu

#     USAGE

#     While in PyMOL

#     run spectrumbar.py

#     spectrumbar (RGB_Colors,radius=1.0,name=spectrumbar,head=(0.0,0.0,0.0),tail=(10.0,0.0,0.0),length=10.0, ends=square)

#     Parameter     Preset         Type     Description
#     RGB_Colors    [1.0,1.0,1.0]  N/A      RGB colors can be specified as a
#                                           triplet RGB value or as PyMOL
#                                           internal color name (i.e. red)
#     radius        1.0            float    Radius of cylindrical spectrum bar
#     name          spectrumbar    string   CGO object name for spectrum bar
#     head          (0.0,0.0,0.0)  float    Starting coordinate for spectrum bar
#     tail          (10.0,0.0,0.0) float    Ending coordinate for spectrum bar
#     length        10.0           float    Length of spectrum bar
#     ends          square         string   For rounded ends use ends=rounded

#     Examples:

#     spectrumbar red, green, blue
#     spectrumbar 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0

#     The above two examples produce the same spectrumbar!

#     spectrumbar radius=5.0
#     spectrumbar length=20.0

#     """

#     rgb=[1.0, 1.0, 1.0]
#     name="spectrumbar"
#     radius=1.0
#     ends="square"
#     x1=0
#     y1=0
#     z1=0
#     x2=10
#     y2=0
#     z2=0
#     num=re.compile('[0-9]')
#     abc=re.compile('[a-z]')

#     for key in kwargs:
#         if (key == "radius"):
#             radius = float(kwargs["radius"])
#         elif (key == "name"):
#             name=kwargs["name"]
#         elif (key == "head"):
#             head=kwargs["head"]
#             head=head.strip('" []()')
#             x1,y1,z1=map(float,head.split(','))
#         elif (key == "tail"):
#             tail=kwargs["tail"]
#             tail=tail.strip('" []()')
#             x2,y2,z2=map(float,tail.split(','))
#         elif (key == "length"):
#             if (abc.search(kwargs["length"])):
#                 print "Error: The length must be a value"
#                 return
#             else:
#                 x2=float(kwargs["length"]);
#         elif (key == "ends"):
#             ends=kwargs["ends"]
#         elif (key != "_self"):
#             print "Ignoring unknown option \""+key+"\""
#         else:
#             continue

#     args=list(args)
#     if (len(args)>=1):
#         rgb=[]
#     while (len(args)>=1):
#         if (num.search(args[0]) and abc.search(args[0])):
#             if (str(cmd.get_color_tuple(args[0])) != "None"):
#                 rgb.extend(cmd.get_color_tuple(args.pop(0)))
#             else:
#                 return
#         elif (num.search(args[0])):
#             rgb.extend([float(args.pop(0))])
#         elif (abc.search(args[0])):
#             if (str(cmd.get_color_tuple(args[0])) != "None"):
#                 rgb.extend(cmd.get_color_tuple(args.pop(0)))
#             else:
#                 return
#         else:
#             print "Error: Unrecognized color format \""+args[0]+"\""
#             return

#     if (len(rgb) % 3):
#         print "Error: Missing RGB value"
#         print "Please double check RGB values"
#         return

#     dx=x2-x1
#     dy=y2-y1
#     dz=z2-z1
#     if (len(rgb) == 3):
#         rgb.extend([rgb[0]])
#         rgb.extend([rgb[1]])
#         rgb.extend([rgb[2]])
#     t=1.0/(len(rgb)/3.0-1)
#     c=len(rgb)/3-1
#     s=0
#     bar=[]

#     while (s < c):
#         if (len(rgb) >0):
#             r=rgb.pop(0)
#             g=rgb.pop(0)
#             b=rgb.pop(0)
#         if (s == 0 and ends == "rounded"):
#             bar.extend([COLOR, float(r), float(g), float(b), SPHERE, x1+(s*t)*dx, y1+(s*t)*dy, z1+(s*t)*dz, radius])
#         bar.extend([CYLINDER])
#         bar.extend([x1+(s*t)*dx, y1+(s*t)*dy, z1+(s*t)*dz])
#         bar.extend([x1+(s+1)*t*dx, y1+(s+1)*t*dy, z1+(s+1)*t*dz])
#         bar.extend([radius, float(r), float(g), float(b)])
#         if (len(rgb) >= 3):
#             bar.extend([float(rgb[0]), float(rgb[1]), float(rgb[2])])
#             r=rgb[0]
#             g=rgb[1]
#             b=rgb[2]
#         else:
#             bar.extend([float(r), float(g), float(b)])
#         if (s == c-1 and ends == "rounded"):
#             bar.extend([COLOR, float(r), float(g), float(b), SPHERE, x1+(s+1)*t*dx, y1+(s+1)*t*dy, z1+(s+1)*t*dz, radius])
#         s=s+1

#     cmd.delete(name)
#     cmd.load_cgo(bar,name)


#     return
# cmd.extend("spectrumbar",spectrumbar)
