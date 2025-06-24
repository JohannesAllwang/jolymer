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
            # need to spit by number of characters as columns join at times
            if line[0:6].strip() == 'ATOM':
                # get x, y, z coodinates
                coords.append([line[30:38].strip(),
                               line[38:46].strip(),
                               line[46:54].strip()])
                # read element symbol column
                atom_type.append(line[76:78].strip())
                #atom_type.append(line[12])
                # might have to change depending on if you have element symbol i.e. 'C'
    return np.array(coords).astype('float64'), atom_type

def debye_calc(coords, atom_type, q):
    f = []
    print('starting scattering factor calc (f), might take a while')
    # scattering vecor as sin(theta)/lambda
    for item in atom_type: f.append(xraydb.f0(item, q/4/np.pi))
    print('finished scattering factor calc')

    #calculate Debye equation part porting from matlab
    X1, X2 = np.meshgrid(coords[:,0], coords[:,0])
    Y1, Y2 = np.meshgrid(coords[:,1], coords[:,1])
    Z1, Z2 = np.meshgrid(coords[:,2], coords[:,2])

    # get all Euclidian distances
    euc_dist_all = np.sqrt((X1- X2)**2 + (Y1 - Y2)**2+(Z1 - Z2)**2)
    #print(euc_dist_all)
    S = []
    for index, item in enumerate(q): # loop though q
        fq = []
        for f_item in f: fq.append(f_item[index]) # get f for all atoms at that q
        fq1, fq2 = np.meshgrid(fq, fq)
        part1 = np.float64(fq1*fq2*np.sin(item*euc_dist_all)) / np.float64(item*euc_dist_all) # f1*f2*sin(qr)/qr
        S.append(np.sum(np.array(fq)**2) + np.nansum(np.nansum(part1))) #below
        # first part corrects for divide by 0. sin(qr)/qr should be 1 when r=0, second is sum
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

def make_pdb_dict():
    new_dict = {
        "entry": [],
        "atom_num": [],
        "atom_name": [],
        "res": [],
        "chain": [],
        "res_num": [],
        "X": [],
        "Y": [],
        "Z": [],
        "Occ": [],
        "TF": [],
       # "SI": [],
        "Elem": [],
       # "Charge": []
        }
    return new_dict

def read_pdb(file_name):

    all_info = make_pdb_dict()

    with open(file_name) as f:
        for line in f:
            if line[0:6].strip()=='ATOM':
                all_info["entry"].append(line[0:4].strip())
                all_info["atom_num"].append(int(line[6:11].strip()))
                all_info["atom_name"].append(line[12:16].strip())
                all_info["res"].append(line[17:20].strip())
                all_info["chain"].append(line[21].strip())
                all_info["res_num"].append(int(line[22:26].strip()))
                all_info["X"].append(float(line[30:38].strip()))
                all_info["Y"].append(float(line[38:46].strip()))
                all_info["Z"].append(float(line[46:54].strip()))
                all_info["Occ"].append(float(line[54:60].strip()))
                all_info["TF"].append(float(line[60:66].strip()))
                #all_info["SI"].append(line[72:76].strip())
                all_info["Elem"].append(line[76:78].strip())
                #all_info["Charge"] = float(line[78:80].strip())

    return all_info


def sep_sol(res_name, input_dict ):
    out_dic = make_pdb_dict()

    for ele, item in enumerate(input_dict["res"]):
        if item == res_name:
            for key in out_dic:
                out_dic[key].append(input_dict[key][ele])

    return out_dic


def sep_sol_cl(res_name1, resname2, input_dict ):
    out_dic = make_pdb_dict()

    for ele, item in enumerate(input_dict["res"]):
        #print(ele, item)
        if item == res_name1 or item == resname2:
            #print(input_dict["res"][ele])
            for key in out_dic:
                out_dic[key].append(input_dict[key][ele])

    return out_dic


def sep_prot(res_name, input_dict ):
    out_dic = make_pdb_dict()

    for ele, item in enumerate(input_dict["res"]):
        if item != res_name:
            for key in out_dic: out_dic[key].append(input_dict[key][ele])
        elif item == res_name: break

    return out_dic


def euc_dist(xp,yp,zp,xs,ys,zs):
    dist = np.sqrt((xp-xs)**2 + (yp-ys)**2 + (zp-zs)**2)
    return dist


def calc_dist(prot, sol, cut_off):
    solv_shell = make_pdb_dict()

    for k in range(len(prot["X"])):
        for key in solv_shell:
            solv_shell[key].append(prot[key][k])

    for s in range(len(sol["X"])):
        if s % 10000 == 0:
            print(s)

        for p in range(len(prot["X"])):
            dist = euc_dist(prot["X"][p],prot["Y"][p],prot["Z"][p],
                    sol["X"][s],sol["Y"][s],sol["Z"][s])
            if dist < cut_off:
                for key in solv_shell:
                    solv_shell[key].append(sol[key][s])
                break
    return solv_shell


def write_pdb(all_info,fname):

    with open(fname,'w') as f:
        f.close

    f = open(fname,"a")
    for ind in range(len(all_info["entry"])):
        j = []
        j.append(all_info["entry"][ind].ljust(6))#atom#6s
        j.append(str(all_info["atom_num"][ind]).rjust(5))#aomnum#5d
        j.append("  ")
        j.append(all_info["atom_name"][ind].ljust(4))#atomname$#4s
        j.append(all_info["res"][ind].ljust(3))#resname#1s
        j.append(all_info["chain"][ind].rjust(2)) #Astring
        j.append(str(all_info["res_num"][ind]).rjust(4)) #resnum
        j.append(str('%8.3f' % (all_info["X"][ind])).rjust(12)) #x
        j.append(str('%8.3f' % (all_info["Y"][ind])).rjust(8))#y
        j.append(str('%8.3f' % (all_info["Z"][ind])).rjust(8)) #z\
        j.append(str('%6.2f'%(all_info["Occ"][ind])).rjust(6))#occ
        j.append(str('%6.2f'%(all_info["TF"][ind])).ljust(6))#temp
        j.append(all_info["Elem"][ind].rjust(12))#elname

        write_string = ''
        for val in j:
            write_string = write_string+val
        write_string = write_string + "\n"
        f.write(write_string)

    f.close()
    return

def solv_layer(file_name, outname='test', dist=3.0):
    all_info = read_pdb(file_name)

    #sol = sep_sol("SOL", all_info)
    sol_cl = sep_sol_cl("SOL", "CL", all_info)
    prot = sep_prot("SOL", all_info)
    solv_shell = calc_dist(prot, sol_cl, dist) # get all sol atoms less than 3 Angstrom
    write_pdb(sol_cl, f"{outname}-solv-shell.pdb")
    write_pdb(prot, f"{outname}-prot.pdb")
    write_pdb(solv_shell, f"{outname}-solv-shell-{dist}A.pdb")

    return
