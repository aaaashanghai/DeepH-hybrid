# Script for additionally preprocess for DeepH-hybrid
# Coded by ZC Tang @ Tsinghua Univ for DeepH-hybrid support. e-mail: az_txycha@126.com
# Please follow the instructions in "README" to use this script

import numpy as np 
import h5py
import json
import argparse
import os

Bohr2Ang = 0.52918

def modify_h5(h5_in, h5_out, Rxlist, Rylist, Rzlist, max_rc, all_atoms, lat, element_rc, element_info, nao):
    dist = np.zeros(len(all_atoms))
    if h5_in == h5_out:
        print("Please don't try overwriting {} for safety concern!".fromat(h5_in))
        exit()
    with h5py.File(h5_in, 'r') as S_original_f:
        with h5py.File(h5_out, 'w', libver='latest') as S_new_f:
            S_new = {}
            for Rx in Rxlist:
                for Ry in Rylist:
                    for Rz in Rzlist:
                        mirror_atoms = all_atoms.copy()
                        mirror_atoms += (Rx * lat[0,:])[None, :] + (Ry * lat[1,:])[None, :] + (Rz * lat[2,:])[None, :]
                        for ia in range(len(all_atoms)):
                            dist[:] = np.linalg.norm((mirror_atoms-(all_atoms[ia,:])[None,:]),axis=1)
                            if len(np.where(dist<max_rc * Bohr2Ang * 2)[0]) > 0:
                                for ja in np.where(dist<max_rc * Bohr2Ang * 2)[0]:
                                    this_key = "[{}, {}, {}, {}, {}]".format(Rx,Ry,Rz,ia+1,ja+1)
                                    t_rc = (element_rc[element_info[ia]] + element_rc[element_info[ja]]) * Bohr2Ang
                                    if dist[ja] > t_rc:
                                        continue
                                    if this_key not in S_original_f.keys():
                                        S_new[this_key] = np.zeros((nao[ia+1],nao[ja+1]))
                                    else:
                                        S_new[this_key] = np.array(S_original_f[this_key],dtype=np.float32)
            for key,value in S_new.items():
                S_new_f[key] = value

def modify_DeepH_hybrid(input_path, element_rc, only_S):
    # generate Rxlist, Rylist, Rzlist according to reciprocal cell and maximal cutoff
    os.chdir(input_path)
    max_rc = max(element_rc.values())
    rlat = np.transpose(np.loadtxt("rlat.dat"))
    nRx = (int(np.ceil(max_rc*Bohr2Ang/np.pi*np.linalg.norm(rlat[0,:]))) +1) * 2 - 1
    nRy = (int(np.ceil(max_rc*Bohr2Ang/np.pi*np.linalg.norm(rlat[1,:]))) +1) * 2 - 1
    nRz = (int(np.ceil(max_rc*Bohr2Ang/np.pi*np.linalg.norm(rlat[2,:]))) +1) * 2 - 1
    Rxlist = np.arange(nRx)-int((nRx-1)/2)
    Rylist = np.arange(nRy)-int((nRy-1)/2)
    Rzlist = np.arange(nRz)-int((nRz-1)/2)

    all_atoms = np.transpose(np.loadtxt("site_positions.dat"))
    nao = {} # orbital number of every site
    element_info = np.loadtxt("element.dat")
    with open("orbital_types.dat", 'r') as ot_f:
        this_ia = 1
        line = ot_f.readline()
        while line:
            line = line.strip().split()
            this_nao = 0
            for itype in range(len(line)):
                this_nao += int(line[itype]) * 2 + 1
            nao[this_ia] = this_nao
            this_ia += 1
            line = ot_f.readline()
    lat = np.transpose(np.loadtxt("lat.dat"))

    modify_h5("overlaps.h5","overlaps_refined.h5", Rxlist, Rylist, Rzlist, max_rc, all_atoms, lat, element_rc, element_info, nao)
    os.system("mv overlaps_refined.h5 overlaps.h5")
    if not only_S:
        modify_h5("hamiltonians.h5","hamiltonians_refined.h5", Rxlist, Rylist, Rzlist, max_rc, all_atoms, lat, element_rc, element_info, nao)
        os.system("mv hamiltonians_refined.h5 hamiltonians.h5")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Modifications for DeepH-hybrid')
    parser.add_argument(
        '-i','--input_dir', type=str, default='./',
        help='path of directory to be modifyed'
        )
    parser.add_argument(
        '-c','--config', type=str, default='rc_config.json',
        help='path of the config file for modification of the cut-off radius'
        )
    parser.add_argument(
        '-S','--only_S', type=int, default=0
        )
    args = parser.parse_args()

    input_path = args.input_dir
    config_path = args.config
    only_S = bool(args.only_S)

    with open(config_path, 'r') as config_f:
        element_rc_raw = json.load(config_f)
        element_rc = {}
        for (key, value) in element_rc_raw.items():
            element_rc[int(key)] = value
    
    modify_DeepH_hybrid(input_path, element_rc, only_S)