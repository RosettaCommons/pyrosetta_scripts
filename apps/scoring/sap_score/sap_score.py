#!/usr/bin/env python

# This program accepts arguments like this:

#./sap_score.py pdb1.pdb pdb2.pdb pdb3.pdb
# or
#./sap_score.py -in:file:silent my.silent

# This program calculates total sap_score of a protein
#  It also outputs a model where the per-atom sap score is in the protein b-factors
#  So be sure to open the output in pymol to see what it is flagging

# Brian Coventry (bcov@uw.edu) RosettaCon 2020

# But more realistically 
# 
# Developability index: a rapid in silico tool for the screening of antibody aggregation propensity.


# Also, this script isn't written very well/efficiently. Working on getting this into Rosetta



import os
import sys
import math

from pyrosetta import *
from pyrosetta.rosetta import *

import numpy as np
from collections import defaultdict
import time
import argparse
import itertools
import subprocess
import time

# This came in the same folder
# or you can get it from https://github.com/bcov77/npose
import voxel_array


init("-corrections:beta_nov16  -in:file:silent_struct_type binary -keep_input_scores false -mute all"
    # " -holes:dalphaball /home/bcov/dev_rosetta/main/source/external/DAlpahBall/DAlphaBall.gcc"
    )


parser = argparse.ArgumentParser()
parser.add_argument("-in:file:silent", type=str, default="")
parser.add_argument("pdbs", type=str, nargs="*")


args = parser.parse_args(sys.argv[1:])

pdbs = args.pdbs
silent = args.__getattribute__("in:file:silent")




alpha = "ACDEFGHIKLMNPQRSTVWY"

seq = ""
for letter in alpha:
    seq += "AA%sAA"%letter

sasa_pose = pose_from_sequence(seq)
scorefxn = get_fa_scorefxn()
all_sub = core.select.residue_selector.TrueResidueSelector().apply(sasa_pose)
protocols.toolbox.pose_manipulation.repack_these_residues(all_sub, sasa_pose, scorefxn)


def get_per_atom_sasa(pose, probe_size=1.1):
    atoms = core.id.AtomID_Map_bool_t()
    atoms.resize(pose.size())
    for i in range(1, pose.size()+1):
        atoms.resize( i, pose.residue(i).natoms(), True)
    surf_vol = core.scoring.packing.get_surf_vol( pose, atoms, probe_size)
    # print(surf_vol.tot_surf)
    # print(surf_vol.surf(2, 1))  # this is per atom sasa (residue 2, atom 1)
    return surf_vol

def get_per_atom_sasa2(pose, probe_size=1.1):
    sasas = core.id.AtomID_Map_double_t()
    rsd_sasa = utility.vector1_double()
    core.scoring.calc_per_atom_sasa(pose, sasas, rsd_sasa, probe_size, False)
    return sasas


surf_vol = get_per_atom_sasa2(sasa_pose)

max_sasa = {}
for i in range(len(alpha)):
    resnum = i*5+3
    letter = alpha[i]

    sasa = 0
    res = sasa_pose.residue(resnum)
    assert(res.name1() == letter)
    for atno in range(1, res.natoms()+1):
        if ( res.atom_is_backbone(atno) ):
            continue
        sasa += surf_vol(resnum, atno)

    max_sasa[letter] = sasa

# Development of hydrophobicity parameters to analyze proteins which bear post- or cotranslational modifications
# then you subtract 0.5 from scaled
hydrophobicity = {
    "A": 0.116,
    "C": 0.18,
    "D": -0.472,
    "E": -0.457,
    "F": 0.5,
    "G": 0.001,
    "H": -0.335,
    "I": 0.443,
    "K": -0.217,
    "L": 0.443,
    "M": 0.238,
    "N": -0.264,
    "P": 0.211,
    "Q": -0.249,
    "R": -0.5,
    "S": -0.141,
    "T": -0.05,
    "V": 0.325,
    "W": 0.378,
    "Y": 0.38,
}



# script_dir = os.path.dirname(os.path.realpath(__file__))
# xml = script_dir + "/py_xml/remove_superfluous_nonpolar.xml"


# objs = protocols.rosetta_scripts.XmlObjects.create_from_file(xml)


scorefxn = core.scoring.ScoreFunctionFactory.create_score_function("beta_nov16")


def my_rstrip(string, strip):
    if (string.endswith(strip)):
        return string[:-len(strip)]
    return string


def from_vector(vec):
    xyz = np.array([0, 0, 0]).astype(float)
    xyz[0] = vec.x
    xyz[1] = vec.y
    xyz[2] = vec.z
    return xyz



the_locals = None


# Developability index: a rapid in silico tool for the screening of antibody aggregation propensity.
def sap_score(pose, name_no_suffix, out_score_map, out_string_map, suffix):


    # R from the paper
    R = 5

    surf_vol = get_per_atom_sasa2(pose)


    # get the per res base stats

    res_max_sasa = [None]
    res_hydrophobicity = [None]

    for resnum in range(1, pose.size()+1):
        letter = pose.residue(resnum).name1()
        res_max_sasa.append(max_sasa[letter])
        res_hydrophobicity.append(hydrophobicity[letter])


    # make the things required to find 5A neighbors 

    idx_to_atom = []
    xyzs = []
    atom_sasa = []

    for resnum in range(1, pose.size()+1):
        res = pose.residue(resnum)
        for at in range(1, res.natoms()+1):
            if ( res.atom_is_backbone(at) ):
                continue
            xyzs.append(from_vector(res.xyz(at)))
            idx_to_atom.append([resnum, at])
            atom_sasa.append(surf_vol(resnum, at))
            # if ( atom_sasa[-1] > 1000 ):
            #     print("AARF", atom_sasa[-1], resnum, at, res.natoms(), len(atom_sasa), surf_vol.vol(resnum, at))
            # print("AAAA", atom_sasa[-1], resnum, at, res.natoms(), len(atom_sasa), surf_vol.vol(resnum, at))

    # print("TOTAL: ", surf_vol.tot_surf)

    # return None
    atom_sasa = np.array(atom_sasa)
    idx_to_atom = np.array(idx_to_atom)
    xyzs = np.array(xyzs)

    resl = 1

    low = np.min(xyzs, axis=0) - R*2 - resl*2
    high = np.max(xyzs, axis=0) + R*2 + resl*2

    print("Making neighbor grid")
    clashgrid = voxel_array.VoxelArray(low, high, np.array([resl]*3), object)
    for idx, _ in enumerate(clashgrid.arr.flat):
        clashgrid.arr.flat[idx] = []


    for ixyz, xyz in enumerate(xyzs):
        indices = clashgrid.indices_within_x_of(R+resl, xyz)
        for index in indices:
            # print(clashgrid.arr[tuple(index)])
            clashgrid.arr[tuple(index)].append(ixyz)


    atom_grid_indices = clashgrid.floats_to_indices(xyzs)

    sap_scores = []

    pdb_info = pose.pdb_info()


    for iatom in range(len(xyzs)):
        xyz = xyzs[iatom]
        resnum, at = idx_to_atom[iatom]
        grid_index = atom_grid_indices[iatom]

        grid_list = np.array(list(clashgrid.arr[tuple(grid_index)]))

        distances = np.linalg.norm( xyzs[grid_list] - xyz, axis=-1)

        idx_within_R = grid_list[distances <= R]

        atoms_within_R = idx_to_atom[idx_within_R]
        resnums = np.unique(atoms_within_R[:,0])

        atom_score = 0
        for ot_resnum in resnums:
            ats_idx = idx_within_R[atoms_within_R[:,0] == ot_resnum]
            res_sasa = np.sum(atom_sasa[ats_idx])

            res_score = res_sasa / res_max_sasa[ot_resnum] * res_hydrophobicity[ot_resnum]
            if ( res_score > 1000 ):
                print(ot_resnum, pose.residue(ot_resnum).name1(), res_sasa, res_max_sasa[ot_resnum], res_hydrophobicity[ot_resnum])

            atom_score += res_score

        pdb_info.bfactor(resnum, at, atom_score)
        sap_scores.append(atom_score)

    sap_scores = np.array(sap_scores)

    sap_score = np.sum( sap_scores[sap_scores > 0])
    print("sap score: %.1f"%sap_score)



    out_score_map['sap_score'] = sap_score




    return pose



############### BEGIN MAIN FUNCTION ###########################

if ( silent != "" ):
    sfd_in = rosetta.core.io.silent.SilentFileData(rosetta.core.io.silent.SilentFileOptions())
    sfd_in.read_file(silent)

    pdbs = list(sfd_in.tags())

    sfd_out = core.io.silent.SilentFileData( "out.silent", False, False, "binary", core.io.silent.SilentFileOptions())


num = -1
for pdb in pdbs:
    t0 = time.time()
    print("Attempting pose: " + pdb)

    # try:
    for k in [1]:
        if ( silent == "" ):
            pose = pose_from_file(pdb)
        else:
            pose = Pose()
            sfd_in.get_structure(pdb).fill_pose(pose)

        name_no_suffix = my_rstrip(my_rstrip(os.path.basename(pdb), ".gz"), ".pdb")

        sfd = core.io.raw_data.ScoreFileData("score.sc")

        score_map = std.map_std_string_double()
        string_map = std.map_std_string_std_string()


        out_pose = sap_score(pose, name_no_suffix, score_map, string_map, "")


        core.io.raw_data.ScoreMap.add_arbitrary_score_data_from_pose( pose, score_map)
        core.io.raw_data.ScoreMap.add_arbitrary_string_data_from_pose( pose, string_map)

        sfd.write_pose( pose, score_map, name_no_suffix, string_map)
        if (out_pose != None):


            # pdb_info = core.pose.PDBInfo(pose)
            # pose.pdb_info(pdb_info)
            if ( silent == "" ):
                out_pose.dump_pdb(name_no_suffix + ".pdb")
            else:
                struct = sfd_out.create_SilentStructOP()
                struct.fill_struct(out_pose, name_no_suffix)
                sfd_out.add_structure(struct)


        seconds = int(time.time() - t0)

        print("protocols.jd2.JobDistributor: " + name_no_suffix + " reported success in %i seconds"%seconds)

    # except Exception as e:
    #     print("Error!!!")
    #     print(e)



if ( silent != "" ):
    sfd_out.write_all("out.silent", False)




















