#!/usr/bin/python

import os
import sys


from pyrosetta import *
from rosetta import *

init("-use_truncated_termini")

if (len(sys.argv) != 3):
    print("./copy_chain_2_aminos.py source_pdb dest_pdb")

src_pose = pose_from_pdb(sys.argv[1])
dest_pose = pose_from_pdb(sys.argv[2])




src_poses = src_pose.split_by_chain()
dest_poses = dest_pose.split_by_chain()

src = src_poses[2]
dest = dest_poses[2]

if (src.size() != dest.size()):
    print("Error: chain 2 not same size in each pose")
    sys.exit(1)

for i in range(src.size()):
    i = i + 1

    res = src.residue(i)
    cln = res.clone()

    dest.replace_residue(i, cln, True)





pose = dest_poses[1]
for po in list(dest_poses)[1:]:
    pose.append_pose_by_jump(po, 1)



name = sys.argv[2]

had_gz = False
if name.endswith(".gz"):
    had_gz = True
    name = name[:-3]

if (name.endswith(".pdb")):
    name = name[:-4]

name2 = sys.argv[1]
if name2.endswith(".gz"):
    name2 = name2[:-3]
if (name.endswith(".pdb")):
    name2 = name2[:-4]

name += "_chain_2_from_" + name2

name = name + ".pdb"

if (had_gz):
    name = name + ".gz"

pose.dump_pdb(name)
