#!/usr/bin/python

import os
import sys


from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.protocols.interface.HowDesign import *

init("-use_truncated_termini")

if (len(sys.argv) != 3):
    print("./chain_2_rmsd.py pdb1 pdb2")

src_pose = pose_from_pdb(sys.argv[1])
dest_pose = pose_from_pdb(sys.argv[2])



src_poses = src_pose.split_by_chain()
dest_poses = dest_pose.split_by_chain()

src = src_poses[2]
dest = dest_poses[2]

if (src.size() != dest.size()):
    print("Error: chain 2 not same size in each pose")
    sys.exit(1)

if (src_poses[1].size() != dest_poses[1].size()):
    print("Error: chain 1 not same size in each pose")
    sys.exit(1)


chain_sel = core.select.residue_selector.ChainSelector("2")
subset = chain_sel.apply(src_pose)

print("Unaligned rmsd: %.3f"%(subset_CA_rmsd(src_pose, dest_pose, subset, False)))
print("  Aligned rmsd: %.3f"%(subset_CA_rmsd(src_pose, dest_pose, subset, True)))




