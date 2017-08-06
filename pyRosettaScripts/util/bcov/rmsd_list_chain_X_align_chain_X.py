#!/usr/bin/python
from __future__ import division

import os
import sys


from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.protocols.interface.HowDesign import *

init("-use_truncated_termini -mute core.conformation.Conformation -mute core.import_pose.import_pose")

if (len(sys.argv) < 4):
    print("./rmsd_list_chain_X_align_chain_X.py list.list pdb.pdb A B")
    sys.exit(1)

the_list = sys.argv[1]
the_pdb = sys.argv[2]
rmsd_chain = sys.argv[3]
align_chain = sys.argv[4]


native = pose_from_file(the_pdb)

rmsd_sel = core.select.residue_selector.ChainSelector(rmsd_chain)
align_sel = core.select.residue_selector.ChainSelector(align_chain)


rmsd_subset = rmsd_sel.apply(native)
align_subset = align_sel.apply(native)

align_seq_poss = core.select.get_residues_from_subset( align_subset );

f = open(the_list)
lines = f.readlines()
f.close()

for line in lines:
    line = line.strip()
    if (len(line) == 0):
        continue
    pose = pose_from_file(line)
    aligned_rmsd = protocols.toolbox.pose_manipulation.superimpose_pose_on_subset_CA(
        pose, native, align_seq_poss, 0)
    if (aligned_rmsd > 0.5):
        print("Error: file: %s has aligned rmsd %.3f"%(line, aligned_rmsd))
        sys.exit(1)

    rmsd = subset_CA_rmsd( pose, native, rmsd_subset, False )

    print "%s %.1f"%(line, rmsd)





