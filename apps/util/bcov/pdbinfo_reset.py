#!/usr/bin/python

import os
import sys


from pyrosetta import *
from rosetta import *

init("-use_truncated_termini -extra_res_fa /Users/brian/Documents/baker/from/daniel/h3/H3i.fa.mm.params")

pose = pose_from_pdb(sys.argv[1])


a = core.pose.PDBInfo(pose)

pose.pdb_info(a)



name = sys.argv[1]

had_gz = False
if name.endswith(".gz"):
    had_gz = True
    name = name[:-3]

if (name.endswith(".pdb")):
    name = name[:-4]

name = name + ".pdbinfo_reset.pdb"

if (had_gz):
    name = name + ".gz"

pose.dump_pdb(name)