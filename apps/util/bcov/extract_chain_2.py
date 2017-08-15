#!/usr/bin/python

import os
import sys


from pyrosetta import *
from rosetta import *

init("-use_truncated_termini")

pose = pose_from_pdb(sys.argv[1])




poses = pose.split_by_chain()

pose = poses[2]


name = sys.argv[1]

had_gz = False
if name.endswith(".gz"):
    had_gz = True
    name = name[:-3]

if (name.endswith(".pdb")):
    name = name[:-4]

name = name + ".chain_2.pdb"

if (had_gz):
    name = name + ".gz"

print("Saved as: %s"%name)

pose.dump_pdb(name)
