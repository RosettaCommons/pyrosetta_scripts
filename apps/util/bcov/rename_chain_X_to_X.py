#!/usr/bin/python

import os
import sys


from pyrosetta import *
from rosetta import *

init("-use_truncated_termini -ex1 -ex2aro -extra_res_fa /Users/brian/Documents/baker/from/daniel/h3/H3i.fa.mm.params")

pose = pose_from_file(sys.argv[1])

chain1 = sys.argv[2]
chain2 = sys.argv[3]

pdb_info = pose.pdb_info()

for i in range(pose.size()):
    i = i + 1
    if (pdb_info.chain(i) == chain1):
        pdb_info.chain(i, chain2)

pose.pdb_info(pdb_info)




name = sys.argv[1]

had_gz = False
if name.endswith(".gz"):
    had_gz = True
    name = name[:-3]

if (name.endswith(".pdb")):
    name = name[:-4]

name = name + ".chain_%s_to_%s.pdb"%(chain1, chain2)

if (had_gz):
    name = name + ".gz"

pose.dump_pdb(name)
