#!/usr/bin/python
from __future__ import division

import os
import sys


from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.protocols.interface.HowDesign import *

params_path = "/Users/brian/Documents/baker/from/daniel/h3/H3i.fa.mm.params"
if (not os.path.exists(params_path)):
    params_path = params_path.replace("/Users/brian/Documents/baker", "/home/bcov")

init("-use_truncated_termini -mute core.conformation.Conformation -mute core.import_pose.import_pose -extra_res_fa %s"%params_path)

the_pdb = sys.argv[1]
chain = sys.argv[2]


pose = pose_from_file(the_pdb)

scorefxn = get_fa_scorefxn()
scorefxn(pose)


chain_sel = core.select.residue_selector.ChainSelector(chain)

# neighborhood = rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
# neighborhood.set_distance(float(radii))
# neighborhood.set_focus_selector(chain_sel)
# neighborhood.set_include_focus_in_subset(False)


subset = chain_sel.apply(pose)
# residues = rosetta.core.select.get_residues_from_subset( 
#      )

unit = numeric.xyzVector_double_t(1, 0, 0)
unit2 = numeric.xyzVector_double_t(10000, 0, 0)
rigid_body_move(unit, 0, unit2, pose, subset)





name = sys.argv[1]

had_gz = False
if name.endswith(".gz"):
    had_gz = True
    name = name[:-3]

if (name.endswith(".pdb")):
    name = name[:-4]

name = name + ".chain_%s_away.pdb"%(chain)

if (had_gz):
    name = name + ".gz"

pose.dump_pdb(name)

