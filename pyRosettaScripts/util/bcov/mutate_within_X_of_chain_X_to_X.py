#!/usr/bin/python
from __future__ import division

import os
import sys


from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.protocols.interface.HowDesign import *

init("-use_truncated_termini -ex1 -ex2aro -extra_res_fa /Users/brian/Documents/baker/from/daniel/h3/H3i.fa.mm.params")

if (len(sys.argv) < 4):
    print("./mutate_within_X_of_chain_X_to_X.py pdb.pdb radius chain ALA")
    sys.exit(1)

the_pdb = sys.argv[1]
radius = float(sys.argv[2])
the_chain = sys.argv[3]
the_aa = sys.argv[4]


pose = pose_from_file(the_pdb)
scorefxn = get_fa_scorefxn()
scorefxn(pose)



protocol = '''<ROSETTASCRIPTS><RESIDUE_SELECTORS>
    <Neighborhood  name="rs" distance="%.3f" include_focus_in_subset="false" >
            <Chain chains="%s"/>
        </Neighborhood>
    </RESIDUE_SELECTORS>
    <MOVERS>
        <MakePolyX name="do_it" residue_selector="rs" aa="%s" keep_pro="0" keep_gly="0"/>
    </MOVERS>
    <PROTOCOLS>
        <Add mover_name="do_it"/>
    </PROTOCOLS>
    </ROSETTASCRIPTS>
'''%(radius, the_chain, the_aa)


protocol_tag = utility.tag.Tag.create(protocol)
parser = protocols.rosetta_scripts.RosettaScriptsParser()

mover = parser.generate_mover_for_protocol(pose, False, protocol_tag)

mover.apply(pose)



name = the_pdb

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

name += "within_%.1f_of_%s_to_%s"%(radius, the_chain, the_aa) + name2

name = name + ".pdb"

if (had_gz):
    name = name + ".gz"

pose.dump_pdb(name)








