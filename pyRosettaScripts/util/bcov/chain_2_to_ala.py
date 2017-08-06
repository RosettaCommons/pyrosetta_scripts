#!/usr/bin/python

import os
import sys


from pyrosetta import *
from rosetta import *

init("-use_truncated_termini -ex1 -ex2aro -extra_res_fa /Users/brian/Documents/baker/from/daniel/h3/H3i.fa.mm.params")

pose = pose_from_file(sys.argv[1])

amino_to_use = "ALA"
if (len(sys.argv) >= 3 and sys.argv[2] != ""):
    amino_to_use = sys.argv[2]

number = 2
if (len(sys.argv) >= 4 and sys.argv[3] != ""):
    number = int(sys.argv[3])



poses = pose.split_by_chain()


poly = protocols.simple_moves.MakePolyXMover(amino_to_use, False, False, True)

poly.apply(poses[number])


pose = poses[1]
for po in list(poses)[1:]:
    pose.append_pose_by_jump(po, 1)


chain2 = core.select.residue_selector.ChainSelector(str(number))
subset = chain2.apply(pose)









tf = rosetta.core.pack.task.TaskFactory()

tf.push_back(rosetta.core.pack.task.operation.InitializeFromCommandline())
tf.push_back(rosetta.core.pack.task.operation.IncludeCurrent())
tf.push_back(rosetta.core.pack.task.operation.NoRepackDisulfides())
tf.push_back(rosetta.core.pack.task.operation.RestrictToRepacking())

# don't repack the rest of the protein
tf.push_back(rosetta.core.pack.task.operation.OperateOnResidueSubset(
    rosetta.core.pack.task.operation.PreventRepackingRLT(), chain2, True ))
                                                                     #        ^^^^
                                                                     #invert selection


# Have to convert the task factory into 
task = tf.create_task_and_apply_taskoperations( pose )


pack_mover = rosetta.protocols.simple_moves.PackRotamersMover(get_fa_scorefxn(), task)
pack_mover.apply(pose)









name = sys.argv[1]

had_gz = False
if name.endswith(".gz"):
    had_gz = True
    name = name[:-3]

if (name.endswith(".pdb")):
    name = name[:-4]

name = name + ".%s.pdb"%amino_to_use.lower()

if (had_gz):
    name = name + ".gz"

pose.dump_pdb(name)
