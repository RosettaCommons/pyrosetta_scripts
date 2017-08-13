#!/usr/bin/python


from pyrosetta import *
from rosetta import *
import sys


init("-beta_nov16_cart -ex1 -ex2aro")


name = sys.argv[1]

pose = pose_from_pdb(name)


scorefxn = get_fa_scorefxn()


min_mover = protocols.simple_moves.MinMover()
min_mover.cartesian(True)
min_mover.min_type("lbfgs_armijo_nonmonotone")
min_mover.tolerance(0.000001)
min_mover.score_function(scorefxn)



mm = core.kinematics.MoveMap()
mm.set_chi(True)
mm.set_bb(False)
mm.set_jump(False)

min_mover.movemap(mm)


print ("Score before: %.1f"%scorefxn(pose))

min_mover.apply(pose)

print ("Score after: %.1f"%scorefxn(pose))



name = name[:-4] + ".min.pdb"
pose.dump_pdb(name)


