#!/usr/bin/python

from pyrosetta import *
from rosetta import *


length = int(sys.argv[1])

init()

pose = pose_from_sequence("A"*length)

for i in range(1, length+1):
    pose.set_phi(i, -57)
    pose.set_psi(i, -47)
    pose.set_omega(i, 180)


pose.dump_pdb("ideal_helix_%i.pdb"%length)




