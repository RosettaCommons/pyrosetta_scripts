#!/home/sheffler/venv/david/bin/python
from xyzMath import *
import pyrosetta
import random

pyrosetta.init()
p=pyrosetta.rosetta.core.import_pose.pose_from_file("origin_complex.pdb")
#print p.residue(4).atom_index('CG')
#print p.residue(4).atom_index('CE2')
print 'OH=' +str(p.residue(4).atom_index('OH'))
print 'CZ=' +str(p.residue(4).atom_index('CZ'))
#print '2HB=' +str(p.residue(2).atom_index('2HB'))
#print '1HG=' +str(p.residue(2).atom_index('1HG'))
#print '2HG=' +str(p.residue(2).atom_index('2HG'))
#print '1HD=' +str(p.residue(2).atom_index('1HD'))
#print '2HD=' +str(p.residue(2).atom_index('2HD'))
#print 'CD=' +str(p.residue(2).atom_index('CD'))
#print 'CD=' +str(p.residue(2).atom_index('CD'))
#print p.residue(2).atom_index('1HD')
#center = p.xyz(pyrosetta.rosetta.core.id.AtomID(7,4)) + p.xyz(pyrosetta.rosetta.core.id.AtomID(10,4)) + p.xyz(pyrosetta.rosetta.core.id.AtomID(11,4))
#print center[0]
#print x
#print y
#print z
