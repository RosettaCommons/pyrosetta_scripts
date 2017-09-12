#!/home/sheffler/venv/david/bin/ipython
from pdb_utils_noclass import *
#from __future__ import print_function
import pyrosetta
import math
import random
import collections
from rosetta import *
from rosetta.numeric import xyzVector_double_t as V3
from rosetta.numeric import xyzMatrix_double_t as M3
import rosetta.core.pack.rotamer_set
import pickle
import string
from sys import argv
from rosetta.protocols.protein_interface_design.movers import TryRotamers
from rif.geom.ray_hash import RayRay10dHash
from rif.geom import Ray
import rotamer_utilities
import hash_subclass

#from user specified input pdb file with amino acid docked on target molecule (protein backbone,small molecule, surface, etc)
#and constraint file, generate hash table by MonteCarlo starting with input pdb file using atom_pair_constraints as score function. the user
#for now should also specify the numbers of 3 atoms in the target from which the coordinate frame will be derived; the first should be
#close to the center. in the input pdb file, the target should be first, then the single docked amino acid.  the residue number and atom name
#of the contacting atoms in target and amino acid must also be input

# user supplied inputs:
# input_file
# Rosetta atom pair constraint function with constraints specified (or just supply constraints)
# origin_atom, x_atom, y_atom
# target_jump_res, target_jump_atom, aa_jump_atom
# new_hash_file     new hash file name

def get_starting_pose_from_pdb():
    p=rosetta.core.import_pose.pose_from_file(input_pdb_file)
    q=move_jump(p.clone())
    return q

# move jump to contacting atom pairs in target and aa for sampling in MC
def move_jump(p):
    chains=p.split_by_chain()
    target=chains[1]
    aa=chains[2]
    target.append_pose_by_jump(aa,target_jump_res,'target_jump_atom','aa_jump_atom')
    return target

##define score function here, or pass one in!
def score_fxn():
    sf = pyrosetta.get_score_function()
    sf.set_weight(core.scoring.fa_dun, 0.0)    
    sf.set_weight(core.scoring.fa_atr, 0.0)    
    sf.set_weight(core.scoring.fa_rep, 1.0)    
    sf.set_weight(core.scoring.fa_sol, 0.0)    
    sf.set_weight(core.scoring.fa_intra_rep, 0.0)    
    sf.set_weight(core.scoring.fa_elec, 0.0)
    sf.set_weight(core.scoring.hbond_bb_sc, 0.0)
    sf.set_weight(core.scoring.hbond_lr_bb, 0.0)    
    sf.set_weight(core.scoring.hbond_sr_bb, 0.0)    
    sf.set_weight(core.scoring.hbond_sc, 0.0)    
# if using atom_pair constraints
# sf.set_weight(core.scoring.atom_pair_constraint,1.0)
# need to add function to add constraints to pose, for example:
# I can generalize this if you give me an example constraint file format
# atomA=[331,"CA"] #ResNDX and ATOM
# atomB=[442,"CA"]
# distance_AB=5.0 #A
# dist_constraintA=pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(pyrosetta.rosetta.core.id.AtomID(pose.residue(atomA[0]).atom_index(atomA[1]),atomA[0]),
# pyrosetta.rosetta.core.id.AtomID(ff_pose.residue(atomB[0]).atom_index(atomB[1]),atomB[0]), pyrosetta.rosetta.core.scoring.func.HarmonicFunc( distance_AB, 0.5 )  #Distance, SD)
# pose.add_constraint(dist_constraintA)
    return sf   

init_pyrosetta()
gn=get_starting_pose_from_pdb()
j = gn.jump(1)  # only one jump in this system
sf =  score_fxn()

T=0.2
total_cycles=10000
equilibration_cycles=0.01*total_cycles
sample_interval=2
dump_interval=total_cycles*.2
dump_it=0

old_trans = V3(j.get_translation())
old_rot = M3(j.get_rotation())
low_rot=old_rot
low_trans=old_trans
old_E=sf(gn)
low_E=old_E

# set hash parameters and setup hash.  WILL CHANGE DEPENDING ON GEOMETRY
resl=0.5
lever=3.0
hc=make_hash_store_rot(resl,lever)

for i in range(total_cycles):
    j.gaussian_move(1, 0.25, 1.25) 
    gn.set_jump(1, j)
    new_E=sf(gn)
    delta=new_E-old_E
    accept=1
    if delta > 0:
        if math.exp(-delta/T) < random.random():
            accept=0
    if accept:
        old_E=new_E
        old_trans = V3(j.get_translation())
        old_rot = M3(j.get_rotation())
        if i>equilibration_cycles  and new_E < -0.3:
            hc.update_table(gn)

        if dump_it> dump_interval:
            gn.dump_pdb('t%s.pdb'%i)
            dump_it=0
        else:
            dump_it=dump_it+1

        if low_E > new_E:
            low_E = new_E
            low_rot=old_rot
            low_trans=old_trans
            print('low accept: %s %s'%(low_E,low_trans))
    else:
        j.set_rotation(old_rot)
        j.set_translation(old_trans)
        gn.set_jump(1,j)

print('low energy: %s'%low_E)
j.set_rotation(low_rot)
j.set_translation(low_trans)
gn.set_jump(1,j)
gn.replace_residue(4,rot0,True)
gn.dump_pdb('final.pdb')
pickle.dump(hc.dd,open(hash_table_name,"wb"))
