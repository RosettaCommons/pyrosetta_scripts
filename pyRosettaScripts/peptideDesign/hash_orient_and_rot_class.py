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



def get_ala_asn_from_pdb():
    p=rosetta.core.import_pose.pose_from_file("asn_bb_c.pdb")
    q=move_jump(p.clone())
    return q

def move_jump(p):
    chains=p.split_by_chain()
    ala=chains[1]
    asn=chains[2]
    ala.append_pose_by_jump(asn,2,'O','ND2')
    return ala

def sc_bb_energy_fxn():
    sf = pyrosetta.get_score_function()
# do print(sf) and look for weights to see what you can modify
    sf.set_weight(core.scoring.fa_dun, 0.0)    
    sf.set_weight(core.scoring.fa_atr, 0.0)    
#    sf.set_weight(core.scoring.fa_rep, 0.0)    
    sf.set_weight(core.scoring.fa_sol, 0.2)    
    sf.set_weight(core.scoring.fa_intra_rep, 0.0)    
    sf.set_weight(core.scoring.fa_elec, 0.0)
    sf.set_weight(core.scoring.hbond_bb_sc, 2.0)
    sf.set_weight(core.scoring.hbond_lr_bb, 0.0)    
    sf.set_weight(core.scoring.hbond_sr_bb, 0.0)    
    sf.set_weight(core.scoring.hbond_sc, 0.0)    
    return sf   



init_pyrosetta()
dir=string.split(argv[1])[0]
gn=get_ala_asn_from_pdb()
print(gn)
j = gn.jump(1)  # only one jump in this system
sf =  sc_bb_energy_fxn()
old_trans = V3(j.get_translation())
old_rot = M3(j.get_rotation())
low_rot=old_rot
low_trans=old_trans
old_E=sf(gn)
low_E=old_E

T=0.2
total_cycles=10000
#total_cycles=1000000
total_cycles=int(total_cycles)
equilibration_cycles=0.01*total_cycles
#equilibration_cycles=0.1*total_cycles
sample_interval=2
data=[]
dump_interval=total_cycles*.2
dump_it=0

resl=0.5
lever=3.0
hc=make_hash_store_rot(resl,lever)

for phi in range(-150,-100,3):
 for psi in range(90,150,3):
   gn.set_phi(2,float(phi))
   gn.set_psi(2,float(psi))
   for i in range(total_cycles):
    j.gaussian_move(1, 0.25, 1.25) 
    gn.set_jump(1, j)
    new_E=sf(gn)
    delta=new_E-old_E # +0.05*(new_d - old_d)
    accept=1
    if delta > 0:
        if math.exp(-delta/T) < random.random():
            accept=0
    if accept:
#        print('accept %s %s %s %s %s '%(n_accept,new_E,old_E,low_E,delta))
        old_E=new_E
        old_trans = V3(j.get_translation())
        old_rot = M3(j.get_rotation())
        if i>equilibration_cycles  and new_E < -0.3:  #augment table here.  now write out all orientations with energy better than cutoff
#        if i>equilibration_cycles and n_accept >= sample_interval:  #augment table here
#            print(' accept:  update table ')

            n_accept=0
            hc.update_table(gn)
#            data.append(update_table(gn,rots,dd,nadded))
        else:
            n_accept=n_accept+1

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
orient_rots(gn.residue(4),rots)
rot0=rots[1]
gn.replace_residue(4,rot0,True)
gn.dump_pdb('final.pdb')
#update_table(gn,rots,dd,100)
#test_hash(s)
#outfile=open('asn_sc_bb.dat','w')

# for dat_from_one_pose in data:
#   for dat in dat_from_one_pose:
#     bb_rays=dat[0]  # hbonding backbone C O N H
#     bb_rots=dat[1]  # inverse rotamer N, CA, C
#     frame = rosetta.core.kinematics.Stub(bb_rots[0],bb_rots[1],bb_rots[2])
#     Catom=frame.global2local(bb_rays[0])
#     Oatom=frame.global2local(bb_rays[1])
#     Natom=frame.global2local(bb_rays[2])
#     Hatom=frame.global2local(bb_rays[3])
#     COray = Ray(Catom,Oatom)
#     NHray = Ray(Natom,Hatom)
#     k=binner.get_key(COray,NHray)
#     s.add(k)
#     outfile.write(('%s %s %s %s %s %s %s \n'%(bb_rays[0],bb_rays[1],bb_rays[2],bb_rays[3],bb_rots[0],bb_rots[1],bb_rots[2])))
#for k in dd.keys():
#    print k,dd[k]
pickle.dump(hc.dd,open("%s/bb_asn_dd_%s"%(dir,resl),"wb"))
#    for x in bb_rays:
#        for y in bb_rots:
#            print(x,y)
        
