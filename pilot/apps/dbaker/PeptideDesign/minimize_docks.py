#!/home/sheffler/venv/david/bin/python
from sys import argv
from pyrosetta import *
from pdb_utils_noclass import *
from pyrosetta.toolbox import pose_from_rcsb
from repeat_utils import *
import random,string
from math import *
from xyzMath import *
#from two_sided_design_pyr import *
import pickle
from rif.geom.ray_hash import RayRay10dHash
from rif.geom import Ray
from rosetta.numeric import xyzVector_double_t as V3
from rosetta.numeric import xyzMatrix_double_t as M3
from hash_subclass import *

# read in list of pdbs of peptides docked on pdbs making at least one Asn-bb hbond
# rerun hash lookup, this time keeping information on all rotamers in hash (not just first one)
# for each rotamer, generate pose with all repeat units having this rotamer at the equivalent position
# rigid body minimization, perhaps with sc_bb hbond term ramped up
# rigid body and internal coordinate minimization

def sc_bb_energy_fxn():
    sf = pyrosetta.get_score_function()
# do print(sf) and look for weights to see what you can modify
    sf.set_weight(core.scoring.fa_dun, 0.0)    
    sf.set_weight(core.scoring.fa_atr, 0.5)    
#    sf.set_weight(core.scoring.fa_rep, 0.0)    
    sf.set_weight(core.scoring.fa_sol, 0.2)    
    sf.set_weight(core.scoring.fa_intra_rep, 0.0)    
    sf.set_weight(core.scoring.fa_elec, 0.0)
    sf.set_weight(core.scoring.hbond_bb_sc, 2.0)
    sf.set_weight(core.scoring.hbond_lr_bb, 0.0)    
    sf.set_weight(core.scoring.hbond_sr_bb, 0.0)    
    sf.set_weight(core.scoring.hbond_sc, 0.0)    
    return sf

def sc_bb_only_fxn():
    sf = pyrosetta.get_score_function()
    sf = pyrosetta.create_score_function("empty")
    # do print(sf) and look for weights to see what you can modify
    sf.set_weight(core.scoring.fa_dun, 0.0)
    sf.set_weight(core.scoring.fa_atr, 0.0)
#    sf.set_weight(core.scoring.fa_rep, 0.0)
    sf.set_weight(core.scoring.fa_sol, 0.0)
    sf.set_weight(core.scoring.fa_intra_rep, 0.0)
    sf.set_weight(core.scoring.fa_elec, 0.0)
    sf.set_weight(core.scoring.hbond_bb_sc, 2.0)
    sf.set_weight(core.scoring.hbond_lr_bb, 0.0)
    sf.set_weight(core.scoring.hbond_sr_bb, 0.0)
    sf.set_weight(core.scoring.hbond_sc, 0.0)
    return sf

def set_move_map(total_residue):
    mm=pyrosetta.rosetta.core.kinematics.MoveMap()
    min_chi_vector=pyrosetta.rosetta.utility.vector1_bool(total_residue)
    min_bb_vector=pyrosetta.rosetta.utility.vector1_bool(total_residue)

    for ires in xrange(total_residue):
        min_chi_vector[ires+1]=False  #Disable SC
        min_bb_vector[ires+1]=False   #disable BB
    mm.set_chi( min_chi_vector );
    mm.set_bb( min_bb_vector );
    mm.set_jump(1, True);  #Set the jump to True!!!
    return mm

def set_move_map_free_chi(total_residue):
    mm=pyrosetta.rosetta.core.kinematics.MoveMap()
    min_chi_vector=pyrosetta.rosetta.utility.vector1_bool(total_residue)
    min_bb_vector=pyrosetta.rosetta.utility.vector1_bool(total_residue)

    for ires in xrange(total_residue):
        min_chi_vector[ires+1]=True  #Disable SC
        min_bb_vector[ires+1]=False   #disable BB
    mm.set_chi( min_chi_vector );
    mm.set_bb( min_bb_vector );
    mm.set_jump(1, True);  #Set the jump to True!!!
    return mm

def translate_repeat(seq_pos,repeat_length,n):
    return seq_pos + n*repeat_length

def get_repeat_lengths():
    rp_list=map(string.split,open('repeat+xtl.list','r').readlines())
    lengths={}
    for rp_line in rp_list:
        b=rp_line[0].index('.')
        base=rp_line[0][:b]
        lengths[base]=int(rp_line[1])
    return lengths

init_pyrosetta()
input_file =argv[1]
pdb_list=map(string.split,open(input_file,'r').readlines())
#repeat_length=42
#total_residue=168
nrepeats = 4
sf=sc_bb_energy_fxn()
sc_only=sc_bb_only_fxn()

#mm_chi=set_move_map_free_chi(total_residue)
#min_mover_chi=pyrosetta.rosetta.protocols.simple_moves.MinMover(mm_chi,sf,  "lbfgs_armijo_nonmonotone",.001,True,False,False)
print 'before hash input'
h=use_hash_rot(1.0,3.0,"bb_asn_dict_combo_1.0_rot0")
print 'after hash input'
mr = protocols.simple_moves.MutateResidue()
mr.set_res_name('ASN')
lengths=get_repeat_lengths()
for pdb_f in pdb_list:
    pdb_file=pdb_f[0]
    base=string.split(pdb_file,'/')[-1]
    DHR=string.split(base,'_')[0]
    p=rosetta.core.import_pose.pose_from_file(pdb_file)

#    print sc_only(p),sf(p),' hbond and solv+hb energy of starting complex'
    starting_complex_sc_only=sc_only(p)

    pept,prot=p.clone().split_by_chain()
    starting_pept_sc_only=sc_only(pept)
#    print sc_only(prot),sf(prot),' hbond and solv+hb energy of repeat protein'
    starting_prot_sc_only=sc_only(prot)
    print (starting_complex_sc_only - starting_pept_sc_only - starting_prot_sc_only),'starting hbond e accross interface '
    nhb,rots=h.get_all_rots(pept,prot)
    total_residue=prot.size()
    if total_residue % 4 != 0:
        print 'WARNING: ',total_residue,' not multiple of 4'
        for rp in lengths.keys():
            if DHR in rp:
                print 'assigning repeat length of ', rp,'to ', pdb_file,lengths[rp]
                repeat_length=lengths[rp]
    else:
        repeat_length=total_residue/4
    mm=set_move_map(total_residue)
    min_mover=pyrosetta.rosetta.protocols.simple_moves.MinMover(mm,sf,  "lbfgs_armijo_nonmonotone",.001,True,False,False)
    for seq_pos in rots.keys():
# identify all positions of form seq_pos+n*repeat_length
        first_pos=seq_pos%repeat_length
        equiv_pos=[first_pos + n*repeat_length for n in range(nrepeats)] 

        best_E=1000.
        best_pose=prot.clone()
        for rot in rots[seq_pos]:
            t=prot.clone()
            for pos in equiv_pos:
                mr.set_target(pos)
                mr.apply(t)
                t.set_chi(1,pos,rot[0]*10.)
                t.set_chi(2,pos,rot[1]*10.)
            protE=sc_only(t)
            t.append_pose_by_jump(pept.clone(),1)
 #           print pdb_file,len(rots[seq_pos]),rot,sc_only(t),sf(t),'before min'               
            min_mover.apply(t)
#            min_mover_chi.apply(t)
            new_E=sc_only(t)
 #           print pdb_file,len(rots[seq_pos]),rot,new_E,sf(t)
            if new_E < best_E:
                best_E=new_E
                best_pose=t.clone()
                best_protE=protE
        if best_E < starting_complex_sc_only and best_E-starting_pept_sc_only-best_protE < -15 :
            best_pose.dump_pdb('%s_rb_min'%base)
            print 'dumped file %s %s %s '%(base,starting_complex_sc_only,best_E-starting_pept_sc_only-protE)
                #now minimize
