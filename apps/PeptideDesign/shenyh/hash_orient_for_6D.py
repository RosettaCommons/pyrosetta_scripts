#!/software/miniconda3/envs/pyrosetta3/bin/python
import pyrosetta
import argparse
import math
import random
import collections
from rosetta import *
from rosetta.numeric import xyzVector_double_t as V3
from rosetta.numeric import xyzMatrix_double_t as M3
import pickle
from rosetta.protocols.protein_interface_design.movers import TryRotamers
from rif.geom.ray_hash import RayRay10dHash
from rif.geom import Ray
from hash_subclass_for_6D import *


def get_ala_asn_from_pdb(p):
    q=move_jump(p.clone())
    return q
#set up the move_jump
def move_jump(p):
    chains=p.split_by_chain()
    Surf=chains[1]
    Resi=chains[2]
    #change here
    Surf.append_pose_by_jump(Resi,89,'O3','CD')
    return Surf
#set up score function
def sc_energy_fxn2():
    sf = pyrosetta.get_score_function()
    #print sf
    sf.set_weight(core.scoring.fa_sol, 0.1)
    sf.set_weight(core.scoring.fa_rep, 0.55)
    sf.set_weight(core.scoring.fa_elec, 1.0)
    sf.set_weight(core.scoring.atom_pair_constraint, 1.0)
    sf.set_weight(core.scoring.angle_constraint, 1.0)
    sf.set_weight(core.scoring.dihedral_constraint, 1.0)
    return sf

#set up the parameters
options={}
options['input_pdb']='A011_2by5_CleaveSurface__GLU_highA.pdb'
options['energy_function']='-empty'
options['constraint_file']='A011_Glu_highA.cst'
options['extra_res_fa']='CAL.params CO3.params'
options['jump']=120
options['Gaus_trans_mag']=float(0.15)
options['Gaus_rot_mag']=float(1.5)
options['Energy_cut']=-3235
options['Resi_num']=121
options['Resi_name']='GLU'
options['hashtable_name']='CaCO3_GLU_hashtable'
options['target_phi_psi_list'] = [(-118.5, 118.3)]
options['orient_atoms'] = ['OE1','CD','OE2']
options['create_ray_atoms_num']=[83,89,89]
options['create_ray_atoms_name']=['Ca2p','C1','O1']
options['chi_angle_num']=3


if __name__ == '__main__':
    print('[ATTENTION] Please make sure you have set up move_jump, all the parameters and the score function in this script!')
    pdb_jump = int(options['jump'])
    opts=[]
    opts.append(options['energy_function'])
    if options['extra_res_fa']:
        opts.append('-extra_res_fa')
        opts.append(options['extra_res_fa'])
    if options['constraint_file']:
        opts.append('-cst_file')
        opts.append(options['constraint_file'])
    print(options)
    print(opts)
    pyrosetta.init(extra_options=' '.join(opts))
    p=pyrosetta.rosetta.core.import_pose.pose_from_file(options['input_pdb'])
    gn=get_ala_asn_from_pdb(p)
    gn.dump_pdb('origin_complex.pdb')
    if options['constraint_file']:
        pyrosetta.rosetta.core.scoring.constraints.add_constraints_from_cmdline_to_pose(gn)
    j=gn.jump(pdb_jump)
    #sample_resi = gn.residue(int(options['Resi_num']))
    print('[ATTENTION] Please check residue numbers in the constraint file!')
    sf_sys = sc_energy_fxn2()
    old_trans = V3(j.get_translation())
    old_rot = M3(j.get_rotation())
    low_rot=old_rot
    low_trans=old_trans
    old_E=sf_sys(gn)
    low_E=old_E
    print(old_E)
    print(sf_sys)
    T = 0.2
    total_cycle = int(1000)
    cart_resl = 1.0
    ang_resl = 15.0
    n_accept = 0
    hc = make_hash_store_rot(cart_resl , ang_resl, options['target_phi_psi_list'],options['orient_atoms'])
    for i in range(total_cycle):
        rand=random.random()
        k = 1
        if rand <0.5 :
            k = -1
        j.gaussian_move(k,options['Gaus_trans_mag'],options['Gaus_rot_mag'])
        gn.set_jump(pdb_jump,j)
        #sample_resi = gn.residue(int(options['Resi_num']))
        #name='test_gaussian_' + str(i) + '.pdb'
        #gn.dump_pdb(name)
        new_E = sf_sys(gn)
        print(new_E)

        delta = new_E - old_E
        accept = 1
        if delta > 0 :
            if math.exp(-delta/T) <random.random():
                accept = 0
        if accept :
            name='accept_gaussian_' + str(i) + '.pdb'
            gn.dump_pdb(name)
            print(new_E)
            old_E = new_E
            old_trans = V3(j.get_translation())
            old_rot = M3(j.get_rotation())
            if low_E > new_E:
                low_E = new_E
                low_rot = old_rot
                low_trans = old_trans

            if new_E < float(options['Energy_cut']) :
                #put it into hash
                bb_rays = []
                for num in range(3):
                    bb_rays.append(V3(gn.residue(int(options['create_ray_atoms_num'][num])).xyz(options['create_ray_atoms_name'][num])))
                #sample_resi = gn.residue(int(options['Resi_num']))
                Mark = hc.update_table_chi(bb_rays,gn.residue(int(options['Resi_num'])),options['Resi_name'],options['chi_angle_num'])
                print(str(Mark))
                n_accept = n_accept + 1
                if str(Mark) == 'True':
                    break;
        else:
            j.set_rotation(old_rot)
            j.set_translation(old_trans)
            gn.set_jump(pdb_jump,j)

print('low energy: %s'%low_E)
j.set_rotation(low_rot)
j.set_translation(low_trans)
gn.set_jump(pdb_jump,j)
rot1 = make_hash(cart_resl,ang_resl,options['target_phi_psi_list'],options['orient_atoms']).orient_rots(gn.residue(int(options['Resi_num'])), options['Resi_name'])[1]
gn.dump_pdb('finalpre.pdb')
gn.replace_residue(int(options['Resi_num']),rot1,False)
gn.dump_pdb('final.pdb')
pickle.dump(hc.dd, open(options['hashtable_name'], "wb"))

