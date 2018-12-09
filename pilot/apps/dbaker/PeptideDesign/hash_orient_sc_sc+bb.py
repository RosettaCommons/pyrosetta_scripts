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
import rif.hash
from hash_subclass_sc_sc+bb import *

## Use Monte Carlo to generate look up tables for low energy interactions (eg bidentate hbonds) between a sidechain in "sc_res" and a sidechain+backbone in "scbb_res".
## Must provide an input pdb file with sc_res as a single residue and scbb_res in the 2nd position of a tripeptide, and specify a pair of interacting heavy atoms for
## jump between the two (sc_atom and scbb_jump_atom).  Also must specify the torsion degrees of freedom of scbb_res that influence the interaction geometry--some subset
## of phi, psi, chi1 and chi2  (scbb_torsions).  These will be explicitly varied in the Monte Carlo along with the rigid body degrees of freedom.  The hash key variables are the rigid body transform between the backbones of sc_res and scbb_res and whichever of phi and psi of scbb_res are in scbb_torsions.  For efficiency, the inverse rotamers of sc_res
## and their backbone frames are generated at the beginning and then kept fixed during the MC.  By default, what is stored in hash are the chi angles of sc_res, and any chi angles## in scbb_torsions.  Can additional store residue names for case where hash contains multiple residue pairs ("store_res_names_in_hash")
# the input pdb should have two chains.  one chain is a single residue contributing the sidechain interactions (sc_res) and the other a tripeptide with its 2nd residue contributing the sidechain+backbone interactions (scbb_res); specify the residue number of sc_res in the input file using sc_resN
.   
def get_options:
    #set up the parameters
    options={}
#required inputs
    options['input_pdb']='asn_bb_c.pdb'
    options['sc_resN']=1
    options['scbb_resN']=3
    options['hash_function']= rif.hash.RosettaStubTwoTorsionHash(phi_resl=15, cart_resl=1, ang_resl=11)
## use this if interaction only depends on one torsion    rif.hash.RosettaStubTorsionHash(phi_resl=15, cart_resl=1, ang_resl=11)
## use this if interaction depends on no torsions (ie, sidechain-sidechain) rif.hash.RosettaStubTorsionHash(cart_resl=1, ang_resl=11)
    options['scbb_torsions']=['phi','psi']  # these are the torsions that affect hbond geometry and will be varied in MC
    options['scbb_jump_atom']='O'
    options['sc_jump_atom']='ND2'
    options['scbb_torsions']= {'phi': {'interval': float(2.5), 'nsamples': int(20) }, 'psi' : {'interval': float(2.5) , 'nsamples': int(20) }}
#required output specifications
    options["store_res_names_in_hash"]= 0
    options['hashtable_name']='ASN_BB'
#the below don't have to be set on a case by case basis    
    options['target_phi_psi_list'] = ['helical','sheet']
    options['create_stub_atoms_num']=[3,3,3]
    options['create_stub_atoms_name']=['N','CA','C']
    options['temperature']=float(0.2)
    options['MC_cycles']=1000
    options['energy_function']='-empty'
    options['constraint_file']=''
    options['extra_res_fa']=''
    options['jump']=1
    options['Gaus_trans_mag']=float(0.25)
    options['Gaus_rot_mag']=float(1.5)
    options['Energy_cut']=-5
    options['orient_atoms'] = ['N','CA','C']
    return options


def get_start_pdb(input_pdb,sc_jump_atom,scbb_jump_atom,sc_resN,scbb_resN):
    p=pyrosetta.rosetta.core.import_pose.pose_from_file(options['input_pdb'])
    names=named_tuple('names',['sc_resN','scbb_resN'])
    name_data=named_tuple(pose.residue(sc_resN).name1(),pose.residue(scbb_resN).name1())
    chains=p.split_by_chain()
    if sc_resN==1:
        sc_res=chains[1]
        scbb_res=chains[2]
    else:
        sc_res=chains[2]
        scbb_res=chains[1]
    scbb_res.append_pose_by_jump(sc_res,2,scbb_jump_atom,sc_jump_atom)  # sc_res appended by jump to 2nd residue of scbb_res
    ft=scbb_res.fold_tree()
    ft.reorder(sc_resN)
    pose.fold_tree(ft)  #fold tree rooted on sc_res so it and its inverse rotamers stay fixed during MC
    return scbb_res,name_data

def get_tor_combos_to_sample(start_tor_vals, Nsamples_per_tors, Tors_sample_interval):
# generate all combinations of values of tors_vars
    combs=[]
    Ntors=len(start_tor_vals)
    Ncombs=1
    for i in range(Ntors):
        Ncombs=Ncombs*Nsamples_per_tors[i]

    print Ncombs
    for i in range(Ncombs):
        combs.append([])
    remain=Ncombs    
    for i in range(Ntors):
        remain=remain/Nsamples_per_tors[i]
        iter=0
        state=0
        for j in range(Ncombs):
	    combs[j].append(start_tor_vals[i]+((state-(Nsamples_per_tors[i]-1)/2)*Tors_sample_interval[i]))
            iter=iter+1
            if iter == remain:
                state=state+1
                if state == Nsamples_per_tors[i]:
                    state=0
                iter=0
    return combs

def update_pose_tors(pose,tor_ids,resN,combo):
    for i in len(tor_ids):
        if tor_ids[i] == 'phi':
            pose.set_phi(resN,combo[i])
        else:
            if tor_ids[i]=='psi':
                pose.set_psi(resN,combo[i])
            else if:
                tor_ids[i]=='chi1':
                    pose.set_chi(1,resN,combo[i])  
                else if:
                    tor_ids[i]=='chi2':
                        pose.set_chi(2,resN,combo[i])
                else:
                    print('WARNING TOR_ID NOT FOUND!')
    return

def get_tors_start_vals(pose,tor_ids,resN):
    start_vals=[]
    for i in len(tor_ids):
        if tor_ids[i] == 'phi':
            start_vals.append(gn(resN).phi())
        else:
            if tor_ids[i]=='psi':
                start_vals.append(gn(resN).psi())j
             else if:
                tor_ids[i]=='chi1':
                    start_vals.append(gn.residue(resN).chi(1))
                else if:
                    tor_ids[i]=='chi2':
                    start_vals.append(gn.residue(resN).chi(2))
                 else:
                    print('WARNING TOR_ID NOT FOUND!')
    return start_vals

#set up score function
def hb_energy_fxn():
    sf = pyrosetta.get_score_function()
# do print(sf) and look for weights to see what you can modify
    sf.set_weight(core.scoring.fa_dun, 0.5)
    sf.set_weight(core.scoring.fa_atr, 0.0)
    sf.set_weight(core.scoring.fa_rep, 1.0)
    sf.set_weight(core.scoring.fa_sol, 0.2)
    sf.set_weight(core.scoring.fa_intra_rep, 0.0)
    sf.set_weight(core.scoring.fa_elec, 0.5)
    sf.set_weight(core.scoring.hbond_bb_sc, 2.0)
    sf.set_weight(core.scoring.hbond_lr_bb, 0.0)
    sf.set_weight(core.scoring.hbond_sr_bb, 0.0)
    sf.set_weight(core.scoring.hbond_sc, 2.0)
#    sf.set_weight(core.scoring.fa_sr_bb, 0.0)
    return sf

if __name__ == '__main__':
    options=get_options()
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
    pdb_jump = int(options['jump'])
    sc_resN=int(options['sc_resN'])
    scbb_resN=int(options['scbb_resN'])
    gn,name_data=get_start_pdb(options['input_pdb'],options['sc_jump_atom'], options['scbb_jump_atom'],sc_resN,scbb_resN)
    gn.dump_pdb('re_rooted.pdb')
    if options['constraint_file']:
        pyrosetta.rosetta.core.scoring.constraints.add_constraints_from_cmdline_to_pose(gn)
    j=gn.jump(pdb_jump)
    sc_res=gn.residue(sc_resN)
    sf_sys = hb_energy_fxn()
    old_trans = V3(j.get_translation())
    old_rot = M3(j.get_rotation())
    low_rot=old_rot
    low_trans=old_trans
    old_E=sf_sys(gn)
    low_E=old_E
    print('Starting energy',old_E)
    T = options['temperature']
    thresholdE_for_hash=float(options['thresholdE_for_hash'])
    total_cycle = int(options['MC_cycles'])
    dump_interval=total_cycle*0.1
    dump_it=0
    n_accept = 0
    # pass in sc_res at initialization so can generate its inverse rotamers and frames, which will stay fixed through MC
    hc=make_hash_sc_scbb(options['hash_function'], options['scbb_torsions'], options['store_res_names_in_hash'], options['orient_atoms'], options['target_phi_psi_list'], sc_res,name_data):

    Nsamples_per_tors_list=[]
    Tors_sample_interval_list=[]
    tor_info=options['scbb_torsions']
    tor_list=tor_info.keys()
    for tor in tor_list:
        Nsamples_per_tors_list.append(tor_info[tor]['nsamples'])
        Tors_sample_interval_list.append(tor_info[tor]['interval'])
    start_vals=get_tors_start_vals(gn,tor_list,scbb_resN)
    combos=get_tor_combos_to_sample(start_vals, Nsamples_per_tors_list, Tors_sample_interval_list)
    tot_combos=len(combos)
    for j in range(tot_combos):
     update_pose_tors(gn,tor_list,scbb_resN,combos[j]):
     for i in range(total_cycle):
        rand=random.random()
        k = 1
#        if rand <0.5 :
#            k = -1
        j.gaussian_move(k,options['Gaus_trans_mag'],options['Gaus_rot_mag'])
        gn.set_jump(pdb_jump,j)
 
        new_E = sf_sys(gn)
        print(new_E)
        delta = new_E - old_E

        accept = 1
        if delta > 0 :
            if math.exp(-delta/T) <random.random():
                accept = 0

        if accept :
            print('accept %s %s %s %s %s '%(n_accept,new_E,old_E,low_E,delta))
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

            if new_E < thresholdE_for_hash :
                hc.update_table(gn,scbb_resN)
                n_accept = n_accept + 1

                if dump_it> dump_interval:
                    gn.dump_pdb('t%s.pdb'%i)
                    dump_it=0
        

        else:
            j.set_rotation(old_rot)
            j.set_translation(old_trans)
            gn.set_jump(pdb_jump,j)

        dump_it=dump_it+1

print('low energy: %s'%low_E)
#j.set_rotation(low_rot)
#j.set_translation(low_trans)
#gn.set_jump(pdb_jump,j)
#rot1 = make_hash(cart_resl,ang_resl,options['target_phi_psi_list'],options['orient_atoms']).orient_rots(gn.residue(int(options['Resi_num'])), options['Resi_name'])[1]
#gn.dump_pdb('finalpre.pdb')
#gn.replace_residue(int(options['Resi_num']),rot1,False)
#gn.dump_pdb('final.pdb')
pickle.dump(hc.dd, open(options['hashtable_name'], "wb"))

