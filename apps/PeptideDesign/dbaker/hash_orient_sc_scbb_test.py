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
# IMPORTANT: you can't test anything in this repo that depends on tools from outside it!
####import rif.hash
from hash_subclass_sc_scbb import *
import copy


## Use Monte Carlo to generate look up tables for low energy interactions (eg bidentate hbonds) between a sidechain in "sc_res" and a sidechain+backbone in "scbb_res".
## Must provide an input pdb file with sc_res as a single residue and scbb_res in the 2nd position of a tripeptide, and specify a pair of interacting heavy atoms for
## jump between the two (sc_atom and scbb_jump_atom).  Also must specify the torsion degrees of freedom of scbb_res that influence the interaction geometry--some subset
## of phi, psi, chi1 and chi2  (scbb_torsions).  These will be explicitly varied in the Monte Carlo along with the rigid body degrees of freedom.  The hash key variables are the rigid body transform between the backbones of sc_res and scbb_res and whichever of phi and psi of scbb_res are in scbb_torsions.  For efficiency, the inverse rotamers of sc_res
## and their backbone frames are generated at the beginning and then kept fixed during the MC.  By default, what is stored in hash are the chi angles of sc_res, and any chi angles## in scbb_torsions.  Can additional store residue names for case where hash contains multiple residue pairs ("store_res_names_in_hash")
# the input pdb should have two chains.  one chain is a single residue contributing the sidechain interactions (sc_res) and the other a tripeptide with its 2nd residue contributing the sidechain+backbone interactions (scbb_res); specify the residue number of sc_res in the input file using sc_resN


def get_options():
    #set up the parameters
    options={}
#required inputs
    options['input_pdb']='asn_bb_c.pdb'
    options['sc_resN']=4
    options['scbb_resN']=2
    options['hash_function']= rif.hash.RosettaStubTwoTorsionHash(phi_resl=10, cart_resl=1, ori_resl=10, cart_bound=32)
    options['orient_atoms'] = ['OD1','ND2','CG']
    ## use this if interaction only depends on one torsion    rif.hash.RosettaStubTorsionHash(phi_resl=15, cart_resl=1, ang_resl=11)
## use this if interaction depends on no torsions (ie, sidechain-sidechain) rif.hash.RosettaStubTorsionHash(cart_resl=1, ang_resl=11)
    options['scbb_torsions']=['phi','psi']  # these are the torsions that affect hbond geometry and will be varied in MC
    options['scbb_jump_atom']='O'
    options['sc_jump_atom']='ND2'
    options['scbb_torsions']= {'phi': {'interval': float(5.), 'nsamples': int(21) }, 'psi' : {'interval': float(5.) , 'nsamples': int(21) }}
#required output specifications
    options["store_res_names_in_hash"]= 1
    options['hashtable_name']='ASN_BB_10_1_10_5_21_w_resNAME2'
#the below don't have to be set on a case by case basis
    options['target_phi_psi_list'] = ['helical','sheet']
    options['temperature']=float(0.2)
    options['MC_cycles']=4000
    options['dump_interval']=1000  # write out structures for every nth addition to hash
    options['Gaus_trans_mag']=float(0.25)
    options['Gaus_rot_mag']=float(1.5)
    options['thresholdE_for_hash']=-0.2
    options['energy_function']='-empty'
    options['constraint_file']=''
    options['extra_res_fa']=''
    options['jump']=1
    #   options['create_stub_atoms_num']=[3,3,3]
 #   options['create_stub_atoms_name']=['N','CA','C']
    
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
    return options,opts

def get_argparse():
    parser = argparse.ArgumentParser(description='Generate hash table by MC sampling of interaction between sidechain and peptide')
    parser.add_argument('--input_pdb_file', type=str, dest='input_pdb',
                   default=options['input_pdb'],
                   help='input structure with example of interaction')
    parser.add_argument('--MC_cycles', type=int, dest='total_cycles',
                   default=int(options['MC_cycles']),
                   help='number of monte carlo cycles')
    parser.add_argument('--dump_interval', type=int, dest='dump_interval',
                   default=options['dump_interval'],
                   help='how often to dump pdbs')
    parser.add_argument('--hash_file_name',type=str,dest='hashtable_name',
                        default=options['hash_table_name'],
                        help='name of generated hash file')
    return parser

name_data_type=namedtuple('name_data_type',['sc_resName','scbb_resName'])


def get_start_pdb(input_pdb,sc_jump_atom,scbb_jump_atom,sc_resN,scbb_resN):
    p=pyrosetta.rosetta.core.import_pose.pose_from_file(options['input_pdb'])
#    names=namedtuple('names',['sc_resName','scbb_resName'])
    name_data=name_data_type(p.residue(sc_resN).name1(),p.residue(scbb_resN).name1())
    chains=p.split_by_chain()
    if sc_resN==1:
        sc_res=chains[1]
        scbb_res=chains[2]
    else:
        sc_res=chains[2]
        scbb_res=chains[1]
    print(sc_res,scbb_jump_atom,scbb_res,sc_jump_atom)    
    scbb_res.append_pose_by_jump(sc_res,2,scbb_jump_atom,sc_jump_atom)  # sc_res appended by jump to 2nd residue of scbb_res
    ft=scbb_res.fold_tree()
    print('reordering fold tree')
    ft.reorder(sc_resN)
    scbb_res.fold_tree(ft)  #fold tree rooted on sc_res so it and its inverse rotamers stay fixed during MC
    print('done reordering fold tree')
    return scbb_res,name_data

def get_tor_combos_to_sample(start_tor_vals, Nsamples_per_tors, Tors_sample_interval):
# generate all combinations of values of tors_vars
    combs=[]
    Ntors=len(start_tor_vals)
    Ncombs=1
    for i in range(Ntors):
        Ncombs=Ncombs*Nsamples_per_tors[i]

    print(Ncombs)
    for i in range(Ncombs):
        combs.append([])
    remain=Ncombs    
    for i in range(Ntors):
        remain=remain/Nsamples_per_tors[i]
        iter=0
        state=0
        for j in range(Ncombs):
# need to make so start with input values and sample outwards (so keep low energy start of each MC run)
# rather than sampling from -x to +x, would be better to go from 0 to +x, and then 0 to -x
            if state < Nsamples_per_tors[i]/2:
              delta= Tors_sample_interval[i]*state
            else:
              delta= -1.* (state-Nsamples_per_tors[i]/2)*Tors_sample_interval[i]
            combs[j].append(start_tor_vals[i]+ delta)
            iter=iter+1
            if iter == remain:
                state=state+1
                if state == Nsamples_per_tors[i]:
                    state=0
                iter=0
    return combs

def update_pose_tors(pose,tor_ids,resN,combo):
    for i in range(len(tor_ids)):
        if tor_ids[i] == 'phi':
            pose.set_phi(resN,combo[i])
        else:
            if tor_ids[i]=='psi':
                pose.set_psi(resN,combo[i])
            elif tor_ids[i]=='chi1':
                pose.set_chi(1,resN,combo[i])  
            elif tor_ids[i]=='chi2':
                pose.set_chi(2,resN,combo[i])
            else:
                print('WARNING TOR_ID NOT FOUND!')
    return

def get_tors_start_vals(pose,tor_ids,resN):
    start_vals=[]
    for i in range(len(tor_ids)):
        if tor_ids[i] == 'phi':
            start_vals.append(gn.phi(resN))
        elif tor_ids[i]=='psi':
            start_vals.append(gn.psi(resN))
        elif tor_ids[i]=='chi1':
            start_vals.append(gn.residue(resN).chi(1))
        elif tor_ids[i]=='chi2':
            start_vals.append(gn.residue(resN).chi(2))
        else:
            print('WARNING TOR_ID NOT FOUND!')
    return start_vals

#set up score function
def hb_energy_fxn():
    sf = pyrosetta.get_score_function()
# do print(sf) and look for weights to see what you can modify
    sf.set_weight(core.scoring.fa_dun, 0.3)
    sf.set_weight(core.scoring.fa_atr, 0.3)
    sf.set_weight(core.scoring.fa_rep, 0.1)
    sf.set_weight(core.scoring.fa_sol, 0.2)
    sf.set_weight(core.scoring.fa_intra_rep, 0.0)
    sf.set_weight(core.scoring.fa_elec, 0.2)
    sf.set_weight(core.scoring.hbond_bb_sc, 2.0)
    sf.set_weight(core.scoring.hbond_lr_bb, 0.0)
    sf.set_weight(core.scoring.hbond_sr_bb, 0.0)
    sf.set_weight(core.scoring.hbond_sc, 2.0)
#    sf.set_weight(core.scoring.fa_sr_bb, 0.0)
    return sf

if __name__ == '__main__':
    # AMW NOTE: We can't actually use this as a test because it imports "rif.hash" and that isn't included
    # in our current distribution of pyrosetta_cripts?
    pass
    options, opts=get_options()
    pyrosetta.init(extra_options=' '.join(opts))
    pdb_jump = int(options['jump'])
    sc_resN=int(options['sc_resN'])
    scbb_resN=int(options['scbb_resN'])
    T = options['temperature']
    thresholdE_for_hash=float(options['thresholdE_for_hash'])
    gn,name_data=get_start_pdb(input_pdb,options['sc_jump_atom'], options['scbb_jump_atom'],sc_resN,scbb_resN)
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

    dump_it=0
    n_accept = 0
    n_accept_below_threshold=0
    # pass in sc_res at initialization so can generate its inverse rotamers and frames, which will stay fixed through MC
    hasher=options['hash_function']
    hc=make_hash_sc_scbb(hasher, options['scbb_torsions'], options['orient_atoms'], options['target_phi_psi_list'], sc_res,options['store_res_names_in_hash'],name_data)

    Nsamples_per_tors_list=[]
    Tors_sample_interval_list=[]
    tor_info=options['scbb_torsions']
    tor_list=list(tor_info) 
    for tor in tor_list:
        Nsamples_per_tors_list.append(tor_info[tor]['nsamples'])
        Tors_sample_interval_list.append(tor_info[tor]['interval'])
    start_vals=get_tors_start_vals(gn,tor_list,scbb_resN)
    print(tor_list)
    print('starting torsion values:', start_vals)
    gn.dump_pdb('before_mc.pdb')
    combos=get_tor_combos_to_sample(start_vals, Nsamples_per_tors_list, Tors_sample_interval_list)
    tot_combos=len(combos)
    tr=namedtuple('tr',['trans','rot'])
    ori_jump=tr(trans=old_trans,rot=old_rot)
    low_jumps=[]
    low_jumps.append( ori_jump)
    low_trans=old_trans
    low_rot=old_rot
#need to move smoothly in torsion space to avoid MC blowing up (can't lose h-bond completely)
    for icombo in range(tot_combos):
         print('adding to hash for torsion angle combo: ',combos[icombo])   
    #     print('gn before update',gn,'jump1',gn.jump(1))
         update_pose_tors(gn,tor_list,scbb_resN,combos[icombo])
    #     print('gn after update',gn,'jump1',gn.jump(1))
    #     gn.dump_pdb('after_tors_update.pdb')
         E_after_update=sf_sys(gn)
         low_E=100

    #to increase efficiency of MC, start with best jump found so far for new torsion angle combo
         for trans_rot in low_jumps:
             j.set_translation( trans_rot.trans)
             j.set_rotation(trans_rot.rot)
             gn.set_jump(pdb_jump,j)
             new_E=sf_sys(gn)
             print('replacing jump',trans_rot.trans,trans_rot.rot,new_E,low_E)
             if new_E < low_E:
                 low_E=new_E
                 low_trans=trans_rot.trans
                 low_rot=trans_rot.rot
         j.set_translation(low_trans)
         j.set_rotation(low_rot)
         gn.set_jump(pdb_jump,j)
         print('energy after torsion update',E_after_update,'energy with best jump',low_E)
         low_E=100
         for i in range(total_cycle):
            rand=random.random()
            k = 1
            j.gaussian_move(k,options['Gaus_trans_mag'],options['Gaus_rot_mag'])
            gn.set_jump(pdb_jump,j)
            new_E = sf_sys(gn)
            delta = new_E - old_E

            accept = 1
            if delta > 0 :
                if math.exp(-delta/T) <random.random():
                    accept = 0

            if accept :
                n_accept=n_accept+1
                print('accept %s %s %s %s %s %s '%(n_accept,n_accept_below_threshold,new_E,old_E,low_E,delta))
     #           name='accept_gaussian_' + str(i) + '.pdb'
     #           gn.dump_pdb(name)
     #           print(new_E)
                old_E = new_E
                old_trans = V3(j.get_translation())
                old_rot = M3(j.get_rotation())

                if low_E > new_E:
                    low_E = new_E
                    low_rot = old_rot
                    low_trans = old_trans

                if new_E < thresholdE_for_hash :
                    hc.update_hash(gn,scbb_resN)
                    n_accept_below_threshold = n_accept_below_threshold + 1

                    if dump_it> dump_interval:
                        gn.dump_pdb('t%s.pdb'%i)
                        dump_it=0


            else:
                j.set_rotation(old_rot)
                j.set_translation(old_trans)
                gn.set_jump(pdb_jump,j)
            dump_it=dump_it+1
            if i == 500 and low_E > 1.2:
                print('Stopped sampling because in outer space: ',combos[icombo])
                break     ##  if haven't found good interactions by this time, don't waste more

         print('tors combo and low_E: ',combos[icombo],low_E,n_accept,n_accept_below_threshold)
         low_jump=tr(trans=low_trans,rot=low_rot)   
         low_jumps.append(low_jump)
     
        
# AMW NOTE: We can't actually use this as a test because it imports "rif.hash" and that isn't included
# in our current distribution of pyrosetta_cripts?
#print('low energy: %s'%low_E)
#j.set_rotation(low_rot)
#j.set_translation(low_trans)
#gn.set_jump(pdb_jump,j)
##rot1 = make_hash(cart_resl,ang_resl,options['target_phi_psi_list'],options['orient_atoms']).orient_rots(gn.residue(int(options['Resi_num'])), options['Resi_name'])[1]
#gn.dump_pdb('finalpre.pdb')
##gn.replace_residue(int(options['Resi_num']),rot1,False)
##gn.dump_pdb('final.pdb')
pickle.dump(hc.dd, open(hashtable_name, "wb"))

