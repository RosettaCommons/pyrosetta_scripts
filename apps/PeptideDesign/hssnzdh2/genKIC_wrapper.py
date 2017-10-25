#Authors: Parisa Hosseinzadeh
#Date: 17/10/2017
#Description: Expands and runs genKIC on chain A to generate macrocycles with a user defiend length

from __future__ import print_function

import argparse
import math
import os
import sys
import itertools
from itertools import chain
from multiprocessing import Manager, Pool, freeze_support, cpu_count

from pyrosetta import *
from rosetta import *

_DEBUG=False

def gen_kic_mover(p, num_res_added_to_n, num_res_added_to_c, n_cys, c_cys,
                  scorefxn, flank):
    
    pivot_res = list(chain.from_iterable([(p.size() -
                                        (num_res_added_to_c + flank) + 1,
                                           'CA'),
                                        (n_cys, 'CA'),
                                        (num_res_added_to_n + flank, 'CA')]))
    terminus_ranges = ((1, num_res_added_to_n + flank + 1),
                       (p.size() - (num_res_added_to_c + flank) + 1,
                        p.size() + 1))
                           
    gk = protocols.generalized_kinematic_closure.GeneralizedKIC()
    gk.set_selector_type('lowest_energy_selector')
    gk.set_selector_scorefunction(scorefxn)
    gk.set_closure_attempts(int(1E4))
    gk.set_min_solution_count(10)
                                          
    gk.add_perturber('randomize_alpha_backbone_by_rama')
    gk.set_perturber_custom_rama_table('flat_symm_dl_aa_ramatable')

    for terminus_range in reversed(terminus_ranges):
    #for res_num in reversed(range(*terminus_range)):
        for res_num in range(*terminus_range):
            gk.add_loop_residue(res_num)
            if res_num not in (n_cys, c_cys):
                # add residue to perturber
                gk.add_residue_to_perturber_residue_list(res_num)

    for res_num in pivot_res:
        if type(res_num) != int or res_num in (n_cys, c_cys):
            continue
            # The pivot residues are not necessarily in good regions of
            # Ramachandran space, so we should filter by the rama_prepro energy
            # of pivot positions.
        gk.add_filter('rama_prepro_check')
        gk.set_filter_resnum(res_num)
        gk.set_filter_rama_cutoff_energy(2.0)

    gk.add_filter('loop_bump_check')
    gk.close_bond(n_cys, 'N',
                  c_cys, 'C',
                  0, '', 0, '',  # optional params -- use default values
                  1.32,
                  114,
                  123,
                  180.,
                  False, False)

    if _DEBUG:
        print('GeneralizedKIC: {}'.format(pivot_res))
    
    gk.set_pivot_atoms(*pivot_res)
    
    '''
    gk.set_preselection_mover(_setup_preselection_mover(p, scorefxn,
                                                        terminus_ranges,
                                                        flank))
    '''
    #gk.set_preselection_mover(_setup_fast_design(p, scorefxn, terminus_ranges,flank))
        
    gk.apply(p)
    return gk.get_last_move_status()

def main(argv):
    
    parser = argparse.ArgumentParser(description='Program')
    parser.add_argument('-i', '--input', action='store', type=str,
                    required=True,
                    help='input target pdb')
    parser.add_argument('-l', '--length', action='store', type=int,
                    required=True,
                    help='max length of the loop, here peptide')
    parser.add_argument('-r', '--resn', action='store', type=str,
                    default='GLY',
                    help='residue type to append')
    parser.add_argument('-n', '--nstruct', action='store', type=int,
                    default=1,
                    help='how many times to run each KIC run')
    parser.add_argument('-c', '--chain', action='store', type=int,
                        default=1,
                        help='what chain number is the one I am extending')

    args = parser.parse_args()

    #initiaiting Rosetta
    init(extra_options='-in:file:fullatom true -mute all -write_all_connect_info -extra_res_fa /home/parisah/work/ncAA/SHA.params')
    scorefxn = get_score_function()
    
    # get the pose from target and scaffold
    p_in=rosetta.core.import_pose.pose_from_file(args.input)

    in_fname = args.input
    include_initial_termini_in_loop = True
    base_fname = in_fname.split('.')[0]
    out_put_fname = '{dir}/{base}_N-{n_add}_C-{c_add}_{numrun}.pdb'
    
    if (args.length + int(include_initial_termini_in_loop) < 3):
        print('the loop needs to be at least 3 residues', file=sys.stderr)
        sys.exit(1)

    #defining the residue I want to append
    chm = rosetta.core.chemical.ChemicalManager.get_instance()
    rts = chm.residue_type_set( 'fa_standard' )
    res = rosetta.core.conformation.ResidueFactory.create_residue(rts.name_map(args.resn))
            
    for loop in range (3,args.length+1):
        out_dir = '{}_genKIC_{}'.format(base_fname,loop)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        for runs in range (0,args.nstruct):
            p=Pose()
            for resNo in rrange(p_in.size()):
                if (p_in.residue(resNo).chain() == args.chain):
                    p.append_residue_by_bond(p_in.residue(resNo),False)
            num=p.size()+1
            for ir in rrange(p.size()):
                if ( p.residue(ir).has_variant_type(core.chemical.UPPER_TERMINUS_VARIANT)):
                    core.pose.remove_variant_type_from_pose_residue( p, core.chemical.  UPPER_TERMINUS_VARIANT, ir)
                if ( p.residue(ir).has_variant_type(core.chemical.LOWER_TERMINUS_VARIANT)):
                    core.pose.remove_variant_type_from_pose_residue( p, core.chemical.LOWER_TERMINUS_VARIANT, ir)
                if ( p.residue(ir).has_variant_type(core.chemical.CUTPOINT_LOWER)):
                        core.pose.remove_variant_type_from_pose_residue( p, core.chemical.CUTPOINT_LOWER, ir)
                if ( p.residue(ir).has_variant_type(core.chemical.CUTPOINT_UPPER)):
                        core.pose.remove_variant_type_from_pose_residue( p, core.chemical.CUTPOINT_UPPER, ir)
            if _DEBUG:
                p.dump_pdb('pose_clone.pdb')
        
            # setting residues to prepend to N-term and setting omega to 180
            N_add=loop//2
            for i in range(0, N_add):
                p.prepend_polymer_residue_before_seqpos(res, 1, True)
            if _DEBUG:
                p.dump_pdb('prepend_test.pdb')
            for res_no in range(1,N_add+1):
                p.set_omega(res_no,180.)
            # setting residues to append to C-term and setting omega to 180
            C_add=loop-N_add
            for ir in rrange(p.size()):
                if ( p.residue(ir).has_variant_type(core.chemical.UPPER_TERMINUS_VARIANT)):
                    core.pose.remove_variant_type_from_pose_residue( p, core.chemical.UPPER_TERMINUS_VARIANT, ir)
                if ( p.residue(ir).has_variant_type(core.chemical.LOWER_TERMINUS_VARIANT)):
                    core.pose.remove_variant_type_from_pose_residue( p, core.chemical.LOWER_TERMINUS_VARIANT, ir)
                if ( p.residue(ir).has_variant_type(core.chemical.CUTPOINT_LOWER)):
                    core.pose.remove_variant_type_from_pose_residue( p, core.chemical.CUTPOINT_LOWER, ir)
                if ( p.residue(ir).has_variant_type(core.chemical.CUTPOINT_UPPER)):
                    core.pose.remove_variant_type_from_pose_residue( p, core.chemical.CUTPOINT_UPPER, ir)
            for i in range(0,C_add):
                p.append_residue_by_bond(res, True)
            if _DEBUG:
                p.dump_pdb('append_test.pdb')
            for res_no in range((p.size()-C_add)-1, p.size()+1):
                p.set_omega(res_no,180.)
            # declaring the bond before genKIC call
            to_close=(1,p.size())
            pcm=protocols.cyclic_peptide.PeptideCyclizeMover()
            pcm.apply(p)
            if _DEBUG:
                p.dump_pdb('bonded.pdb')
            # calling genKIC mover from Brian's script
            st=gen_kic_mover(p,N_add,C_add,to_close[0],to_close[1],scorefxn,int(include_initial_termini_in_loop))
            if st == protocols.moves.MoverStatus.MS_SUCCESS:
                p_fin=Pose()
                for resi in rrange(p.size()):
                    p_fin.append_residue_by_bond(p.residue(resi),False)
                for resNo in rrange(p_in.size()):
                    if (p_in.residue(resNo).chain() != args.chain):
                        if (resNo == num):
                            p_fin.append_residue_by_jump(p_in.residue(resNo),p.size(),'','',True)
                        else:
                            p_fin.append_residue_by_bond(p_in.residue(resNo),False)
                db = protocols.cyclic_peptide.DeclareBond()
                db.set(to_close[0],'N',to_close[1],'C',False,False,0,0,True)
                db.apply(p_fin)
                p_fin.dump_pdb(out_put_fname.format(dir=out_dir,
                                    base=base_fname,
                                    n_add=N_add,
                                    c_add=C_add,
                                    numrun=runs+1))
            else:
                if _DEBUG:
                    print ("no successful solutions found")


if __name__ == '__main__':
    freeze_support()
    main(sys.argv)
