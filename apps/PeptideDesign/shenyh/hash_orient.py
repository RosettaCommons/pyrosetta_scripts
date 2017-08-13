#!/home/sheffler/venv/david/bin/python
import pyrosetta
import argparse
import math
import random
import collections
import pyrosetta.rosetta
from pyrosetta.rosetta import *
from pyrosetta.rosetta.numeric import xyzVector_double_t as V3
from pyrosetta.rosetta.numeric import xyzMatrix_double_t as M3
import pickle
from pyrosetta.rosetta.protocols.protein_interface_design.movers import TryRotamers
from rif.geom.ray_hash import RayRay10dHash
from rif.geom import Ray
from hash_subclass import *
import sys

def get_ala_asn_from_pdb(p):
    q=move_jump(p.clone())
    return q
#set up the move_jump
def move_jump(p):
    chains=p.split_by_chain()
    Surf=chains[1]
    Resi=chains[2]
    #change here
    Surf.append_pose_by_jump(Resi,2,'O3','CG')
    return Surf
#set up score function
def sc_energy_fxn2():
    sf = pyrosetta.get_score_function()
    sf.set_weight(core.scoring.fa_sol, 0.1)
    sf.set_weight(core.scoring.fa_rep, 0.55)
    sf.set_weight(core.scoring.fa_elec, 1.0)
    sf.set_weight(core.scoring.atom_pair_constraint, 1.0)
    sf.set_weight(core.scoring.angle_constraint, 1.0)
    sf.set_weight(core.scoring.dihedral_constraint, 1.0)
    return sf



def main( ReplacementArgs=[] ):
    # Build argument string with arguments passed to main from testing script
    if len(ReplacementArgs):
        if type(ReplacementArgs) == str:
            ReplacementArgs = list(ReplacementArgs.split())
        # set sys.argv as if passed from commandline 
        sys.argv = ReplacementArgs[:]

    print '[ATTENTION] Please make sure you have set up move_jump, all the parameters and the score function in this script!'
    print 'sys.argv:'
    print sys.argv

    #set up the parameters
    # Should be arguements:
    argument_parser = argparse.ArgumentParser(description=' Arguments for generating rotamer hashes with hash_orient.py' )
    argument_parser.add_argument('-pdb', type=str, help='Input pdb file', required = True )
    argument_parser.add_argument('-cst', default=False,   type=str, help='constraint file describing residue interactions to add to hash' )
    argument_parser.add_argument('-score_function', default='empty',   type=str, help='starting score function (will be modified)' )
    argument_parser.add_argument('-params', default=False,    type=str, nargs="+", help='params files' )
    argument_parser.add_argument('-Gaus_trans_mag', default=0.15,   type=float, help='magnitude of Gaussian translation moves' )
    argument_parser.add_argument('-Gaus_rot_mag', default=1.5,   type=float, help='magnitude of Gaussian rotation moves' )
    argument_parser.add_argument('-Energy_cut', default=-50.0,    type=float, help=' Energy cutoff for montecarlo rotamer cloud sampling' )
    argument_parser.add_argument('-cycles', default=100000,  type=int, help='Number of Monte Carlo cycles for rotamer sampling' )
    # Maybe should be arguments
    argument_parser.add_argument('-rotamer_phi_psi_list', type=list, help='' , default = [(-118.5, 118.3), (-77.6, 175.9), (-63.6, 148.9), (63.2, 17.1), (74.2, 9.9), (-105.6, 110.3), (-64.0, 132.3), (69.4, -2.7), (-90.3, 107.8), (-55.3, 120.8)] )
    argument_parser.add_argument('-hashtable_name', type=str, help='' , default = 'CaCO3_hashtable' )
    # Should not be arguments, values should be derived from input files or hashes
    argument_parser.add_argument('-jump', type=int, default=4,   help='jump to sample rigid body movements though' )
    argument_parser.add_argument('-Resi_num', default=5,   type=int, help='' )
    argument_parser.add_argument('-Resi_name', default='ASP',   type=str, help='' )
    argument_parser.add_argument('-orient_atoms', type=str, nargs="+", help='' , default = ['OD1','CG','OD2'] )
    argument_parser.add_argument('-create_ray_atoms_num', type=str, nargs="+", help='' , default = [1,2,2,4] )
    argument_parser.add_argument('-create_ray_atoms_name', type=str, nargs="+", help='' , default = ['Ca2p','C1','O1','O1'] )
    argument_parser.add_argument('-chi_angle_num', type=int, help='chi_angle_num' , default = 2 )
    options = argument_parser.parse_args()

    pdb_jump = int(options.jump)
    pyrosetta_opts=[]
    pyrosetta_opts.append( '-{0}'.format(options.score_function) )
    if options.params:
        pyrosetta_opts.append('-extra_res_fa')
        pyrosetta_opts.extend(options.params)
    if options.cst:
        pyrosetta_opts.append('-cst_file')
        pyrosetta_opts.append(options.cst)
    

    pyrosetta.init(extra_options=' '.join(pyrosetta_opts))
    p=pyrosetta.rosetta.core.import_pose.pose_from_file(options.pdb)
    gn=get_ala_asn_from_pdb(p)
    gn.dump_pdb('origin_complex.pdb')
    if options.cst:
        pyrosetta.rosetta.core.scoring.constraints.add_constraints_from_cmdline_to_pose(gn)
    j=gn.jump(pdb_jump)
    # sample_resi = gn.residue(int(options.Resi_num))
    print '[ATTENTION] Please check residue numbers in the constraint file!'
    sf_sys = sc_energy_fxn2()
    print sf_sys

    old_trans = V3(j.get_translation())
    old_rot = M3(j.get_rotation())
    low_rot=old_rot
    low_trans=old_trans
    old_E=sf_sys(gn)
    low_E=old_E
    print(old_E)
    print(sf_sys)
    T = 0.2
    # total_cycle = int(100000)
    resl = 0.4
    lever = 3.0
    n_accept = 0
    print make_hash_store_rot.__module__
    hc = make_hash_store_rot(resl, lever, options.rotamer_phi_psi_list, options.orient_atoms )
    for i in range(options.cycles):
        rand=random.random()
        k = 1
        if rand <0.5 :
            k = -1
        j.gaussian_move(k, options.Gaus_trans_mag, options.Gaus_rot_mag)
        gn.set_jump(pdb_jump,j)
        # sample_resi = gn.residue(int(options.Resi_nums))
        #name='test_gaussian_' + str(i) + '.pdb'
        #gn.dump_pdb(name)
        new_E = sf_sys(gn)
        print new_E

        delta = new_E - old_E
        accept = 1
        if delta > 0 :
            if math.exp(-delta/T) <random.random():
                accept = 0
        if accept :
            name='accept_gaussian_' + str(i) + '.pdb'
            gn.dump_pdb(name)
            # print type(gn)
            # print gn
            # print 'phi, psi:', gn.phi(5), gn.psi(5), '|{0}|'.format(gn.sequence())
            # print new_E
            old_E = new_E
            old_trans = V3(j.get_translation())
            old_rot = M3(j.get_rotation())
            if low_E > new_E:
                low_E = new_E
                low_rot = old_rot
                low_trans = old_trans

            if new_E < float(options.Energy_cut) :
                #put it into hash
                bb_rays = []
                for num in range(4):
                    bb_rays.append(V3(gn.residue(int(options.create_ray_atoms_num[num])).xyz(options.create_ray_atoms_name[num])))
                #sample_resi = gn.residue(int(options['Resi_num']))
                Mark = hc.update_table_chi(bb_rays,gn.residue(int(options.Resi_num)), options.Resi_name, options.chi_angle_num)
                print str(Mark)
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
    rot1 = make_hash(resl, lever, options.rotamer_phi_psi_list, options.orient_atoms).orient_rots(gn.residue(int(options.Resi_num)), options.Resi_name)[1]
    gn.dump_pdb('finalpre.pdb')
    gn.replace_residue(int(options.Resi_num),rot1,False)
    gn.dump_pdb('final.pdb')
    pickle.dump(hc.dd, open(options.hashtable_name, "wb"))


if __name__ == '__main__':
   sys.exit(main())

