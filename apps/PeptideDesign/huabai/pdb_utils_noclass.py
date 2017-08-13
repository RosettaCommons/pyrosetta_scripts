#!/work/sheffler/anaconda2/bin/ipython
from pyrosetta import rosetta,init,Pose
#from rosetta.core.io.silent import SilentFileData as SilentFileData
from rosetta.core.io.silent import SilentFileData, SilentFileOptions
from rosetta.protocols.cyclic_peptide import *

def init_pyrosetta():
    init("-mute all  -ignore_unrecognized_res")

def input_silent(filename):
    hbond_list={}
    torsion_list={}
#    sfd=SilentFileData()
    sfd = SilentFileData(SilentFileOptions())
    sfd.read_file(filename)
    p=Pose()
    sf_tags=sfd.tags()
    for tag in sf_tags:
        silent_struct=sfd.get_structure(tag)
        silent_struct.fill_pose(p)
        torsion_list[tag]=get_torsions(p)
        hbond_list[tag]=find_hbonds(p)
    return torsion_list,hbond_list

def input_silent_score_seq(filename):
    hbond_list={}
    torsion_list={}
    score_list={}
    seq_list={}
   #  sfd=SilentFileData()
    sfd = SilentFileData(SilentFileOptions())
    sfd.read_file(filename)
    p=Pose()
    sf_tags=sfd.tags()
    for tag in sf_tags:
        silent_struct=sfd.get_structure(tag)
        silent_struct.fill_pose(p)
        en = silent_struct.energies()[1]
        score_list[tag]=en.value()
        seq_list[tag]=[p.residue(i).name() for i in range(1,p.size()+1)]
        torsion_list[tag]=get_torsions(p)
        hbond_list[tag]=find_hbonds(p)
    return torsion_list,hbond_list,score_list,seq_list

def load_pose(filename):
    p=rosetta.core.import_pose.pose_from_file(filename)
    return p
    
def get_torsions(p):
    torsions=[]
    for i in range(p.total_residue()):
        torsions.append( (p.phi(i+1),p.psi(i+1),p.omega(i+1)))
    return torsions

def score_pose(p):
    scorefxn_tal_name="beta_nov15"
    scorefxn = rosetta.core.scoring.ScoreFunctionFactory.create_score_function(scorefxn_tal_name)
    pc=PeptideCyclizeMover()
    # help (pc)
    pc.apply(p)
    energymethodoptions=scorefxn.energy_method_options()
    energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(True)
    scorefxn.set_energy_method_options( energymethodoptions );
    return scorefxn(p)


def find_hbonds(p):
    x=score_pose(p)
    #p.dump_pdb('temp.pdb')
    p.update_residue_neighbors();
    hbond_set = rosetta.core.scoring.hbonds.HBondSet()
    rosetta.core.scoring.hbonds.fill_hbond_set(p, False, hbond_set)
    hbond_set.setup_for_residue_pair_energies(p , 
                                                      False, 
                                                    False);
    result=[]
    for i_hbond in xrange(hbond_set.nhbonds() ):
            test_hbond=hbond_set.hbond(i_hbond+1)
            if (test_hbond.acc_atm_is_backbone() and test_hbond.don_hatm_is_backbone()):
              #  print test_hbond
             #   print p.residue(1)
                accRes=test_hbond.acc_res()
                donRes=test_hbond.don_res()

                atomAname=p.residue(test_hbond.acc_res()).atom_name(test_hbond.acc_atm())
                atomBname=p.residue(test_hbond.don_res()).atom_name(test_hbond.don_hatm())
                result.append([int(accRes), 
                              int(donRes),
                              str(atomAname),
                              str(atomBname)])
    return result
    
    # def __init__(self):
    #     if not self.b_init:
    #         self.init('-beta_nov15')
    #         self.b_init=True
    #         scorefxn_tal_name="talaris2013"
    #         self.scorefxn_tal2013_vanilla = self.pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function(scorefxn_tal_name)
    #         energymethodoptions=self.scorefxn_tal2013_vanilla.energy_method_options()
    #         energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(True)
    #         self.scorefxn_tal2013_vanilla.set_energy_method_options( energymethodoptions );

    #     print 'Rosetta Initialized'
        
        
