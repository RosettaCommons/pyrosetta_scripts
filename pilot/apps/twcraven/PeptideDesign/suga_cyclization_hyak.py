import random
import itertools
import sys
import math
import numpy as np

import pyrosetta
pyrosetta.init('-database /suppscr/baker/twcraven/database -ex1 -ex2 -ex3 -ex4 -score:symmetric_gly_tables true -set_weights atom_pair_constraint 1 dihedral_constraint 1 angle_constraint 1 mm_twist 0 -out:file:output_virtual false -mute all')#

class cyc_pep_predict_suga:
  def __init__(self, res_length, n_struct):
    self.chm = pyrosetta.rosetta.core.chemical.ChemicalManager.get_instance()
    self.rts = self.chm.residue_type_set( 'fa_standard' )
    self.res_length = res_length
    self.build_pose()
    self.declare_bond()
    self.create_constraint()
    self.sfxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function('mm_std')
    self.n_struct = n_struct

  def randomize_all_phipsi(self):
    for i in range(1, len(self.pose)):
      sm = pyrosetta.rosetta.core.scoring.ScoringManager.get_instance()
      rama = sm.get_RamaPrePro()
      torsions = pyrosetta.rosetta.utility.vector1_double()
      rama.random_mainchain_torsions(self.pose.conformation(), self.pose.residue(i).type(), self.pose.residue(i+1).type(), torsions)

      for j in range(1, len(self.pose.residue(i).mainchain_torsions())):
        self.pose.set_torsion( pyrosetta.rosetta.core.id.TorsionID(i, pyrosetta.rosetta.core.id.BB, j), torsions[j] )


    #Need to randomize the 'phi' and 'psi' of the last residue (CYY in this case) or else you will introduce artifacts.
    self.pose.conformation().set_torsion_angle(pyrosetta.rosetta.core.id.AtomID(self.pose.residue(self.res_length).atom_index('C'), self.res_length), pyrosetta.rosetta.core.id.AtomID(self.pose.residue(self.res_length).atom_index('CA'), self.res_length), pyrosetta.rosetta.core.id.AtomID(self.pose.residue(self.res_length).atom_index('N'), self.res_length), pyrosetta.rosetta.core.id.AtomID(self.pose.residue(self.res_length - 1).atom_index('C'), self.res_length - 1), random.uniform(0, 360)*0.0174532925)
    self.pose.conformation().set_torsion_angle(pyrosetta.rosetta.core.id.AtomID(self.pose.residue(self.res_length).atom_index('O'), self.res_length), pyrosetta.rosetta.core.id.AtomID(self.pose.residue(self.res_length).atom_index('C'), self.res_length), pyrosetta.rosetta.core.id.AtomID(self.pose.residue(self.res_length).atom_index('CA'), self.res_length), pyrosetta.rosetta.core.id.AtomID(self.pose.residue(self.res_length).atom_index('N'), self.res_length), random.uniform(0, 360)*0.0174532925)


  def build_pose(self):
    self.pose = pyrosetta.rosetta.core.pose.Pose()
    for i in range(1, self.res_length):
      self.pose.append_residue_by_bond(pyrosetta.rosetta.core.conformation.ResidueFactory.create_residue( self.rts.name_map('GLY')), True)
    self.pose.append_residue_by_bond(pyrosetta.rosetta.core.conformation.ResidueFactory.create_residue( self.rts.name_map('CYY_p:CtermProteinFull')), True)
    for i in range(1, self.res_length):
      self.pose.set_omega(i, 180.0)


  def pack_suga_rotamer(self):
    packer_task = pyrosetta.rosetta.core.pack.task.TaskFactory.create_packer_task(self.pose)
    packer_task.initialize_extra_rotamer_flags_from_command_line()
    packer_task.restrict_to_repacking()
    pack = pyrosetta.rosetta.protocols.simple_moves.PackRotamersMover(self.sfxn, packer_task)
    pack.apply(self.pose)

  def declare_bond(self):
    self.pose.conformation().declare_chemical_bond(1, 'N', self.res_length, 'CE')

  def create_constraint(self):
    dist_cst = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(pyrosetta.rosetta.core.id.AtomID(self.pose[-1].atom_index('V1'), len(self.pose)),
                                                                             pyrosetta.rosetta.core.id.AtomID(self.pose[1].atom_index('N'), 1),
                                                                             pyrosetta.rosetta.core.scoring.func.HarmonicFunc(0.0, 0.05))
    tors_cst = pyrosetta.rosetta.core.scoring.constraints.DihedralConstraint(pyrosetta.rosetta.core.id.AtomID(self.pose[-1].atom_index('CD'), len(self.pose)),
                                                                            pyrosetta.rosetta.core.id.AtomID(self.pose[-1].atom_index('CE'), len(self.pose)),
                                                                            pyrosetta.rosetta.core.id.AtomID(self.pose[1].atom_index('N'), 1),
                                                                            pyrosetta.rosetta.core.id.AtomID(self.pose[1].atom_index('CA'), 1),
                                                                            pyrosetta.rosetta.core.scoring.func.CircularHarmonicFunc(3.14, 0.01))
    ang1_cst = pyrosetta.rosetta.core.scoring.constraints.AngleConstraint(pyrosetta.rosetta.core.id.AtomID(self.pose[-1].atom_index('CD'), len(self.pose)),
                                                                          pyrosetta.rosetta.core.id.AtomID(self.pose[-1].atom_index('CE'), len(self.pose)),
                                                                          pyrosetta.rosetta.core.id.AtomID(self.pose[1].atom_index('N'), 1),
                                                                          pyrosetta.rosetta.core.scoring.func.HarmonicFunc(2.01, 0.01))
    ang2_cst = pyrosetta.rosetta.core.scoring.constraints.AngleConstraint(pyrosetta.rosetta.core.id.AtomID(self.pose[-1].atom_index('CE'), len(self.pose)),
                                                                          pyrosetta.rosetta.core.id.AtomID(self.pose[1].atom_index('N'), 1),
                                                                          pyrosetta.rosetta.core.id.AtomID(self.pose[1].atom_index('CA'), 1),
                                                                          pyrosetta.rosetta.core.scoring.func.HarmonicFunc(2.14, 0.01))
    self.pose.add_constraint(dist_cst)
    self.pose.add_constraint(tors_cst)
    self.pose.add_constraint(ang1_cst)
    self.pose.add_constraint(ang2_cst)

    #These two constraints are unfortunately neccessary to ensure the closed amide gets cleaned up properly.
    dhf = pyrosetta.rosetta.core.scoring.func.HarmonicFunc( 180.0*0.0174533 , 0.01 )
    dhc = pyrosetta.rosetta.core.scoring.constraints.DihedralConstraint( pyrosetta.rosetta.core.id.AtomID(self.pose.residue(1).atom_index('H'), 1) , pyrosetta.rosetta.core.id.AtomID(self.pose.residue(1).atom_index('N'), 1) , pyrosetta.rosetta.core.id.AtomID(self.pose.residue(1).atom_index('CA'), 1) , pyrosetta.rosetta.core.id.AtomID(self.pose.residue(self.res_length).atom_index('CE'), self.res_length) , dhf )
    self.pose.add_constraint( dhc )

    dhf = pyrosetta.rosetta.core.scoring.func.HarmonicFunc( 180.0*0.0174533 , 0.1 )
    dhc = pyrosetta.rosetta.core.scoring.constraints.DihedralConstraint( pyrosetta.rosetta.core.id.AtomID(self.pose.residue(1).atom_index('N'), 1) , pyrosetta.rosetta.core.id.AtomID(self.pose.residue(self.res_length).atom_index('CD'), self.res_length) , pyrosetta.rosetta.core.id.AtomID(self.pose.residue(self.res_length).atom_index('CE'), self.res_length) , pyrosetta.rosetta.core.id.AtomID(self.pose.residue(self.res_length).atom_index('OZ'), self.res_length) , dhf )
    self.pose.add_constraint( dhc )

  def relax(self):
    frm = pyrosetta.rosetta.protocols.relax.FastRelax(pyrosetta.get_score_function())
    frm.ramp_down_constraints(False)
    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_bb_true_range(1, len(self.pose))
    mm.set_chi_true_range(1, len(self.pose))
    frm.set_movemap(mm)
    frm.apply(self.pose)
    declare_bond(self.pose)


  def relax_full(self):
    frm = pyrosetta.rosetta.protocols.relax.FastRelax(pyrosetta.get_score_function())
    frm.ramp_down_constraints(False)
    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_bb(True)
    mm.set_chi(True)
    frm.set_movemap(mm)
    frm.apply(self.pose)
    self.declare_bond()


  def genkic_backbone(self):

    #Need to randomize N-terminal H-N-CA-C torsion. If you don't you will introduce artifacts.
    self.pose.conformation().set_torsion_angle(pyrosetta.rosetta.core.id.AtomID(self.pose.residue(1).atom_index('H'), 1), pyrosetta.rosetta.core.id.AtomID(self.pose.residue(1).atom_index('N'), 1), pyrosetta.rosetta.core.id.AtomID(self.pose.residue(1).atom_index('CA'), 1), pyrosetta.rosetta.core.id.AtomID(self.pose.residue(1).atom_index('C'), 1), random.uniform(0, 360)*0.0174532925)

    genkic = pyrosetta.rosetta.protocols.generalized_kinematic_closure.GeneralizedKIC()

    genkic.add_perturber('randomize_backbone_by_rama_prepro')
    
    cyclization_point_start = res_length
    cyclization_point_end = 1
    anchor_res_min = 2
    anchor_res_max = res_length - 1

    anchor_res = np.random.randint(anchor_res_min, high= anchor_res_max, size = 1)[0]
    first_loop_res =  anchor_res + 1
    last_loop_res = anchor_res - 1 

    middle_loop_res = np.random.randint(cyclization_point_end, high=cyclization_point_start - 3, size = 1)[0]

    if middle_loop_res == last_loop_res:
      middle_loop_res += 3
    elif middle_loop_res == anchor_res:
      middle_loop_res +=2
    elif middle_loop_res == first_loop_res:
      middle_loop_res +=1
    if middle_loop_res > res_length :
      middle_loop_res = middle_loop_res - res_length 

    for i in range(first_loop_res, cyclization_point_start + 1):
      genkic.add_loop_residue(i)
    for i in range (cyclization_point_end, last_loop_res + 1):
      genkic.add_loop_residue(i)

    i = anchor_res + 1
    while (i!=anchor_res):
      if i > cyclization_point_start:
        i = cyclization_point_end
      elif i == anchor_res:
        i = i
      else:
        if i != self.res_length:
          genkic.add_residue_to_perturber_residue_list(1, i)
      i += 1

    genkic.add_filter('loop_bump_check')

    genkic.close_bond(len(self.pose), 'CE', 1, 'N', len(self.pose), 'CD', 1, 'CA', 1.328, 116.206, 121.702, 180, False, False)

    genkic.set_pivot_atoms( first_loop_res, 'CA', middle_loop_res, 'CA', last_loop_res, 'CA' );
    genkic.set_closure_attempts(2000)
    genkic.set_selector_type('lowest_energy_selector')
    genkic.apply(self.pose)

    #Clean up N-terminal H-N-CA-CE and OZ-CE-CD-SG torsions and bond angles as these gets screwed up with KIC for some reason. This is definitely not the most elegant way to do this, but its a quick and dirty way for the moment.

    self.pose.conformation().set_bond_angle( pyrosetta.rosetta.core.id.AtomID( self.pose.residue(1).atom_index('H') , 1 ) , pyrosetta.rosetta.core.id.AtomID( self.pose.residue(1).atom_index('N') , 1 ) , pyrosetta.rosetta.core.id.AtomID( self.pose.residue(1).atom_index('CA') , 1 ), 119.122*0.0174532925 )

    min_e = 10**10; best_HN = 0.0
    for i in range(0, 360):
      self.pose.conformation().set_torsion_angle(pyrosetta.rosetta.core.id.AtomID(self.pose.residue(1).atom_index('H'), 1), pyrosetta.rosetta.core.id.AtomID(self.pose.residue(1).atom_index('N'), 1), pyrosetta.rosetta.core.id.AtomID(self.pose.residue(1).atom_index('CA'), 1), pyrosetta.rosetta.core.id.AtomID(self.pose.residue(1).atom_index('C'), 1), i*0.0174532925)
      current_score = self.sfxn(self.pose)
      if current_score < min_e:
        min_e = current_score
        best_HN = i

    self.pose.conformation().set_torsion_angle(pyrosetta.rosetta.core.id.AtomID(self.pose.residue(1).atom_index('H'), 1), pyrosetta.rosetta.core.id.AtomID(self.pose.residue(1).atom_index('N'), 1), pyrosetta.rosetta.core.id.AtomID(self.pose.residue(1).atom_index('CA'), 1), pyrosetta.rosetta.core.id.AtomID(self.pose.residue(1).atom_index('C'), 1), best_HN*0.0174532925)

    min_e = 10**10; best_O = 0.0
    for i in range(0, 360):
      self.pose.conformation().set_torsion_angle(pyrosetta.rosetta.core.id.AtomID(self.pose.residue(self.res_length).atom_index('OZ'), self.res_length), pyrosetta.rosetta.core.id.AtomID(self.pose.residue(self.res_length).atom_index('CE'), self.res_length), pyrosetta.rosetta.core.id.AtomID(self.pose.residue(self.res_length).atom_index('CD'), self.res_length), pyrosetta.rosetta.core.id.AtomID(self.pose.residue(self.res_length).atom_index('SG'), self.res_length), i*0.0174532925)
      current_score = self.sfxn(self.pose)
      if current_score < min_e:
        min_e = current_score
        best_O = i

    self.pose.conformation().set_torsion_angle(pyrosetta.rosetta.core.id.AtomID(self.pose.residue(self.res_length).atom_index('OZ'), self.res_length), pyrosetta.rosetta.core.id.AtomID(self.pose.residue(self.res_length).atom_index('CE'), self.res_length), pyrosetta.rosetta.core.id.AtomID(self.pose.residue(self.res_length).atom_index('CD'), self.res_length), pyrosetta.rosetta.core.id.AtomID(self.pose.residue(self.res_length).atom_index('SG'), self.res_length), best_O*0.0174532925)


  def closure_check(self):
    #Just check that the N-terminal 'N' and CYY side-chain 'CE' are within 2.0 Angstoms.
    N_coords = self.pose.residue(1).xyz(self.pose.residue(1).atom_index('N'))
    C_coords = self.pose.residue(self.res_length).xyz(self.pose.residue(self.res_length).atom_index('CE'))

    if self.dist(N_coords, C_coords) < 2.0:
        return True
    else:
        return False

  def dist(self, x, y):
    return math.sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2] - y[2])**2)

  def output_pose(self):
    if self.closure_check():
      final_score = self.sfxn(self.pose)
      self.pose.center()
      self.pose.dump_pdb('/suppscr/baker/twcraven/suga_cycles/' + str(self.res_length) + 'mers/' + 'suga_cycle_' + str(self.res_length) + 'mer_' + str(final_score) + '.pdb')
  def run(self):
    for i in range(n_struct):
      #print 'Cyclizing peptide number ', i
      self.randomize_all_phipsi()
      self.pack_suga_rotamer()
      self.genkic_backbone()
      #self.relax_full()
      self.output_pose()


n_struct = int(sys.argv[2])
res_length = int(sys.argv[1])
cyclic_peptide = cyc_pep_predict_suga(res_length, n_struct)
cyclic_peptide.run()
  
