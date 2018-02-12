import glob, sys, math
import operator
import numpy as np
import os
from numpy import *

from rosetta import *
from pyrosetta import *

init('-override_rsd_type_limit -ex1 -ex2 -use_input_sc -ignore_unrecognized_res -mute all -bbdep_omega false -corrections::beta_nov16')

class macrocycle_designifier:
	def __init__(self, file):
		print file
		self.chm = rosetta.core.chemical.ChemicalManager.get_instance()
		self.rts = self.chm.residue_type_set( 'fa_standard' )

		self.sfxn_beta = create_score_function('beta_nov16')
		self.sfxn_mm = create_score_function('mm_std')

		self.rot_set_aro = rosetta.protocols.toolbox.task_operations.LimitAromaChi2_RotamerSetOperation(110.0, 70.0)
		self.rot_set_aro.include_trp(True)



		self.InterfaceAnalyzer = protocols.analysis.InterfaceAnalyzerMover(1, True, self.sfxn_beta, False, False, False, False)
		#help(self.InterfaceAnalyzer)
		
		#Pose should be set up so that A is macrocycle with stub and B is protein target. AB denotes the complex of chains A and B.
		self.AB = pose_from_pdb(file)
		self.AB_clone = self.AB.clone()
		self.split_AB = self.AB.split_by_chain()

		#Determine baseline dg_cross and sc values. These values will be used to determine which derivates are the best. 

		self.InterfaceAnalyzer.apply_const(self.AB)
		#self.starting_sc = self.InterfaceAnalyzer.get_sc###
		self.starting_sc_dg_cross = self.InterfaceAnalyzer.get_crossterm_interface_energy()###
		self.starting_score = self.sfxn_beta(self.AB)
		######Optimize the Stub#######

		#Find the stub residue number. Typically backbone libraries are scanned using the matchifier and are output with either GLY or ALA (w/wo N_methylation) and a single stub.
		#for i in range(1, len(self.split_AB[1])+1):
		#	if self.split_AB[1].residue(i).name().count('CYY') != 0:
		self.stub_residue = 5


		#Minmover for just side-chains.
		self.movemap_sc = MoveMap()
		self.movemap_sc.set_bb(False)
		self.movemap_sc.set_chi(True)
		self.movemap_sc.set_chi(self.stub_residue, False)
		self.movemap_sc.set_jump(False)
		self.minmover_sc_beta = pyrosetta.rosetta.protocols.simple_moves.MinMover(self.movemap_sc, self.sfxn_beta, "linmin", 0.001, True)
		self.minmover_sc_mm = pyrosetta.rosetta.protocols.simple_moves.MinMover(self.movemap_sc, self.sfxn_mm, "linmin", 0.001, True)

		#Minmover for both side-chains and jump 1 (across interface).
		self.movemap_scj = MoveMap()
		self.movemap_scj.set_bb(False)
		self.movemap_scj.set_chi(True)
		self.movemap_scj.set_chi(self.stub_residue, False)
		self.movemap_scj.set_jump(True)
		self.minmover_scj_beta = pyrosetta.rosetta.protocols.simple_moves.MinMover(self.movemap_scj, self.sfxn_beta, "linmin", 0.001, True)
		self.minmover_scj_mm = pyrosetta.rosetta.protocols.simple_moves.MinMover(self.movemap_scj, self.sfxn_mm, "linmin", 0.001, True)

		#Explode examplars into full set of derivatives
		residue_list = ['DPRO', 'ALA', 'PRO', 'GLY', 'MET', 'TYR', 'PHE', 'SER', 'THR', 'HIS', 'ASN', 'GLN', 'ASP', 'GLU', 'VAL', 'LEU', 'ILE', 'ARG', 'LYS']
		#residue_list = ['DPRO', 'ALA', 'GLY', 'PRO']
		self.residues_at_each_position = []
		for i in range(2, len(self.split_AB[1])):
			self.residues_at_each_position.append([])
			for mutation in residue_list:
				self.AB_clone = self.AB.clone()
				#print i, mutation
				self.AB_clone.replace_residue(i, rosetta.core.conformation.ResidueFactory.create_residue( self.rts.name_map(mutation)), True)
				self.pack_clone(i)
				#self.minmover_sc_beta.apply(self.AB_clone)
				self.InterfaceAnalyzer.apply_const(self.AB_clone)
				energies = self.AB_clone.energies()
				weights = [core.scoring.ScoreType(s) for s in range(1, int(core.scoring.end_of_score_type_enumeration) + 1) if self.sfxn_beta.weights()[core.scoring.ScoreType(s)]]
				residue_weighted_energies_matrix = [[energies.residue_total_energies(ii)[w] * self.sfxn_beta.weights()[w] for ii in range(1, self.AB_clone.total_residue() + 1)] for w in weights]
				#new_sc = self.InterfaceAnalyzer.get_sc - self.starting_sc
				new_dg_cross = self.InterfaceAnalyzer.get_crossterm_interface_energy() -  self.starting_sc_dg_cross
				#print len(residue_weighted_energies_matrix)
				#print len(residue_weighted_energies_matrix[i])
				new_rama_prepro = residue_weighted_energies_matrix[-3][i-1]
				new_total_score = self.sfxn_beta(self.AB_clone) - self.starting_score

				self.residues_at_each_position[i-2].append([mutation, new_dg_cross, new_rama_prepro, new_total_score, 10*new_dg_cross + new_rama_prepro + new_total_score])


			self.residues_at_each_position[i-2].sort(key=operator.itemgetter(-1), reverse=False)

			#print self.residues_at_each_position[i-2]

		for i in range(2, len(self.split_AB[1])):
			self.AB.replace_residue(i, rosetta.core.conformation.ResidueFactory.create_residue( self.rts.name_map(self.residues_at_each_position[i-2][0][0])), True)
		self.pack()
		for i in range(0,100):
			self.minmover_sc_beta.apply(self.AB)
		for i in range(0,100):
			self.minmover_scj_beta.apply(self.AB)
		self.InterfaceAnalyzer.apply_const(self.AB)
		self.InterfaceAnalyzer.add_score_info_to_pose(self.AB)

		xx = file.split('_')

		file = file.replace(xx[-1], str(self.sfxn_beta(self.AB)) + '.pdb')

		self.AB.conformation().declare_chemical_bond(1, 'N', self.stub_residue, 'CE')

		self.AB.dump_pdb('../suga_SKP1_complexes_designed/' + file)
		#sys.exit(1)

	def pack_clone(self, cur_res):

		res_length = self.stub_residue
		interface_residues = [2118-2095+res_length+1, 2121-2095+res_length+1, 2128-2095+res_length+1, 2132-2095+res_length+1, 2145-2095+res_length+1, 2154-2095+res_length+1, 2158-2095+res_length+1, 2146-2095+res_length+1]

		pack_l = []
		for l in range(1, len(self.AB_clone.sequence())):
			if (l in interface_residues and l != 1) or l == cur_res:
				pack_l.append(True)
			else:
				pack_l.append(False)

		pack_l = Vector1(pack_l)
		

		pack_task = rosetta.core.pack.task.TaskFactory.create_packer_task(self.AB_clone)
		pack_task.restrict_to_residues(pack_l)
		pack_task.initialize_extra_rotamer_flags_from_command_line()
		pack_task.restrict_to_repacking()
		pack_task.append_rotamerset_operation(self.rot_set_aro)
		pack = rosetta.protocols.simple_moves.PackRotamersMover(self.sfxn_beta, pack_task)
		pack.apply(self.AB_clone)

	def pack(self):
		res_length = self.stub_residue
		pack_l = []
		for l in range(1, len(self.AB.sequence())):
			if l != res_length and l != 1:
				pack_l.append(True)
			else:
				pack_l.append(False)

		pack_l = Vector1(pack_l)
		

		pack_task = rosetta.core.pack.task.TaskFactory.create_packer_task(self.AB)
		pack_task.restrict_to_residues(pack_l)
		pack_task.initialize_extra_rotamer_flags_from_command_line()
		pack_task.restrict_to_repacking()
		pack_task.append_rotamerset_operation(self.rot_set_aro)
		pack = rosetta.protocols.simple_moves.PackRotamersMover(self.sfxn_beta, pack_task)
		pack.apply(self.AB)



split_files = []
f = open(sys.argv[1], 'r')
for line in f:
	split_files.append(line.strip())
f.close()

os.chdir('suga_SKP1_complexes')

for ele in split_files:
	macrocycle_designifier(ele)