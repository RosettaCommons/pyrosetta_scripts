#!/usr/bin/python

import argparse
from pyrosetta import *
from pyrosetta.rosetta import *

from rosetta.core.pose import *
from rosetta.core.import_pose import *
from rosetta.protocols.rna.denovo import *

"""
Run rna_denovo through a PyRosetta interface. Embed this in a larger python script! Go to town.
"""

def get_argparse():
	parser = argparse.ArgumentParser(description='Run denovo through a PyRosetta interface. Embed this in a larger python script! Go to town.')
	parser.add_argument('--native_pdb', dest='native_pdb', default='', help='Native PDB structure')
	parser.add_argument('--starting_pdbs', dest='starting_pdbs', default='', nargs='+', help='starting PDB(s)')
	parser.add_argument('--align_pdb', dest='align_pdb', default='', help='PDB to which to align')
	parser.add_argument('--nstruct', dest='nstruct', default=1, type=int, help='Number of structures to produce')
	parser.add_argument('--silent_file', dest='silent_file', default='default.out', help='Silent file to fill')
	parser.add_argument('--fasta', dest='fasta', required=True, help='fasta file')
	return parser

def denovo(starting_pdbs, native_pdb, align_pdb, fasta, nstruct, silent_file):
	chm = rosetta.core.chemical.ChemicalManager.get_instance()
	rts = chm.residue_type_set('fa_standard')
	sfxn = core.scoring.ScoreFunctionFactory.create_score_function("stepwise/rna/rna_res_level_energy4.wts")

	# Obtain 'native' and 'align' pose. The 'align' pose is sometimes a subset of the native
	# pose if the starting configuration has multiple freely moving parts.
	if align_pdb == '': align_pdb = native_pdb
	native_pose = core.import_pose.get_pdb_with_full_model_info(native_pdb, rts)

	# Import the starting pose. We don't need to grab the full model info here yet -- we will
	# actually be constructing a more complicated full model info for this Pose because it will
	# be changing over the course of the simulation
	input_pdbs = rosetta.utility.vector1_std_string()
	fasta_files = rosetta.utility.vector1_std_string()
	for pdb in starting_pdbs: input_pdbs.append(pdb)
	# Though it's a vector, you can only give one fasta file. An API oddity!
	fasta_files.append(fasta)
	rna_de_novo_setup = RNA_DeNovoSetup()
	
	rna_de_novo_setup.set_fasta_files(fasta_files) 
	rna_de_novo_setup.set_minimize_rna(True)
	rna_de_novo_setup.set_input_pdbs(input_pdbs) 
	rna_de_novo_setup.set_native_pose(native_pose) 
	rna_de_novo_setup.set_align_pdb(align_pdb) 
	rna_de_novo_setup.set_nstruct(nstruct)
	rna_de_novo_setup.set_silent_file(silent_file)
	
	rna_de_novo_setup.initialize_from_command_line()
	
	pose = rna_de_novo_setup.pose()
	rna_de_novo_protocol = RNA_DeNovoProtocol(rna_de_novo_setup.options(), rna_de_novo_setup.rna_params())
	rna_de_novo_protocol.set_native_pose(rna_de_novo_setup.native_pose())
	rna_de_novo_protocol.apply(pose)

if __name__ == '__main__':
	parser = get_argparse()
	args = parser.parse_args()
	init()
	denovo(args.starting_pdbs, args.native_pdb, args.align_pdb, args.fasta, args.nstruct, args.silent_file)

