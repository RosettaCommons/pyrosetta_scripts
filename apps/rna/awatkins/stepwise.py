#!/usr/bin/python

import argparse
from pyrosetta import *
from pyrosetta.rosetta import *

from rosetta.protocols.stepwise.monte_carlo import *
from rosetta.protocols.stepwise.monte_carlo.options import *
from rosetta.protocols.stepwise.monte_carlo.submotif import *


"""
Run stepwise through a PyRosetta interface. Embed this in a larger python script! Go to town.
"""

def get_argparse():
	parser = argparse.ArgumentParser(description='Run stepwise through a PyRosetta interface. Embed this in a larger python script! Go to town.')
	parser.add_argument('--native_pdb', dest='native_pdb', default='', help='Native PDB structure')
	parser.add_argument('--starting_pdbs', dest='starting_pdbs', default='', nargs='+', help='starting PDB(s)')
	parser.add_argument('--align_pdb', dest='align_pdb', default='', help='PDB to which to align')
	parser.add_argument('--nstruct', dest='nstruct', default=1, type=int, help='Number of structures to produce')
	parser.add_argument('--silent_file', dest='silent_file', default='default.out', help='Silent file to fill')
	parser.add_argument('--fasta', dest='fasta', required=True, help='fasta file')
	return parser

def stepwise(starting_pdbs, native_pdb, align_pdb, fasta, nstruct, silent_file):
	chm = rosetta.core.chemical.ChemicalManager.get_instance()
	rts = chm.residue_type_set('fa_standard')
	sfxn = core.scoring.ScoreFunctionFactory.create_score_function("stepwise/rna/rna_res_level_energy4.wts")

	# Obtain 'native' and 'align' pose. The 'align' pose is sometimes a subset of the native
	# pose if the starting configuration has multiple freely moving parts.
	if align_pdb == '': align_pdb = native_pdb
	native_pose = core.import_pose.get_pdb_with_full_model_info(native_pdb, rts)
	align_pose = core.import_pose.get_pdb_with_full_model_info(align_pdb, rts)

	input_poses = rosetta.utility.vector1_std_shared_ptr_core_pose_Pose_t()
	# Import the starting pose. We don't need to grab the full model info here yet -- we will
	# actually be constructing a more complicated full model info for this Pose because it will
	# be changing over the course of the simulation
	for pdb in starting_pdbs:
		input_poses.append(pose_from_pdb(pdb))

	# Load up the object that constructs the full model info. We need to provide it the full
	# ('target') sequence as well as, potentially, any other input poses that will be combined
	# with the above pose (in this case, there aren't any). Then we read a bunch of default
	# options from the command line defaults. We could go on to configure the builder further
	# with new values for the options we care most about, or we could provide it an 
	# OptionsCollection, if we so chose.
	builder = rosetta.core.import_pose.FullModelPoseBuilder()
	builder.set_fasta_file(fasta)
	builder.set_input_poses(input_poses)
	builder.initialize_further_from_options()
	pose = builder.build()

	# If there is a van der Waals screen operative (read the command line; in
	# our case, there isn't) add it to the Pose. This helps cheaply avoid clashes
	# with a Pose that is too expensive to directly hold in memory.
	protocols.scoring.fill_vdw_cached_rep_screen_info_from_command_line(pose)
	core.import_pose.initialize_native_and_align_pose(native_pose, align_pose, rts, pose)

	# You may want to activate a PyMOL observer here.

	stepwise_monte_carlo = StepWiseMonteCarlo(sfxn)
	stepwise_monte_carlo.set_native_pose(native_pose)
	#stepwise_monte_carlo.set_move( test_move )
	
	options = StepWiseMonteCarloOptions()
	options.initialize_from_command_line()
	stepwise_monte_carlo.set_options(options)
	if (options.from_scratch_frequency() > 0.0 or const_full_model_info(pose).other_pose_list().size() > 0) and not sfxn.has_nonzero_weight(rosetta.core.scoring.other_pose):
		sfxn.set_weight(rosetta.core.scoring.other_pose, 1.0)
	
	stepwise_monte_carlo.set_submotif_library(SubMotifLibrary(rts, options.lores(), options.use_first_jump_for_submotif(), options.exclude_submotifs()))

	# If we do this, then when we call initialize it tries to call the version from RNA_JobDistributor that is pure
	# virtual. I'm not sure how we could let this app run either CSA or MonteCarlo...
	# (jd3 save us)
	stepwise_job_distributor = protocols.rna.setup.RNA_MonteCarloJobDistributor(stepwise_monte_carlo, silent_file, nstruct)

	stepwise_job_distributor.set_native_pose(native_pose);

	# We are doing an apical loop
	stepwise_job_distributor.set_superimpose_over_all(False)
	stepwise_job_distributor.initialize(pose);
	while stepwise_job_distributor.has_another_job():
		stepwise_job_distributor.apply(pose)

if __name__ == '__main__':
	parser = get_argparse()
	args = parser.parse_args()
	init()
	stepwise(args.starting_pdbs, args.native_pdb, args.align_pdb, args.fasta, args.nstruct, args.silent_file)
