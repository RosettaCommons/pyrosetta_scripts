#@file: flex_ddG_monomer.py
#@author: Rebecca Alford (ralford3#@jhu.edu)
#@brief: Calculate the ddG of mutation within a subunitusing backrub sampling

from pyrosetta import * 
from rosetta.protocols.simple_moves import * 
from rosetta.protocols.rosetta_scripts import *
from rosetta.core.scoring import * 

import sys, os
import math
import numpy as np
import sqlite3

from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )


def main( args ): 

	parser = OptionParser( usage="usage: %prog --pdb sample.pdb --resfile resfile.txt" )
	parser.set_description(main.__doc__)

	parser.add_option('--pdb', '-p', 
		action="store",
		help="Path of PDB file")

	parser.add_option('--resfile', '-r', 
		action="store",
		help="Path of mutation resfile")

	(options, args) = parser.parse_args(args=args[1:])
	global Options
	Options = options 

	if ( not Options.pdb or not Options.resfile ): 
		sys.exit( "Must provide flags --pdb and --resfile! Exiting..." )

	# Initialize PyRosetta from options
	init( extra_options="-in:ignore_unrecognized_res -restore_talaris_behavior")

	# Initialize pose from PDB file
	pose = pose_from_file( Options.pdb )

	# Initialize movers from Rosetta Scripts XML code
	objs = initialize_from_xml( Options.resfile )

	# Run monomer flex ddG protocol
	run_monomer_flex_ddG( pose, objs )

#################################################################################
# @brief Sample the native and mutated conformations using backrub + minimization with atom 
# pair constraints. Return copies of the native, minimized pose and the mutated + backrub pose
def run_monomer_flex_ddG( pose, objs ):

	### Get both fa score functions (w and w/o csts) 
	fa_sfxn = objs.get_score_function( "fa_talaris2014" )
	fa_sfxn_cst = objs.get_score_function( "fa_talaris2014_cst" )

	### Minimize the pose with constraints wrt the mutation site
	objs.get_mover( "addcst" ).apply( pose )
	fa_sfxn_cst.score( pose )
	objs.get_mover( "neighbor_shell_storer" ).apply( pose )
	objs.get_mover( "minimize" ).apply( pose )
	objs.get_mover( "clearcst" ).apply( pose )

	# Apply backrub protocol
	objs.get_mover( "backrub" ).apply( pose )

#################################################################################
# @brief Initialize Rosetta objects from the flex_ddG XML script provided in Barlow 
# et al. BioRxiv 2017 including scorefxns, task operations, residue selectors and movers
def initialize_from_xml( resfile ): 

	padded_resfile = '"' + resfile + '"'
	xmlstr = """
		<SCOREFXNS>
			<ScoreFunction name="fa_talaris2014" weights="talaris2014"/>
			<ScoreFunction name="fa_talaris2014_cst" weights="talaris2014">
				<Reweight scoretype="atom_pair_constraint" weight="1.0"/>
				<Set fa_max_dis="9.0"/>
			</ScoreFunction>
		</SCOREFXNS>		
		<TASKOPERATIONS> 
			<ReadResfile name="res_mutate" filename=""" + padded_resfile + """/>
		</TASKOPERATIONS>
		<RESIDUE_SELECTORS> 
			<Task name="resselector" fixed="0" packable="0" designable="1" task_operations="res_mutate"/>
			<Neighborhood name="bubble" selector="resselector" distance="8.0"/>
			<PrimarySequenceNeighborhood name="bubble_adjacent" selector="bubble" lower="1" upper="1" />
			<StoredResidueSubset name="restore_neighbor_shell" subset_name="neighbor_shell"/>
			<Not name="everythingelse" selector="restore_neighbor_shell"/>
		</RESIDUE_SELECTORS>
		<TASKOPERATIONS>
			<OperateOnResidueSubset name="repackonly" selector="restore_neighbor_shell">
				<RestrictToRepackingRLT/>
			</OperateOnResidueSubset>
			<OperateOnResidueSubset name="norepack" selector="everythingelse">
				<PreventRepackingRLT/>
			</OperateOnResidueSubset>
			<UseMultiCoolAnnealer name="multicool" states="6"/>
			<ExtraChiCutoff name="extrachizero" extrachi_cutoff="0"/>
			<InitializeFromCommandline name="commandline_init"/>
			<RestrictToRepacking name="restrict_to_repacking"/>
		</TASKOPERATIONS>
		<MOVERS>
			<StoreResidueSubset name="neighbor_shell_storer" subset_name="neighbor_shell" 
			residue_selector="bubble_adjacent"/>
			
			<AddConstraintsToCurrentConformationMover name="addcst" use_distance_cst="1" 
			coord_dev="0.5" min_seq_sep="0" max_distance="9" CA_only="1" bound_width="0.0" 
			cst_weight="0.0"/>
			<ClearConstraintsMover name="clearcst"/>

			<MinMover name="minimize" scorefxn="fa_talaris2014_cst" chi="1" bb="1" 
			type="lbfgs_armijo_nonmonotone" tolerance="0.000001" max_iter="5000" 
			abs_score_convergence_threshold="1.0"/>
			
			<PackRotamersMover name="repack" scorefxn="fa_talaris2014" 
			task_operations="commandline_init,repackonly,norepack,multicool"/>
			<PackRotamersMover name="mutate" scorefxn="fa_talaris2014" 
			task_operations="commandline_init,res_mutate,norepack,multicool"/>

			<ReportToDB name="dbreport" batch_description="interface_ddG" database_name="ddG2.db3">
				<ScoreTypeFeatures/>
				<ScoreFunctionFeatures scorefxn="fa_talaris2014"/>
				<StructureScoresFeatures scorefxn="fa_talaris2014"/>
			</ReportToDB>

			<ReportToDB name="structreport" batch_description="interface_ddG_struct" database_name="struct.db3">
				<PoseConformationFeatures/>
				<PdbDataFeatures/>
				<JobDataFeatures/>
				<ResidueFeatures/>
				<PoseCommentsFeatures/>
				<ProteinResidueConformationFeatures/>
				<ResidueConformationFeatures/>
			</ReportToDB>

			<SavePoseMover name="save_wt_bound_pose" restore_pose="0" reference_name="wt_bound_pose"/>
			<SavePoseMover name="save_backrub_pose" restore_pose="0" reference_name="backrubpdb"/>
			<SavePoseMover name="restore_backrub_pose" restore_pose="1" reference_name="backrubpdb"/>

			<ParsedProtocol name="finish_ddg_post_backrub">
				<Add mover_name="save_backrub_pose"/>
				<Add mover_name="dbreport"/>
				
				<Add mover_name="repack"/>
				
				<Add mover_name="addcst"/>
				<Add mover_name="minimize"/>
				<Add mover_name="clearcst"/>
				
				<Add mover_name="save_wt_bound_pose"/>
				<Add mover_name="dbreport"/>

				<Add mover_name="restore_backrub_pose"/>

				<Add mover_name="mutate"/>

				<Add mover_name="addcst"/>
				<Add mover_name="minimize"/>
				<Add mover_name="clearcst"/>

				<Add mover_name="dbreport"/>

			</ParsedProtocol>

			<BackrubProtocol name="backrub" mc_kt="1.2" ntrials="35000" pivot_residue_selector="restore_neighbor_shell" 
			task_operations="restrict_to_repacking,commandline_init,extrachizero" recover_low="0" trajectory_stride="2500" 
			trajectory_apply_mover="finish_ddg_post_backrub"/>

		</MOVERS>
		"""

	xmlobjs = XmlObjects.create_from_string( xmlstr )
	return xmlobjs

if __name__ == "__main__" : main(sys.argv)

