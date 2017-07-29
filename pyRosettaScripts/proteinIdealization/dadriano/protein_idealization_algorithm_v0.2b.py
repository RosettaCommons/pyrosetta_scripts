#!/usr/bin/env python

##DON'T Remove these comments
#DASM Protein Idealizer (not an official name)
#This is a novel method that uses Statistical Clustered Information from the vall to idealize SS and build de-novo loops"
#Author: Daniel-Adriano Silva
#Ver: internal release v0.1b
#Date:  Nov/6/2014
#Intended use: Fold-It, general profiling of protein's design quality
#Citation: To the date this is still unpublished work, if you inted to citate contact the authors. 
#Contact: Daniel-Adriano Silva <dadriano@gmail.com> or  David Baker <dabaker@u.washington.edu>
#  "***THIS IS AN INTERNAL DEVELOPER VERSION. DO NOT REDISTRIBUTE***"
#  "DASM Protein Idealizer (not an official name)"
#  "Version: internal release v0.1b, Nov/6/2014"
#  "This is a novel method that uses Statistical Clustered Information from the vall to idealize SS and build de-novo loops"
#  "Authors: Daniel-Adriano Silva, Enrique Marcos and David Baker. 2014. Rosetta Commons & The Baker lab."
#  "Questions: dadriano@gmail.com"
#  "Documentation: \"Does not exist yet\""
#  "Copyright "The Baker Lab & Rosetta Commons (2014)"
#  "***THIS IS AN INTERNAL DEVELOPER VERSION. DO NOT REDISTRIBUTE***"
####

#Import libraries
import copy as copy
import os as os
import protein_idealizer as protein_idealizer
import argparse as argparse

#Main function
def main():
	protein_idealizer.licence()
	
	print "Reading user arguments: "
	#Read User options
	# General
	parser = argparse.ArgumentParser()
	parser.add_argument("-1", "--stage1", default=1, type=int)
	parser.add_argument("-2", "--stage2", default=1, type=int)
	parser.add_argument("-d", "--clustering_database_path", default="./4merDB", type=str)
	parser.add_argument("-l", "--lres_clusteringDB_res_suffix", default="0.5", type=str)
	parser.add_argument("-m", "--hres_clusteringDB_res_suffix", default="0.2", type=str)
	parser.add_argument("-r", "--root_work_dir", default="./", type=str)
	parser.add_argument("-p", "--b_make_plots", default=1, type=int)
	parser.add_argument("-b", "--b_debug", default=0, type=int)
	# Stage1 
	parser.add_argument("-a", "--stage1_input_pdb_name", default="stage1.pdb" , type=str)
	parser.add_argument("-c", "--stage1_output_prefix", default="stage1" , type=str)
	parser.add_argument("-e", "--stage1_context_pdb_name", default="stage1_context.pdb", type=str)
	parser.add_argument("-f", "--stage1_tmp_dir_name", default="stage1_idealization_tmp", type=str)
	parser.add_argument("-g", "--stage1_target_num_loops", default=100, type=int)
	parser.add_argument("-i", "--stage1_sufficient_num_loops", default=100, type=int)
	parser.add_argument("-j", "--stage1_sufficient_num_loops_for_each_len", default=40, type=int)
	# Stage2
	parser.add_argument("-k", "--stage2_input_pdb_name", default="stage1_result_minimized.pdb", type=str)
	parser.add_argument("-n", "--stage2_context_pdb_name", default="stage1_context.pdb", type=str)
	parser.add_argument("-o", "--stage2_output_prefix", default="stage2", type=str)
	parser.add_argument("-q", "--stage2_opt_cycles", default=20, type=str)
	parser.add_argument("-s", "--stage2_tmp_dir_name", default="stage2_seqOptimization_tmp", type=str)
	parser.add_argument("-t", "--b_use_symmetric_packing", default=0, type=int)
	parser.add_argument("-v", "--pseudo_symm_starts", default=-1, type=int)
	parser.add_argument("-w", "--pseudo_symm_ends", default=-1, type=int)
	
	args = parser.parse_args()
	print "Readed arguments: ", args
	argsdict = vars(args)
	
	print "Setting program's global variables"
	b_run_stage1= bool(argsdict['stage1']==True)
	b_run_stage2= bool(argsdict['stage2']==True)
	b_debug= bool(argsdict['b_debug']==True)
	b_make_plots= bool(argsdict['b_make_plots']==True)
	clustering_database_path=str(argsdict['clustering_database_path'])
	lres_clustering=str(argsdict['lres_clusteringDB_res_suffix'])
	hres_clustering=str(argsdict['hres_clusteringDB_res_suffix'])
	general_out_dir=str(argsdict['root_work_dir'])
	
	stage1_input_pdb_name=str(argsdict['stage1_input_pdb_name'])
	stage1_output_prefix=str(argsdict['stage1_output_prefix'])
	stage1_context_pdb_name=str(argsdict['stage1_context_pdb_name'])
	stage1_tmp_dir_name=str(argsdict['stage1_tmp_dir_name'])
	stage1_target_num_loops=int(argsdict['stage1_target_num_loops'])
	stage1_sufficient_num_loops=int(argsdict['stage1_sufficient_num_loops'])
	stage1_sufficient_num_loops_for_each_len=int(argsdict['stage1_sufficient_num_loops_for_each_len'])
	
	stage2_input_pdb_name=str(argsdict['stage2_input_pdb_name'])
	stage2_output_prefix=str(argsdict['stage2_output_prefix'])
	stage2_context_pdb_name=str(argsdict['stage2_context_pdb_name'])
	stage2_tmp_dir_name=str(argsdict['stage2_tmp_dir_name'])
	stage2_opt_cycles=int(argsdict['stage2_opt_cycles'])
	
	b_use_symmetric_packing=bool(argsdict['b_use_symmetric_packing']==True)
	pseudo_symm_starts=int(argsdict['pseudo_symm_starts'])
	pseudo_symm_ends=int(argsdict['pseudo_symm_ends'])

	print "Done setting global variables"
	
	#Set some options if running both stages
	if (b_run_stage1 and b_run_stage2):
		#Set default names for stage2 input files
		stage2_input_pdb_name="stage1_result_minimized.pdb" 
		stage2_context_pdb_name="stage1_context.pdb"
	
	#Initialize the 4mer clusters database and the fragment matcher
	pide_db=protein_idealizer.database()
	pide_db.read_clustering_databases(clustering_database_path=clustering_database_path, 
									  lres=lres_clustering, 
									  hres=hres_clustering, 
									  tcm_jump_step=3)
									  
	#Instantiate our a.a. dictionary and layer definitions
	pide_pLayers=protein_idealizer.proteinLayers()
	##print pide_pLayers.__dict__
	#OVERWRITE THE default DIC to promote large residues in the Core/Interfaces
	#OptionalNow set our personal layer preferences
	print "Now overwriting the default layer dictionary with our personal preference"
	pide_pLayers.CoreAAdic_extended={'A','F','I','L','M','P','V','W','Y','D','N','S','T'}
	pide_pLayers.CoreAAdic={'F','I','L','M','W','P','V'}   # 'A'
	pide_pLayers.BoundaryAAdic={'A','D','E','F','G','I','K','L','M','N','P','Q','R','S','T','V','W','Y'} # 'H'
	pide_pLayers.SurfAAdic={'A','G','D','E','G','K','N','P','Q','R','S','T'} # 'H'
	pide_pLayers.SepecialAAdic={'C','U','H'} #Special (No in use right now except to complete the pool of a.a., put here anything you doon't wanna use )
	#Check layer dictionary for completeness:
	pide_pLayers.check_aa_layer_dictionaries()
	print "Details:", pide_pLayers.__dict__
	
	#Finally Instantiate the protein idealizer
	pide_idealizer=protein_idealizer.idealizer()
	
	
	if (b_run_stage1):
		print "Starting Stage 1"
		
		#Set stage1 tmp output directories. Generate if they don't exist
		out_stage1_tmp_dir="%s/%s"%(general_out_dir, stage1_tmp_dir_name)
		print "Will output tmp files to path: ", out_stage1_tmp_dir
		if not os.path.exists(out_stage1_tmp_dir):
			print "Generating out directory: ", out_stage1_tmp_dir
			os.makedirs(out_stage1_tmp_dir)
			
		out_stage1_tmp_loop_dir="%s/%s/deNovo_possible_loops"%(general_out_dir, stage1_tmp_dir_name)
		print "Will output possible loops built to path: ", out_stage1_tmp_loop_dir
		if not os.path.exists(out_stage1_tmp_loop_dir):
			print "Generating out directory: ", out_stage1_tmp_loop_dir
			os.makedirs(out_stage1_tmp_loop_dir)
			
		#Read Input PDBs
		print "Reading Stage 1 input files"
		[inputPose_original, 
		 inputPose_idealized_bonds,
		 context_pose,
		 context_pose_centroid,
		 b_has_context_pose] = pide_idealizer.read_input_pdbs_stage1(general_input_pdb=stage1_input_pdb_name, 
																		general_context_pdb_name=stage1_context_pdb_name,
																		general_out_dir=general_out_dir,
																		stage1_tmp_dir=out_stage1_tmp_dir)
																		
		#Split the pose into its SS components
		ssCore_pose_array = pide_idealizer.split_pose_byto_SScore(inputPose_original=inputPose_idealized_bonds, 
																	clustersDB=pide_db, 
																	general_out_dir=general_out_dir,
																	stage1_tmp_dir=out_stage1_tmp_dir)
																	
		#Idealize the SScore backbone using fragments
		print "Now I'll reconstruct the SS core using idealized fragments"
		ssCore_idealized_pose_array = pide_idealizer.idealize_ss_array(pose_array=ssCore_pose_array,
				context_pose_centroid=context_pose_centroid,
				b_has_context_pose=b_has_context_pose,
				clustersDB=pide_db,
				ss_fragment_population_threshold=30,  #100 is default "quite flexible maybe", 10 is allow almost all, 200 is for very ideal structs
				fragment_designability_ratio=0.5,     #0.8 is default "be VERY very carefull when modifing this value", < will afect your ability to design
				rmsd_tolerance=0.8,                   #0.8 is default, < is OK (more strict)
				rmsd_sliding_window_size=8,           #8 is default The shorter, the nearer to input, the larger  the more idealized struct will come out
				stage1_tmp_dir=out_stage1_tmp_dir)
		print "Done idealizing the SS core"
		
		#Backup idealized SS array
		print "Generating a memory backup of %d SS elements: "%len(ssCore_idealized_pose_array)
		ssCore_idealized_pose_array_backup=copy.deepcopy(ssCore_idealized_pose_array)
		print "Done"
		
		#Custom/swap SS ordering 
		#ToDo: Implement an automatic min_path/opt_path finder
		#NOTE: disabled for now
		if False:
			b_use_custom_SS_order=False
			custom_SS_order=[0,1,2,7,6,5,4,3,8,9,10]
			
			idealized_pose_array=[]
			if b_use_custom_SS_order:
				print "Using custom SS order"
				assert(len(idealized_pose_array) <= custom_SS_order)
				print custom_SS_order
				for i in range(len(custom_SS_order)):
					idealized_pose_array.append(ssCore_idealized_pose_array_backup[custom_SS_order[i]])
					idealized_pose_array[-1].dump_pdb("%s/test_reordered_idealSS_%02d_%02d.pdb"%(general_out_dir,i, custom_SS_order[i]))
			else:
				print "Restoring copy of #SS: ", len(ssCore_idealized_pose_array_backup), "using original ordering"
				for ideal_ss_pose in ssCore_idealized_pose_array_backup:
					idealized_pose_array.append(ideal_ss_pose.clone())
		#END disabled code
		
		
		#Build ideal de-novo loop connections
		print "I'll build a colletion of ideal de-novo loops that can reconnect the disjointed SS, be patient, this might take a long time!"
		ideal_loops_array = pide_idealizer.build_loops(idealized_pose_array=copy.deepcopy(ssCore_idealized_pose_array),
							context_pose_centroid=context_pose_centroid,
							b_has_context_pose=b_has_context_pose,
							target_num_loops=100, #Stop condition
							sufficient_num_loops=100, #Suffience condition to succed
							sufficient_num_loops_for_each_len=40, #Suffience condition to succed per loop len
							min_loop_size=0, #Min desired loop size
							max_loop_size=10, #Max desired loop size
							max_num_frag_matches=50,  #N-side seeds number of neighbors
							max_resolution_frag_matches=0.6,  #N-side seeds max matching res
							clustersDB=pide_db,
							out_frag_tmp_dir=out_stage1_tmp_loop_dir)

		#Double check loops solution size (i.e. that we have enough loops)
		print "Num of connections with solution: ", (len(ssCore_idealized_pose_array)-1), "of", (len(ideal_loops_array))
		if ( (len(ssCore_idealized_pose_array)-1) == (len(ideal_loops_array))):
			print "PASS"
		else:
			print "FAILED to find some connections. STOP"
			
		print "Finished building ideal de-novo loops"
		
		#Combine SS and loops and find the best solution(s)
		#Allow alternative solutions (MC style maybe?)
		print "Now I'll try to find suitable combinations of the loops to reconnect the disjointed SSs"
		b_success, connected_result_pose, ss_relative_indexes_2_input = pide_idealizer.combine_SS_and_loops(
																			idealized_pose_array=ssCore_idealized_pose_array,
																			loop_poses_array=ideal_loops_array,
																			entropic_acceptance_ratio_score=0.70,
																			general_out_dir=general_out_dir,
																			out_prefix_name=stage1_output_prefix)
		if b_success==True:
			print "Pose reconnected with novo-ideal-loops sucessfully!"
		else:
			print "FAILED to find solutions to reconnect the pose. Maybe try to find more loops can help?"
			assert (0==1)
			
		#Minimize the result, and make some fancy plots
		[connected_result_pose_minimized, 
		 min_dif_score, 
		 max_dif_score] = pide_idealizer.minimize_calculate_difficulty_and_plot(connected_pose=connected_result_pose, 
																					b_make_plots=b_make_plots,
																					clustersDB=pide_db,
																					general_out_dir=general_out_dir,
																					out_prefix_name=stage1_output_prefix)
		#Note: 1/foo_score = folds-of-goodnes!!! The higher the better, 
		#  better than 10000 is usually required for easy designability purposes
		print "Number of fragment/unique-seq observations ration"
		print "Worst fragment stats (obs num):", 1/max_dif_score
		print "Best fragment stats (obs num):", 1/min_dif_score, 
		print "#**Note: 1/foo_score = folds-of-goodnes!!! The higher the better,\n\
		  worst>~10000 is usually required for easy designability purposes"
		
		print "End of Stage 1"
		
		
	if (b_run_stage2):
		print "Starting Stage 2"
		#Set stage 2 tmp output directory. Generate if it doesn't exist
		stage2_tmp_dir="%s/%s"%(general_out_dir,stage2_tmp_dir_name)
		print "Will output Stage2 (sequence optimization) files to dir_name: ", stage2_tmp_dir_name
		if not os.path.exists(stage2_tmp_dir):
			print "Generating out directory: ", stage2_tmp_dir
			os.makedirs(stage2_tmp_dir)
			
		#Read and process stage2 input PDBs
		stage_2_in_pose = pide_idealizer.read_input_pdbs_stage2(stage2_input_pdb_name=stage2_input_pdb_name,
																  general_out_dir=general_out_dir)
																  
		#Setup symmetric packing booleans
		symm_positions_dic=pide_idealizer.setup_pseudo_symmetric_packing(connected_pose_mini_desig=stage_2_in_pose,
															b_use_symmetric_packing=b_use_symmetric_packing,
															pseudo_symm_starts=pseudo_symm_starts, 
															pseudo_symm_ends=pseudo_symm_ends)
															
		#Calculate layers dictionary and rosetta's packer booleans
		print "Calculating layer's booleans and tasks"
		stage_2_in_pose_layer_dic_array, stage_2_in_pose_layer_design_task_bools=pide_idealizer.calculate_pose_layers(connected_pose_mini_desig=stage_2_in_pose,
														layersDic=pide_pLayers,
														testSetAA=['V'],  #array, a.a. to test
														c_max_Cb_val=[198.00], #array, Constant describing the maximmum Cb SASA per a.a.
														core_cut_off=0.03,
														boundary_cut_off=0.1,
														general_out_dir=general_out_dir,
														stage2_output_prefix=stage2_output_prefix)
		print "Done calculating Layer's booleans and tasks"
														
		#Generate some resfiles that can be used for design with rosetta (outside of the idealizer) later
		print "Generating resfiles for your convinince :)"
		optimized_pikaa=pide_idealizer.calculate_possible_aa_by_res_using_fragments(connected_pose_mini_desig=stage_2_in_pose,
																					this_pose_layer_dic_array=stage_2_in_pose_layer_dic_array,
																					clustersDB=pide_db,
																					min_tar_num_seq_match=10,
																					to_center_cut_off_dist=0.7,
																					general_out_dir=general_out_dir,
																					stage2_output_prefix=stage2_output_prefix)
		print "Done generating resfiles"
		
		#Run the automatic optimization cycle
		#ToDo: Still needs a lot of improvement
		[connected_pose_optimized, 
		 seq_connected_pose_optimized,
		 connected_pose_w_optimized_final_opt_score,
		 list_of_optimized_positions_unique]=pide_idealizer.optimize_sequence_using_clustered_fragments_stats(
															in_pose=stage_2_in_pose,
															clustersDB=pide_db,
															this_pose_layer_dic_array=stage_2_in_pose_layer_dic_array,
															layer_design_task_bools=stage_2_in_pose_layer_design_task_bools,
															design_quality_entropy_cut_off=100000.0, #Adaptive initial target quality ratio
															min_target_design_quality_entropy_cut_off=2.0, #Adaptive final target quality ratio
															design_quality_entropy_look_within_percentage=0.1, #in the scale of 1.0 to 0.0
															out_cycles=stage2_opt_cycles, #AutoOptOption
															in_cycles=5, #AutoOptOption
															b_generate_quality_plot=True, #AutoOptOption
															b_reassign_each_out_cycle=True, #AutoOptOption
															b_allow_alternate_good_and_bad_sequence_regions_design_each_in_cycle=True, #AutoOptOption
															b_reassing_and_minimize_before_output_struct=False, #AutoOptOption
															b_allow_design_neighbors=True, #AutoOptOption
															quality_to_start_rotamer_minimization=100.0, #AutoOptOption
															quality_to_start_structure_relaxation=2.0, #AutoOptOption
															out_cycle_max_minimization_rmsd=1.0, #AutoOptOption
															out_cycle_max_relax_rmsd=1.5, #AutoOptOption
															target_max_cluster_distance=0.8,  #0.8 not fair??? #FineTuneOption
															num_cycles_without_improvement_tol=10, #ex. Use (in_cycles*2), #FineTuneOption
															allow_frag_seq_retry=True, #FineTuneOption
															max_allowed_tar_num_seq_match=30, #Adaptive max number of sequences to match for packing
															min_allowed_tar_num_seq_match=3, #Adaptive min number of sequences to match for packing
															tar_num_seq_match=30, #use e.q. to max_allowed_tar_num_seq_match #Adaptivenumber of sequence to match for packing
															max_allowed_tar_seq_match_min_ratio=0.40, #Adaptive ratio for sequences to match for packing #Don't down much this number, it is better to reduce tar_num_seq_match
															min_allowed_tar_seq_match_min_ratio=0.10, #Adaptive ratio for sequences to match for packing #0.30 works best??
															tar_seq_match_min_ratio=0.40 , #use e.q. to max_allowed_tar_seq_match_min_ratio #Adaptive ratio for sequences to match for packing 
															mc_temperature=500.0,  #mcTemperatureOption Start
															max_mc_temperature=500.0,  #mcTemperatureOption Max
															min_allowed_mc_temperature=0, #mcTemperatureOption Min
															general_out_dir=general_out_dir,
															stage2_tmp_dir_name=stage2_tmp_dir_name,
															locked_aa_indexes_identities_array=[], #Opt. ex: locked_aa_indexes_identities_array.append([2, 'E'])
															b_lock_coreSS_aa=False, #Do lock a.a. that belong to SScore during the optimization?
															b_use_symmetric_packing=b_use_symmetric_packing,
															b_Debug=False )

															
		#Output files and end stage2
		tmp_out_final_pdb_filename="%s/%s_OptResult.pdb"%(general_out_dir, stage2_output_prefix)
		print "Generating output of Stage2's (sequence optimized) to: ", tmp_out_final_pdb_filename
		connected_pose_optimized.dump_pdb(tmp_out_final_pdb_filename)

		print "Details of the result are: "
		print " Score: ", connected_pose_w_optimized_final_opt_score
		print " Opt_positions: ", list_of_optimized_positions_unique
		print " Sequence: ", seq_connected_pose_optimized
		print "End of sequence optimization. You might want to run design now using the loop-locked resfile." 
		print "See ya soon. DASM"
		
	print "DASM Protein Idealizer Protocol Finished (not an official name). But... I'll be back :~)"
	print "Authors: Daniel-Adriano Silva (dadriano), Enrique Marcos and David Baker"
	print protein_idealizer.licence()
	#End

#Call to the main function
if __name__ == "__main__":
    main()




