#The Protein Idealizer (a.k.a. The De Novo Protein Mimic Designer) by D.A.S.
#Authors: Daniel-Adriano Silva, Enrique Marcos, Javier Castellanos and David Baker. 
##Contributions:
### D.-A.S.: Original Idea, Algoritmic and Code Implementation.
### E.M.: Advice on priciples for de novo protein design.
### J.C.: Advice on Automation
### D.B.: Advisor, P.I.

#The Baker Lab. UW, Seattle, USA. 2014-2019"
#Date: Dec/16/2018

print "De Novo Protein Mimic Designer, by D.A.S."

print "\n Sure! There are indeed more than a billion ways to implement the ideas here.\nNevertheless, what is really fundamental (and novel) here is the concept and what it proofs (i.e.'same structure equals same function'\nD.A.S."

print "Citation: \n  Daniel-Adriano Silva*, Shawn Yu*, Umut Y. Ulge*, Jamie B. Spangler, et. al., De novo design of potent and selective mimics of IL-2 and IL-15, Nature, 2019. https://doi.org/10.1038/s41586-018-0830-7"

print "HI D. Let's start!"

from protein_mimic_designer_r2_00b_libs import *
import argparse

def parse_arguments():
    print 
    parser = argparse.ArgumentParser()

    parser.add_argument("-in_pdb_name", default=None,
                        type=str, dest="in_pdb_name")

    parser.add_argument("-in_global_context_path", default=None,
                        type=str, dest="in_global_context_path")

    parser.add_argument("-out_file_dir", dest="out_file_dir", 
                        type=str, default="./")

    parser.add_argument("-keep_intact_chains_list", 
                        type=int, nargs='*', dest="keep_intact_chains_list", default=[])

    parser.add_argument("-target_reduc_size_Nterminus", 
                            type=int, dest="target_reduc_size_Nterminus", default=4)

    parser.add_argument("-target_reduc_size_Cterminus", 
                            type=int, dest="target_reduc_size_Cterminus", default=4)

    parser.add_argument("-target_extra_size_Nterminus", 
                            type=int, dest="target_extra_size_Nterminus", default=12)

    parser.add_argument("-target_extra_size_Cterminus", 
                            type=int, dest="target_extra_size_Cterminus", default=12)

    parser.add_argument("-enable_only_profile_mode", 
                        action='store_true', dest='enable_only_profile_mode', default=False, 
                        help='switch to activate profiling protein only')

    parser.add_argument("-max_num_results_per_topology", 
                        type=int, dest='max_num_results_per_topology', default=20)

    parser.add_argument("-keep_only_top_loops_num_max_during_clustering", 
                        type=int, dest='keep_only_top_loops_num_max_during_clustering', default=100)

    parser.add_argument("-max_num_loops_per_building_hop",
                        type=int, dest='max_num_loops_per_building_hop', default=100)

    parser.add_argument('-force_chain_pairs', action='append',nargs='+', type=str)
    
    parser.add_argument('-clustered_fragments_path', dest='clustered_fragments_path',
                        type=str, default='./dbClusteredFragments_7aa/kcenters_stats.dat_res1.36.gz')
    parser.add_argument('-clustered_fragments_assignments_path', dest='clustered_fragments_assignments_path', 
                        type=str, default='./dbClusteredFragments_7aa/kcenters_assignments.dat_res1.36.gz')

    args=vars(parser.parse_args())
    if args['force_chain_pairs']:
        print args['force_chain_pairs']
        args['force_chain_pairs']=np.asarray([x.split(',') for x in args['force_chain_pairs'][0]], int)
    else:
        args['force_chain_pairs']=[]

    print("Gathering command line arguments...:", args)
    return(args)


def main():
    #initialize the dlooper container
    dlooper_alg=None

    #parse arguments
    args=parse_arguments()
    
    #Numpy options
    np.set_printoptions(precision=2)

    ##Init Rosetta, Database and Scoring
    print "Init Rosetta and Reading Rosetta Database, be patient"
    pyrosetta.init(options="-beta  -ex1 -ex2aro -detect_disulf FALSE -mute all" , #  -corrections::beta -constant_seed
                 extra_options="", 
                 set_logging_handler=True) #Not detecting disulfides... it can lead to interesting cases!!!

    #SFX
    scorefxn_vanilla_loops_name="beta_cart"
    scorefxn_vanilla_loops = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function(scorefxn_vanilla_loops_name)

    scorefxn_ramaHyper_name="beta_cart"
    scorefxn_ramaHyper = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function(scorefxn_ramaHyper_name)
    scorefxn_ramaHyper.set_weight(pyrosetta.rosetta.core.scoring.rama , scorefxn_ramaHyper.get_weight(pyrosetta.rosetta.core.scoring.rama)*3.0)
    scorefxn_ramaHyper.set_weight(pyrosetta.rosetta.core.scoring.hbond_bb_sc , scorefxn_ramaHyper.get_weight(pyrosetta.rosetta.core.scoring.hbond_bb_sc)*3.0)
    scorefxn_ramaHyper.set_weight(pyrosetta.rosetta.core.scoring.hbond_lr_bb , scorefxn_ramaHyper.get_weight(pyrosetta.rosetta.core.scoring.hbond_lr_bb)*3.0)
    scorefxn_ramaHyper.set_weight(pyrosetta.rosetta.core.scoring.hbond_sc , scorefxn_ramaHyper.get_weight(pyrosetta.rosetta.core.scoring.hbond_sc)*3.0)
    #print scorefxn_ramaHyper

    print "Done"



    #ToDo: 
    #1. Add limits for the angle parameters for different SS types!

    #Options for users
    in_pdb_name=args['in_pdb_name'] 
    #Global context?
    in_global_context_path=args['in_global_context_path']
    b_use_global_context=False
    if (in_global_context_path):
        print "Activating global context use:", in_global_context_path
        b_use_global_context=True
    else:
        print "Global context is deactivated since there was no specified input for it"
    max_allowed_num_initial_clashes_with_global_context=1
    #Output 
    out_file_dir=args['out_file_dir']
    out_file_dir_debug=out_file_dir+"/debug"
    if not os.path.exists(out_file_dir):
        print "Out dir:\n %s\ndoesn't exist, generating it"
        os.makedirs(out_file_dir)
    if not os.path.exists(out_file_dir_debug):
        print "Debug out dir:\n %s\ndoesn't exist, generating it"
        os.makedirs(out_file_dir_debug)

    #Hacks
    keep_intact_chains_list=args['keep_intact_chains_list']
    force_chain_pairs=args['force_chain_pairs']     #e.g. [[2,6],
                                                    # [4,5],
                                                    # [6,7],] #Rosetta numbering of chains after splitting
                        
    b_stop_trimming_at_hotspots=True
    only_as_input_conectivity_order_mode=False  #True uses input order, False uses exhaustive ordering
    b_find_incomplete_elements_reconnections=False #For DEBUG is useful
    b_find_incomplete_elements_reconnections_must_contain_elements=[] #[1,2,6,7,10]

    b_generate_final_output_pose_mutating_all_non_hotspots_to_single_absurd_aa=False
    aa_for_mutating_all_non_hotspots_in_absurd_output="F"

    #SS-size
    min_allowed_ss_size={}  #dic
    min_allowed_ss_size['E']=3  #res
    min_allowed_ss_size['H']=5  #res
    min_allowed_ss_size['K']=1  #res (Defines a weird structure that the user wants to keep intact)

    #A fake AA one leter-code that the user will keep intact
    fake_aa_for_design='V'

    #HOTSPOT LABELS
    pdb_info_label_keyword="HOTSPOT"
    pdb_info_label_forceGap_keyword="FORCEGAP"
    pdb_info_label_forceH_keyword="FORCEH"
    pdb_info_label_forceE_keyword="FORCEE"

    #DLooper options
    #The clustered fragments Database:
    clustered_fragments_path=args['clustered_fragments_path']
    clustered_fragments_assignments_path=args['clustered_fragments_assignments_path']
    #Cluster control
    cluster_population_cut_off=100 #For 7aa/kcenters_stats.dat_res1.36, 100 == very good, 150==quite perfect, 200==pristine, etc... depends on the clustering
    check_spans_population_cut=20 #Has to be around 1/2 - 1/5 of "cluster_population_cut_off"

    #ToDo: make this to be automatically determined from the database
    max_rmsd_limit_fragments_spans_check=1.0 #Not more than 1.0? But, maybe is flexible depending on the clustering res
    nc_points_rmsd_dist_threshold=1.2  #1.5 permisive, 1.3 medium, 1.0 very-strict
    nc_points__individual_msd_dist_threshold=1.5 #1.5 permisive, 1.3 medium, 1.0 very-strict
    loop_NCpoints_distance_max_threshold_to_try=1.2 #This is misleading, it is actually CA CA distance right now, 
    allowed_omega_deviation=20.0  #Angle, 180 or 0 +- ; 20.0 OK

    #Cluster control
    min_allowed_loop_len=0 #Values smaller than 1 will be ignored
    max_allowed_loop_len=99
    max_clashes_perLoop_withMinimalSSscaffold=0 #Maybe 1?
    min_aa_prob_for_label=0.03 #P, Before: 0.05 #THIS IS SUPER_IMPORTANT FOR GENERATING the AA profiles used later for design

    #=max_rmsd_limit_fragments_spans_check #max distance to search for fragments similar 
    allowVariateSSsizeBy=1 #deprecated???
    fragment_blind_raddi_cut=0.75
    b_only_close_with_loop=False
    b_generate_loop_loop_closures=True
    break_trimming_at_hotspots=True
    b_debug_print=False
    b_debug_out=True

    #Loop Sampler Reducer options
    clustered_poses_res_threshold=1.0
    max_num_results_per_topology=args['max_num_results_per_topology']  #100 #1000
    fake_aminoacid_letter_for_design_loop='V'
    max_zero_CA_CA_distances=[4.0,3.0,2.0] #A

    clustered_loops_res_threshold=0.8 #Use < 1.0 A, it is selective for the loops or connection regions
    threshold_min_matches_for_clustered_loops=5
    b_make_single_membered_cluster_if_not_friends=True  #True == use clusters with population 1 too.
    maxNumLoops_per_sspair=10000
    maxNumLoops_per_sspair_variation=10000
    keep_only_top_loops_num_max_during_clustering=args['keep_only_top_loops_num_max_during_clustering']  #100 #300

    max_num_loops_per_hop_when_building_topology_filtA=args['max_num_loops_per_building_hop']*2 #20 #minLenght
    max_num_loops_per_hop_when_building_topology_filtB=args['max_num_loops_per_building_hop']   #10 #maxAvgDegree
    

    #Parametric SS size change
    #ToDo: 
    # 1. Different sizes for H and E
    target_reduc_size_Nterminus=args['target_reduc_size_Nterminus'] 
    target_reduc_size_Cterminus=args['target_reduc_size_Cterminus'] 
    target_extra_size_Cterminus=args['target_extra_size_Cterminus'] 
    target_extra_size_Nterminus=args['target_extra_size_Nterminus'] 

    #Parametric SS parameters
    omega_ideal_angle=180.0
    small_rmsd_pert_num=0.01

    #Parameters specific for Hs
    #ToDo: Implement restrictions in the sampled space based on know SS information
    max_perturbation_angle_H_phi=10.0  #15.0
    max_perturbation_angle_H_psi=10.0  #15.0
    max_perturbation_angle_H_omega=0.0
    perturbation_step_H_phi=2.0
    perturbation_step_H_psi=2.0
    perturbation_step_H_omega=1.0 #This parameter is irrelevant since we'll always set max_perturbation_angle_H_omega=0.0
    max_H_phi_curvature_pert=10.0  #15.0
    max_H_psi_curvature_pert=10.0  #15.0
    perturbation_step_H_phi_curvature=2.0
    perturbation_step_H_psi_curvature=2.0
    max_curv_pitch_aa_H_phi=5 #6-8
    max_curv_pitch_aa_H_psi=5 #6-8
    convergence_parametric_RMSD_H=0.5 #Best if <0.4
    max_parametric_RMSD_H=0.6 #Best if <0.5

    ###### START LAYER DESIGN DEFINITIONS #######
    #Dict and parameters used for Layer Profile of sequences
    #Core
    layer_aa_dic={}
    layer_aa_dic["core"]=['A','F','I','L','M','P','V'] #,'W'}   #,'Y'}
    layer_aa_dic["limbo"]=['A','D','E','F','G','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    layer_aa_dic["interface"]=['A','D','E','F','G','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    layer_aa_dic["surface"]=['D','E','G','K','N','P','Q','R','S','T'] #'H',
    #Special
    layer_aa_dic["special"]=['C','U']
    #All
    layer_aa_dic["all"]=['A','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','C']

    test_aa_large='I'
    test_aa_large_ref_sasa=187.0 #raddi=1.7
    test_aa_muta_contact='V'

    #ToDo: update to some kind of neighbor consideration
    layer_sasa_difinitions={}
    layer_sasa_difinitions["core"]=10.0 #% <=10.0-limbo
    layer_sasa_difinitions["limbo"]=25.0 #% <=25.5
    layer_sasa_difinitions["surf"]=40.0 #% <=40.0
    layer_sasa_difinitions["interface"]=1.0 #% <=, This is chainge in SASA with and without context

    layer_contact_distance_difinitions={}
    layer_contact_distance_difinitions["limbo"]=3.0 #Distance to search for contacts with the core.
    ###### END LAYER DESIGN DEFINITIONS #######


    #Other stuff needed
    in_pose_basename=os.path.basename(in_pdb_name)

    unique_run_id=hashlib.sha1(str.encode(str(random.random()))).hexdigest()[0:8]

    print "Done, unique random job ID:", unique_run_id

    #Read poses and make context-hash
    ####
    #START PDBloop here#
    ####
    #ToDo: Add break by SS

    print "Reading poses"
    inpose=pyrosetta.pose_from_pdb(in_pdb_name)
    print "Inpose size:", inpose.total_residue()
    print "Inpose size:", inpose.sequence()
    CA_coordinates=pose_atoms_to_nparray(inpose=inpose, 
                                         target_atoms=["CA"])

    #Find tagged residues
    print "Reading PDB-info residues tagged with the keyword: \"%s\""%(pdb_info_label_keyword)
    in_pdb_info=inpose.pdb_info()
    target_label_residues_ndx=[]
    force_gap_residues=[]
    force_h_residues=[]
    force_e_residues=[]
    for i in range(1, inpose.total_residue()+1):
        if pdb_info_label_keyword in in_pdb_info.get_reslabels(i):
            #Make them 0-index
            target_label_residues_ndx.append(i-1)
        if pdb_info_label_forceGap_keyword in in_pdb_info.get_reslabels(i):
            print "By-PDB-label-forcing a gap in res: %d"%i
            force_gap_residues.append(i)
        if pdb_info_label_forceH_keyword in in_pdb_info.get_reslabels(i):
            print "By-PDB-label-forcing ss H in res: %d"%i
            force_h_residues.append(i)
        if pdb_info_label_forceE_keyword in in_pdb_info.get_reslabels(i):
            print "By-PDB-label-forcing ss E in res: %d"%i
            force_e_residues.append(i)
    target_label_residues_ndx=np.asarray(target_label_residues_ndx)
    if (len(target_label_residues_ndx)<1):
        print "Couldn't find any target residues defined in the PDB with the keyword: \"%s\""%(pdb_info_label_keyword)
    else:
        print "The labeled hotspot residues are: ", target_label_residues_ndx+1

    lineout="select resid "
    for i in target_label_residues_ndx:
        lineout=lineout+"%d+"%(i+1)
    #Debug-pymolline:
    print   lineout


    #Calculate SS spans
    print "Calculating contigous SS spans"
    tmp_aaseq=inpose.sequence()
    #Now calculate contigous SS spans
    [b_has_enough_ss, 
     tmp_dssp_sequence, 
     contiguous_SS_spans, 
     contiguous_SS_spans_type] = calculate_disembodied_ss_spans(inpose=inpose, 
                                                                min_allowed_ss_size_dic=min_allowed_ss_size, 
                                                                discard_positions=force_gap_residues,
                                                                force_h_positions=force_h_residues,
                                                                force_e_positions=force_e_residues,
                                                                b_break_by_chain=True,
                                                               keep_chains=keep_intact_chains_list,
                                                               b_detect_breaks=True,
                                                               propagate_neighbors_NC_ss=False,
                                                               detect_E_bulges=False,#True,
                                                               NC_threshold_dist=1.7,
                                                               bdebug=False)

    print keep_intact_chains_list

    print "SS-seq:", tmp_dssp_sequence
    print "Spans:", contiguous_SS_spans
    print "Types:", contiguous_SS_spans_type

    if (not b_has_enough_ss):
        assert(0==1)

    inpose_split_by_ss=[]
    for ispan in xrange(len(contiguous_SS_spans)):
        tmp_pose=extract_pose_spans_as_chains(inpose = inpose, 
                                                 spans = [contiguous_SS_spans[ispan]]) #rosetta numbers

        #ORI version
        inpose_split_by_ss.append(tmp_pose.clone())

        """
        #DEBUG
        tmp_filename_ori="%s/testSSori_ss%02d.pdb"%(out_file_dir,
                                                     (ispan+1))
        tmp_pose.dump_pdb(tmp_filename_ori)
        inpose_split_by_ss.append(pyrosetta.pose_from_pdb(tmp_filename_ori))

        #ENDDEBUG
        """




        if b_debug_out:
            tmp_filename_ori="%s/debug/%s_testSSori_ss%02d.pdb"%(out_file_dir,
                                                        unique_run_id,
                                                     (ispan+1))
            tmp_pose.dump_pdb(tmp_filename_ori)
            print "\t\tYAY!, SS%d %s pose... OUT: %s" % ((ispan+1),
                                                         contiguous_SS_spans_type[ispan], 
                                                         tmp_filename_ori)

    ##print pose_split_by_ss
    if (args['enable_only_profile_mode']):
        print "!!!\n!!!\n!!!Terminating program!!! (i.e. only_profile_input is activated)\n!!!\n!!!"
        return

    print "Done reading->writting protein parts"
    #print "TERMINATING"
    #return #THIS IS JUST FOR PROFILING SO WE'LL EXIT HERE

    ##dlooper_alg=None
    #Instantiate the fantastic :)) DS looper
    if dlooper_alg is None:
        #some_fallback_operation()
        print "Using databases \n  %s\n  %s:"%(clustered_fragments_path,clustered_fragments_assignments_path)

        print "Initializing the fantastic DLooper!!!"

        dlooper_alg=dlooper( nc_points_rmsd_dist_threshold=nc_points_rmsd_dist_threshold,
                             nc_points_individual_msd_dist_threshold=nc_points__individual_msd_dist_threshold,
                             loop_NCpoints_distance_max_threshold_to_try=loop_NCpoints_distance_max_threshold_to_try,
                             max_rmsd_limit_spans_check=max_rmsd_limit_fragments_spans_check,
                             allowed_omega_deviation=allowed_omega_deviation,
                             cluster_population_cut_off=cluster_population_cut_off, #50=good, 500=suuuuper perfect
                             check_spans_population_cut=check_spans_population_cut, #30=good, 100=suuuuper goood
                             fragment_blind_raddi_cut=fragment_blind_raddi_cut,
                             fake_aminoacid_letter_for_design="A",
                             min_aa_prob_for_label=min_aa_prob_for_label,
                             max_zero_CA_CA_distances=max_zero_CA_CA_distances, #
                             b_only_close_with_loop=b_only_close_with_loop,
                             b_generate_loop_loop_closures=b_generate_loop_loop_closures,
                             b_debug_print=b_debug_print,
                             b_debug_out=b_debug_out,
                             clustered_fragments_path=clustered_fragments_path,
                             clustered_fragments_assignments_path=clustered_fragments_assignments_path,
                             debug_out_path="%s/debug/%s_"%(out_file_dir, unique_run_id))
    else:
        print "Fantastic DLooper exists already, skipping for the sake of your mind :)!!!"

    #Read context
    global_context_pose=None
    global_context_pose_hash=None
    global_context_pose_basename=None
    num_global_context_clashes_initial=0
    if b_use_global_context:
        print "Per your request I am using global context for clash check from: ", in_global_context_path
        print "Note: IF your input is ab-initio clashing... Well you'll get in trouble. I'll better check for you first! WAIT"
        global_context_pose=pyrosetta.pose_from_pdb(in_global_context_path)
        in_context_pose_basename=os.path.basename(in_global_context_path)
        global_context_pose_hash=pyrosetta.rosetta.core.pose.xyzStripeHashPose(global_context_pose,
                                                                    pyrosetta.rosetta.core.pose.PoseCoordPickMode_N_CA_C_CB, #maybe .BB instead of .HVY
                                                                    radius=2.0)
        has_GC_clashes,GC_clash_map=dlooper_alg.check_pm_clashes_hash_vs_pose(global_context_pose_hash,
                                                                            inpose,
                                                                            pick_mode=pyrosetta.rosetta.core.pose.PoseCoordPickMode_N_CA_C_CB)
        num_global_context_clashes_initial=len(list(GC_clash_map))
        if has_GC_clashes:
            print "HAHA! Something is weird, your input is already clashing with the target: %d"%num_global_context_clashes_initial
            if (num_global_context_clashes_initial > max_allowed_num_initial_clashes_with_global_context):
                print "According to max_allowed_num_initial_clashes_with_global_context these are too many clashes. Fix before running this. I'll stop now"
                assert(0==1)
            else:
                print "OK, I'll continue, but this is weird and you should check your pose in the future before comming here :)"
        else:
            print "Alright, I don't see clashes to start with. Let's proceed then"
            num_global_context_clashes_initial=0
    print "All Done"
    #End read-context

    #Calculate the best (grid) SS parameters
    #ToDo:
    #Compare Hotspots to SS-fragments and warn if there is anything outside of the spans
    #ToDo 1: Add Rama space check (it only make's sense to do this inside of the rama space!)
    #ToDo 2: add multithreading 
    #ToDo 3: Make it faster, maybe precalculating stuff
    #ToDo 4: 'E' (strand) SSs equation discovery should happen by pairs or tripplets of bb H-bonded structures
    #ToDo 5: Add resume points.

    print "I'll scan parameters to generate your SS structures. Hold tight, this might take long!"

    best_parameters_per_SSndx={} #collections.default_dic(list)
    for indx in xrange(len(contiguous_SS_spans)):
        ispan=contiguous_SS_spans[indx]
        tmp_in_pose=inpose_split_by_ss[indx].clone()
        print "\nWorking on span (%d of %d):"%(indx+1, len(contiguous_SS_spans)),  ispan, ", of type: ", contiguous_SS_spans_type[indx]
        best_sampled_rmsd_for_parameters=999999.99
        if (contiguous_SS_spans_type[indx] == 'H' ):
            print "Yay!!! Sampling a helix... These ones are usually easy, unless the are really kinked."
            ori_phi_angles=[]
            ori_psi_angles=[]
            ori_omega_angles=[]
            for ires in xrange(tmp_in_pose.total_residue()):
                ori_phi_angles.append(tmp_in_pose.phi(ires+1))
                ori_psi_angles.append(tmp_in_pose.psi(ires+1))
                ori_omega_angles.append(abs(tmp_in_pose.omega(ires+1)))
                ##print ori_phi_angles[-1], ori_psi_angles[-1], ori_omega_angles[-1]
            #Calculate best angles
            assert((max_curv_pitch_aa_H_phi>=0)and(max_curv_pitch_aa_H_psi>=0))
            for i_curv_pitch_aa_H_phi in range(min(1, max_curv_pitch_aa_H_phi), 
                                               max_curv_pitch_aa_H_phi+1):
                for j_curv_pitch_aa_H_psi in range(min(1, max_curv_pitch_aa_H_psi),
                                                   max_curv_pitch_aa_H_psi+1):
                    grid_sampled_phi_psi_vs_rmsd=sample_helix_fitting_parameters(target_pose=tmp_in_pose,
                                                                                avg_phi=np.around(np.average(ori_phi_angles), 
                                                                                                     decimals=1),
                                                                                avg_psi=np.around(np.average(ori_psi_angles), 
                                                                                                     decimals=1),
                                                                                avg_omega=np.around(omega_ideal_angle),  #Harcoded ideal==180.0
                                                                                max_perturbation_angle_phi=max_perturbation_angle_H_phi,
                                                                                max_perturbation_angle_psi=max_perturbation_angle_H_psi,
                                                                                max_perturbation_angle_omega=max_perturbation_angle_H_omega,
                                                                                max_phi_curvature_pert=max_H_phi_curvature_pert,
                                                                                max_psi_curvature_pert=max_H_psi_curvature_pert,
                                                                                perturbation_step_phi=perturbation_step_H_phi,
                                                                                perturbation_step_psi=perturbation_step_H_psi,
                                                                                perturbation_step_omega=perturbation_step_H_omega,
                                                                                perturbation_step_phi_curvature=perturbation_step_H_phi_curvature,
                                                                                perturbation_step_psi_curvature=perturbation_step_H_psi_curvature,
                                                                                curv_phi_pitch_aa=i_curv_pitch_aa_H_phi,
                                                                                curv_psi_pitch_aa=j_curv_pitch_aa_H_psi,
                                                                                convergence_RMSD=convergence_parametric_RMSD_H,
                                                                                fake_aa_for_design=fake_aa_for_design,
                                                                                debug=False)

                    tmp_best_parameters_per_SSndx=grid_sampled_phi_psi_vs_rmsd[grid_sampled_phi_psi_vs_rmsd[:,0].argmin()]
                    #Store the best parameters found
                    if (tmp_best_parameters_per_SSndx[0]<best_sampled_rmsd_for_parameters):
                        best_sampled_rmsd_for_parameters=tmp_best_parameters_per_SSndx[0]
                        best_parameters_per_SSndx[indx]=np.copy(tmp_best_parameters_per_SSndx)
                        print "  Best Parameters found until now: ", best_parameters_per_SSndx[indx]
                    if(best_sampled_rmsd_for_parameters<=convergence_parametric_RMSD_H):
                        break
                if(best_sampled_rmsd_for_parameters<=convergence_parametric_RMSD_H):
                    print "  RMSD (%0.2f) converged to desired parameter <=(%0.2f), stopping sampling"%(best_sampled_rmsd_for_parameters,
                                                                                                        convergence_parametric_RMSD_H)
                    break

            #best_parameters_per_SSndx[indx]=(grid_sampled_phi_psi_vs_rmsd[grid_sampled_phi_psi_vs_rmsd[:,0].argmin()])
            print "*Best Parameters found for ss %d: "%(indx+1), best_parameters_per_SSndx[indx]
            if (best_parameters_per_SSndx[indx][0] > max_parametric_RMSD_H):
                for krand in xrange(2):
                    print "*OhOh!!!: The RMSD value (%0.2f) is too high compared to your limits... You have been warned (twice): " % best_parameters_per_SSndx[indx][0]

        elif(contiguous_SS_spans_type[indx] == 'E'):
            print "Yay!!! Sampling a beta strand... These ones are hard, let's see..."
            print "!!!ERROR... OhOh, not implemented...Contact Daniel... I'll panic now"
            assert(0)
                           
        elif (contiguous_SS_spans_type[indx] == 'K' ):
            print "I'll keep this structure intact :) as per your request"
            best_parameters_per_SSndx[indx]=['K']

        else:
            print "I can't do this kind of SS (%s), sorry dude...."%contiguous_SS_spans_type[indx]
            assert (0==1)

    print "Done"



    #ToDo : add read from disk parameters (i.e. resume behaviour)
    #Write parameters to file
    print "Writing ideal parameters to file"
    tmp_out_parameters_filename="%s/%s_%s.ideal_params"%(out_file_dir,
                                                         unique_run_id,
                                         in_pose_basename[0:-4])
    with open(tmp_out_parameters_filename, 'wb') as handle:
        pickle.dump(best_parameters_per_SSndx, handle)
    print "Done"




    inpose_split_by_ss_optimized_position=[]
    for ipose in inpose_split_by_ss:
        inpose_split_by_ss_optimized_position.append(ipose.clone())

    print("DONE")   


    #ToDo: Change this code, just generate the largest possible size and ake sub-copies of all the possible alternations 
    #Starting by the smallest and checking clashes with the target/context :)
    ordered_SSs_variations=[]
    for indx in best_parameters_per_SSndx.keys():
        #indx is the ss number
        ispan=contiguous_SS_spans[indx]
        print "Working on rebuilding input span: ", ispan
        tmp_pose_ori_ss_segment=inpose_split_by_ss_optimized_position[indx].clone() 
        tmp_pose_ori_ss_segment_size=tmp_pose_ori_ss_segment.total_residue()
        #print "A: ", indx, tmp_pose_ori_ss_segment_size
        #Container for results
        ordered_SSs_variations.append([])

        if (contiguous_SS_spans_type[indx] == 'K'): #Just keep whatever the user input here
            tmp_fake_pose=tmp_pose_ori_ss_segment.clone()


            #Check to remove pyrosetta.rosetta.core.chemical.VariantType.LOWER_TERMINUS_VARIANT variant type (in case)
            for new_res_ndx in xrange(1, (tmp_fake_pose.total_residue()+1)):
                if tmp_fake_pose.residue(new_res_ndx).has_variant_type(pyrosetta.rosetta.core.chemical.VariantType.LOWER_TERMINUS_VARIANT):
                    pyrosetta.rosetta.core.pose.remove_variant_type_from_pose_residue(tmp_fake_pose, pyrosetta.rosetta.core.chemical.VariantType.LOWER_TERMINUS_VARIANT, 
                                                                            new_res_ndx)
                #Check to remove pyrosetta.rosetta.core.chemical.VariantType.UPPER_TERMINUS_VARIANT variant type (in case)
                elif tmp_fake_pose.residue(new_res_ndx).has_variant_type(pyrosetta.rosetta.core.chemical.VariantType.UPPER_TERMINUS_VARIANT):
                    pyrosetta.rosetta.core.pose.remove_variant_type_from_pose_residue(tmp_fake_pose, pyrosetta.rosetta.core.chemical.VariantType.UPPER_TERMINUS_VARIANT, 
                                                                            new_res_ndx)

            print "SS ndx=%d is of type K(of keep)! So, I'll still generate trimmings for this (only)"%(indx+1)
            #Get the pdb_info labels
            tmp_pdb_info=tmp_fake_pose.pdb_info()
            b_previous_was_hotspot_n=False
            for npoint_trim_size in xrange(target_reduc_size_Nterminus+1):
                new_npoint_start_ndx=npoint_trim_size+1 #Rosetta numbers

                if(b_previous_was_hotspot_n and b_stop_trimming_at_hotspots):
                    print "\t\tStop this N- variation trim by the Do-Not-Remove-Hotspot-Rule (var_%d/res_%d)."%( npoint_trim_size, new_npoint_start_ndx)
                    break
                #Check if this n-side is a hotspot
                if (pdb_info_label_keyword in tmp_pdb_info.get_reslabels(new_npoint_start_ndx)):
                        print "STOP AT NEXT ROUND N-side trimming"
                        b_previous_was_hotspot_n=True
                b_previous_was_hotspot_c=False
                for cpoint_trim_size in xrange(target_reduc_size_Nterminus+1):
                    new_cpoint_end_ndx=tmp_fake_pose.total_residue()-cpoint_trim_size #Rosetta numbers

                    if(b_previous_was_hotspot_c and b_stop_trimming_at_hotspots):
                        print "\t\tStop this C- variation trim by the Do-Not-Remove-Hotspot-Rule (var_%d/res_%d)."%( cpoint_trim_size, new_cpoint_end_ndx)
                        break
                    #Check if this c-side is a hotspot
                    if (pdb_info_label_keyword in tmp_pdb_info.get_reslabels(new_cpoint_end_ndx)):
                        print "STOP AT NEXT ROUND C-side trimming"
                        b_previous_was_hotspot_c=True

                    #New pose size
                    target_ss_size_new=tmp_fake_pose.total_residue()-npoint_trim_size-cpoint_trim_size
                    if b_debug_out:
                        print "Generating variation for fragment %d:"%(indx+1), new_npoint_start_ndx, new_cpoint_end_ndx, target_ss_size_new
                    #Create the pose container
                    tmp_fake_seq=""
                    for ires in  np.arange( target_ss_size_new ):
                        tmp_fake_seq+=fake_aa_for_design
                    tmp_fake_pose_final=pyrosetta.pose_from_sequence(tmp_fake_seq)

                    tmp_fake_pose_final.copy_segment(target_ss_size_new,
                                                     tmp_fake_pose,
                                                     1,
                                                     new_npoint_start_ndx)

                    #Copy Labels
                    #ToDO, double check this!!! (Seems good, let's see)
                    for ires in xrange(target_ss_size_new):
                        for pdbinfolabel in tmp_fake_pose.pdb_info().get_reslabels(new_npoint_start_ndx+ires):
                            #print "TeST:", (ires+1), pdbinfolabel #, (ires+1)
                            tmp_fake_pose_final.pdb_info().add_reslabel((ires+1),
                                                                          pdbinfolabel)

                    tmp_ss_type=np.zeros((tmp_fake_pose_final.total_residue()),int)
                    #ToDo: Shall we use the real SS type here?
                    tmp_ss_type.fill(3) #imaginary SS type


                    ordered_SSs_variations[-1].append( [tmp_fake_pose_final.clone(),
                                                        ispan, 
                                                        -npoint_trim_size, #In this case has to be negative because we are cutting
                                                        -cpoint_trim_size, #In this case has to be negative because we are cutting
                                                        tmp_ss_type[:]] )



        elif (contiguous_SS_spans_type[indx] == 'H' or contiguous_SS_spans_type[indx] == 'E'):
            best_RMSD=best_parameters_per_SSndx[indx][0]
            best_phi=best_parameters_per_SSndx[indx][1]
            best_psi=best_parameters_per_SSndx[indx][2]
            best_omega=best_parameters_per_SSndx[indx][3]
            best_phi_curvature_pert=best_parameters_per_SSndx[indx][4]
            best_psi_curvature_pert=best_parameters_per_SSndx[indx][5]
            best_phi_curvature_pert_pitch=best_parameters_per_SSndx[indx][6]
            best_psi_curvature_pert_pitch=best_parameters_per_SSndx[indx][7]


            #Build the largest pose
            target_ss_size_new=ispan[1]-ispan[0]+target_extra_size_Cterminus+target_extra_size_Nterminus+1
            #target_ss_size_oriPose=ispan[1]-ispan[0]+1

            #Create the pose container +2 fake residues
            tmp_fake_seq=""
            for ires in  np.arange( target_ss_size_new+2 ):
                tmp_fake_seq+=fake_aa_for_design
            tmp_fake_pose=pyrosetta.pose_from_sequence(tmp_fake_seq)

            #+2 due to the hanging bonds-residues
            for ires in np.arange(target_ss_size_new):
                internal_num=ires+2
                rosetta_num=ires-target_extra_size_Nterminus+1
                ##print rosetta_num, internal_num
                if ((rosetta_num%best_phi_curvature_pert_pitch)==0):
                    tmp_fake_pose.set_phi(internal_num,
                                          best_phi+best_phi_curvature_pert)
                else:
                    tmp_fake_pose.set_phi(internal_num,
                                          best_phi)
                if((rosetta_num%best_psi_curvature_pert_pitch)==0):
                    tmp_fake_pose.set_psi(internal_num,
                                          best_psi+best_psi_curvature_pert)
                else:
                    tmp_fake_pose.set_psi(internal_num,
                                          best_psi)
                tmp_fake_pose.set_omega(internal_num,
                                        best_omega)




            #Sanity check #1: compare RMSD
            rmsdVal=align_atoms_by_ndxs(tmp_fake_pose, 
                                                     2+target_extra_size_Nterminus, 
                                                     (tmp_pose_ori_ss_segment_size+target_extra_size_Nterminus+1), 
                                                     tmp_pose_ori_ss_segment, #inpose, 
                                                     1, #ispan[0]+1, 
                                                     tmp_pose_ori_ss_segment_size,
                                                     atoms=["CA","C","O","N"]) #ispan[1]+1)

            if (abs(rmsdVal-best_RMSD) <= small_rmsd_pert_num):
                print "\tAll seems OK, the RMSD difference to target is %0.3f (VS expected %0.3f)"%(rmsdVal,best_RMSD)
            else:
                print "\tSomething smells bad here, the RMSD difference to the target structure for this span is too high %0.3f (VS expected %0.3f)"%(rmsdVal,best_RMSD)


            #Copy labeled-residues/hotspots from the original pose
            print "\tChecking Hotspots transfer" ##, target_ss_size_new
            for ires in np.arange(target_ss_size_new):
                new_res_ndx=ires+2

                ori_pose_rosetta_num=ires-target_extra_size_Nterminus+1
                if ( (ori_pose_rosetta_num<1) or 
                     (ori_pose_rosetta_num>tmp_pose_ori_ss_segment_size)): #Check that we are not out of this span :)
                    continue
                tmp_pose_ori_ss_segment_pdb_info=tmp_pose_ori_ss_segment.pdb_info()
                tmp_res_labels=tmp_pose_ori_ss_segment_pdb_info.get_reslabels(ori_pose_rosetta_num)
                #Copy residue's identities
                b_hotspots_added=False
                if(len(tmp_res_labels) > 0):
                    if pdb_info_label_keyword in tmp_res_labels:
                        b_hotspots_added=True
                        #print tmp_res_labels
                        print " Copying labeled residue SC %d to position %d"%( ori_pose_rosetta_num, 
                                                                   new_res_ndx-1)
                        tmp_fake_pose.replace_residue(new_res_ndx,
                                                      tmp_pose_ori_ss_segment.residue(ori_pose_rosetta_num),
                                                      True)
                        #Check to remove pyrosetta.rosetta.core.chemical.VariantType.LOWER_TERMINUS_VARIANT variant type (in case)
                        if tmp_fake_pose.residue(new_res_ndx).has_variant_type(pyrosetta.rosetta.core.chemical.VariantType.LOWER_TERMINUS_VARIANT):
                            pyrosetta.rosetta.core.pose.remove_variant_type_from_pose_residue(tmp_fake_pose, pyrosetta.rosetta.core.chemical.VariantType.LOWER_TERMINUS_VARIANT, 
                                                                                    new_res_ndx)
                        #Check to remove pyrosetta.rosetta.core.chemical.VariantType.UPPER_TERMINUS_VARIANT variant type (in case)
                        elif tmp_fake_pose.residue(new_res_ndx).has_variant_type(pyrosetta.rosetta.core.chemical.VariantType.UPPER_TERMINUS_VARIANT):
                            pyrosetta.rosetta.core.pose.remove_variant_type_from_pose_residue(tmp_fake_pose, pyrosetta.rosetta.core.chemical.VariantType.UPPER_TERMINUS_VARIANT, 
                                                                                    new_res_ndx)
                    #Re-add Reslabels
                    for pdbinfolabel in tmp_res_labels:
                        tmp_fake_pose.pdb_info().add_reslabel(new_res_ndx,
                                                          pdbinfolabel)

            #If we added/copied hotspots make a Sanity check #2: compare RMSD
            if b_hotspots_added:
                print "\tHotspots were transfered, checking RMSD again"
                rmsdVal=align_atoms_by_ndxs(tmp_fake_pose, 
                                                     2+target_extra_size_Nterminus, 
                                                     (tmp_pose_ori_ss_segment_size+target_extra_size_Nterminus+1), 
                                                     tmp_pose_ori_ss_segment, #inpose, 
                                                     1, #ispan[0]+1, 
                                                     tmp_pose_ori_ss_segment_size,
                                                      atoms=["CA","C","O","N"]) #ispan[1]+1)

                if (abs(rmsdVal-best_RMSD) <= small_rmsd_pert_num):
                    print "\tAll seems OK (check #2 after Hotspots), the RMSD difference to target is %0.3f (VS expected %0.3f)"%(rmsdVal,best_RMSD)
                else:
                    print "\tSomething smells bad here (check #2 after Hotspots), the RMSD difference to the target structure for this span is too high %0.3f (VS expected %0.3f)"%(rmsdVal,best_RMSD)
                    #assert (0==1)

            #Make n-terminal size variations
            #ToDO:
            # Correct this, stop at the hotspot!!!
            b_previous_was_hotspot_n=False
            for i_nterm_variation in xrange(target_extra_size_Nterminus+target_reduc_size_Nterminus+1):
                new_pose_rosetta_num_N=i_nterm_variation+2
                size_variation_n=target_extra_size_Nterminus-i_nterm_variation
                if(b_previous_was_hotspot_n and b_stop_trimming_at_hotspots):
                        print "\t\tStop this N- variation trim by the Do-Not-Remove-Hotspot-Rule (var_%d/res_%d)."%( size_variation_n, new_pose_rosetta_num_N-1)
                        break
                #Double check this: new_pose_rosetta_num_N+1
                b_previous_was_hotspot_c=False 
                for j_cterm_variation in xrange(target_extra_size_Cterminus+target_reduc_size_Cterminus+1):
                    new_pose_rosetta_num_C=target_ss_size_new-j_cterm_variation+1
                    size_variation_c=target_extra_size_Cterminus-j_cterm_variation
                    if(b_previous_was_hotspot_c and b_stop_trimming_at_hotspots):
                        print "\t\tStop this -C variation trim by the Do-Not-Remove-Hotspot-Rule (var_%d/res_%d)."%( size_variation_c, new_pose_rosetta_num_C-1)
                        break
 
                    #Check pose size regarding the original SS:
                    ori_ss_size_after_trim=tmp_pose_ori_ss_segment_size+(min(0,size_variation_n))+(min(0,size_variation_c))
                    if (ori_ss_size_after_trim < min_allowed_ss_size[contiguous_SS_spans_type[indx]]):
                        print "\t\tOhoh. I realize that I am making the original (%s) too small (%d vs limit= %d), skipping this variation: N%d_C%d"%(contiguous_SS_spans_type[indx],
                                                                                                      ori_ss_size_after_trim,
                                                                                                      min_allowed_ss_size[contiguous_SS_spans_type[indx]],
                                                                                                      size_variation_n, 
                                                                                                      size_variation_c)
                        continue

                    #Create the finale pose container (without the +2 fake residues)
                    tmp_fake_final_pose_size=new_pose_rosetta_num_C-new_pose_rosetta_num_N+1
                    assert( tmp_fake_final_pose_size > 0)
                    print "\t\tOK, generating pose for size variation: N%d_C%d"%(size_variation_n, size_variation_c)
                    #Pose container
                    tmp_fake_seq_final=""
                    for ires in xrange(tmp_fake_final_pose_size):
                        tmp_fake_seq_final+=fake_aa_for_design
                    tmp_fake_pose_final=pyrosetta.pose_from_sequence(tmp_fake_seq_final)
                    tmp_fake_pose_final.copy_segment(tmp_fake_final_pose_size,
                                                     tmp_fake_pose,
                                                     1,
                                                     new_pose_rosetta_num_N)

                    #Copy Labels
                    #ToDO, double check this!!! (Seems good, let's see)
                    for ires in xrange(tmp_fake_final_pose_size):
                        for pdbinfolabel in tmp_fake_pose.pdb_info().get_reslabels(new_pose_rosetta_num_N+ires):
                            tmp_fake_pose_final.pdb_info().add_reslabel((ires+1),
                                                                          pdbinfolabel)



                    #Check clashes
                    if b_use_global_context:
                        print "\t\tChecking for clashes against the context"
                        [b_hasBBclashes, target_clash_map] = dlooper_alg.check_pm_clashes_hash_vs_pose(global_context_pose_hash, 
                                                                                        tmp_fake_pose_final,
                                                                                        pick_mode=pyrosetta.rosetta.core.pose.PoseCoordPickMode_N_CA_C_CB)
                        #ToDo: The elegant/efficient way to do this would be to prune the size-"Increases" as soon as you hit a clash :_
                        if( b_hasBBclashes ): 
                            if ((len(list(target_clash_map))-num_global_context_clashes_initial) > 0 ):
                                print "\t\tOhNooo... This variation has bb-clashes.. I'll have to skip it"
                                continue
                        else:
                            print "\t\tGood, no bb-clashes detected. Continuing..."
                    #Finally append to the results
                    print "SS assignments by chain: (0=H,1=E,2=L, others=???)"
                    ss_pispred_assignments_by_chain=dlooper_alg.calculate_ss_assignments_by_dssp(in_pose_by_chain=[tmp_fake_pose_final],
                                                                                                 b_flatten_ss_assignments=False )
                    for issAssignment in ss_pispred_assignments_by_chain:
                        print " -", issAssignment
                    ordered_SSs_variations[-1].append( [tmp_fake_pose_final.clone(),
                                                        ispan, 
                                                        size_variation_n, 
                                                        size_variation_c,
                                                        ss_pispred_assignments_by_chain[0]] )
                    #Check Hotspot existence here!
                    #Check for touching hotspots
                    tmp_pdb_info=tmp_fake_pose_final.pdb_info()
                    print tmp_pdb_info.get_reslabels(1), tmp_pdb_info.get_reslabels(tmp_fake_pose_final.total_residue())
                    if (pdb_info_label_keyword in tmp_pdb_info.get_reslabels(1)):
                        print "STOP AT NEXT ROUND N", tmp_fake_pose_final.total_residue()
                        b_previous_was_hotspot_n=True
                    if (pdb_info_label_keyword in tmp_pdb_info.get_reslabels(tmp_fake_pose_final.total_residue())):
                        print "STOP AT NEXT ROUND C", tmp_fake_pose_final.total_residue()
                        b_previous_was_hotspot_c=True

        #Check that we have at least one result in the pool for this SS
        assert (len(ordered_SSs_variations[-1]) > 0)

    print "ALL done"

    #Prunning by clash check of a-given-SS vs all-other-SSs smallest versions minus ends.
    print "Prunning SSs by clashes agains minimal SS core"

    num_ideal_ss=len(ordered_SSs_variations)
    for ndxSSA in xrange(num_ideal_ss):
        print "Prunning by clashes SS#", ndxSSA, ". With ori SS element size=", len(ordered_SSs_variations[ndxSSA])
        comparision_minSS_context_pose=pyrosetta.pose_from_sequence("")
        for ndxSSB in xrange(num_ideal_ss):
            #skip current SS
            if (ndxSSA == ndxSSB):
                continue
            if (comparision_minSS_context_pose.total_residue() == 0):
                #0=pose, 1=span, 2=NextraSize, 3=CextraSize
                tmp_pose=return_pose_copy_minus_ends(ordered_SSs_variations[ndxSSB][-1][0].clone())
                comparision_minSS_context_pose=tmp_pose.clone() #-1 ndx is always the smallest combination in this code
            else:
                tmp_pose=return_pose_copy_minus_ends(ordered_SSs_variations[ndxSSB][-1][0].clone())
                comparision_minSS_context_pose.append_pose_by_jump(tmp_pose, 1)

        comparision_minSS_context_pose_hash = pyrosetta.rosetta.core.pose.xyzStripeHashPose(comparision_minSS_context_pose,
                                                            pyrosetta.rosetta.core.pose.PoseCoordPickMode_N_CA_C_CB,
                                                            radius=2.0)

        tmp_keep_ndxs=[]
        for i_comb in xrange(len(ordered_SSs_variations[ndxSSA])):

            [do_clash,num_clashes]=dlooper_alg.check_pm_clashes_hash_vs_pose(comparision_minSS_context_pose_hash, 
                                                                              ordered_SSs_variations[ndxSSA][i_comb][0],
                                                                              pick_mode=pyrosetta.rosetta.core.pose.PoseCoordPickMode_N_CA_C_CB)
            print ordered_SSs_variations[ndxSSA][i_comb], do_clash, len(list(num_clashes))
            if (len(num_clashes)<3):
                #print i_comb
                tmp_keep_ndxs.append(i_comb)

        assert (len(tmp_keep_ndxs) > 0)

        print "Keeping elements:", tmp_keep_ndxs
        if (len(tmp_keep_ndxs)==1):
            ordered_SSs_variations[ndxSSA]=[itemgetter(*tmp_keep_ndxs)(ordered_SSs_variations[ndxSSA])]
        else:
            ordered_SSs_variations[ndxSSA]=itemgetter(*tmp_keep_ndxs)(ordered_SSs_variations[ndxSSA])

        print " New size =", len(ordered_SSs_variations[ndxSSA])

        #Debug Pose result print
        ##"""
        if b_debug_out:
            print "Out of debug-poses to: ", out_file_dir+"/testSSide*.pdb"
            for i_comb in xrange(len(ordered_SSs_variations[ndxSSA])):
                #print ordered_SSs_variations[ndxSSA][i_comb]
                tmp_filename_new="%s/%s_testSSide_ss%02d_ss%s_varN%02d_varC%02d.pdb"%(out_file_dir+"/debug",
                                                                                   unique_run_id,
                                                                                (ndxSSA+1),      
                                                                            str(ordered_SSs_variations[ndxSSA][i_comb][1]).replace(" ","_"),
                                                                            ordered_SSs_variations[ndxSSA][i_comb][2],
                                                                            ordered_SSs_variations[ndxSSA][i_comb][3])

                ordered_SSs_variations[ndxSSA][i_comb][0].dump_pdb(tmp_filename_new)
                #print "\t\tYAY!, pose... OUT:", tmp_filename_new
        ##"""

    print "All prunning Done"


    #Now Put loops
    #ToDo 1 : add multithreading
    #ToDo 2 : pass a reliable SS definition because if not it is lost without the context. Or fix it somehow else
    #ToDo 3 : use info in ToDo 2.
    #ToDo 4 :Loops of len=0 (i.e. closed with a 2aa loop, will give raise to a large number of identical solutions).
    # Might be fixable by copying the loop from the beggining to the end (join3poses function), but will introduce
    # Further problems in preserving hotspots

    #DEBUG HACKS (tmp):
    #maxNumLoops_per_sspair=3
    #maxNumLoops_per_sspair_variation=3
    #dlooper_alg.b_debug_print=False
    ###dlooper_alg.b_debug_out=False
    ###only_as_input_conectivity_order_mode=False
    #END


    #dlooper_alg.nc_points_rmsd_dist_threshold=1.2
    #dlooper_alg.nc_points_individual_msd_dist_threshold=1.5
    #dlooper_alg.loop_NCpoints_distance_max_threshold_to_try=1.2

    print "Connecting by... pairs??"
    #Create the containers for the results
    num_ideal_ss=len(ordered_SSs_variations)
    looped_results_matrix=np.zeros((num_ideal_ss, num_ideal_ss), int)
    looped_results_matrix.fill(0)
    looped_result_poses_container=[]

    for i in xrange(num_ideal_ss):
        looped_result_poses_container.append([])
        for j in xrange(num_ideal_ss):
            looped_result_poses_container[i].append([])

    sol_counter=0
    for ndxSSA in xrange(num_ideal_ss):
        for ndxSSB in xrange(num_ideal_ss):
            if (ndxSSA == ndxSSB):
                continue
            elif (only_as_input_conectivity_order_mode):
                if (ndxSSA != (ndxSSB-1)):
                    continue
            #Forced Particular connectivities check:
            is_forced_connection_broken=False
            for kndxpair in force_chain_pairs:
                if (((ndxSSA+1) == kndxpair[0]) or 
                    ((ndxSSB+1) == kndxpair[1]) ):
                    print "Connections to this chain are restricted by pair, checking:", kndxpair
                    if not ([ndxSSA+1,ndxSSB+1]==kndxpair):
                        is_forced_connection_broken=True
                        break
            if is_forced_connection_broken:
                print "Skipping search for this pair because only one of the chains of this pair is in the defined in", force_chain_pairs
                continue
            #End forced connectivity

            #Pose hash to check agains the minimal SS structure:
            comparision_minSS_context_pose=pyrosetta.pose_from_sequence("")
            for ndxSSC in xrange(num_ideal_ss):
                #skip current SS
                if ( (ndxSSA == ndxSSC) or (ndxSSB == ndxSSC) ):
                    continue
                if (comparision_minSS_context_pose.total_residue() == 0):
                    #0=pose, 1=span, 2=NextraSize, 3=CextraSize
                    tmp_pose=return_pose_copy_minus_ends(ordered_SSs_variations[ndxSSC][-1][0].clone())
                    comparision_minSS_context_pose=tmp_pose.clone() #-1 ndx is always the smallest combination in this code
                else:
                    tmp_pose=return_pose_copy_minus_ends(ordered_SSs_variations[ndxSSC][-1][0].clone())
                    comparision_minSS_context_pose.append_pose_by_jump(tmp_pose, 1)

            comparision_minSS_context_pose_hash = pyrosetta.rosetta.core.pose.xyzStripeHashPose(comparision_minSS_context_pose,
                                                                                                pyrosetta.rosetta.core.pose.PoseCoordPickMode_N_CA_C_CB,
                                                                                                radius=2.0)


            print "Working on connecting ss_%02d <-> %02d: "%(ndxSSA+1,
                                                              ndxSSB+1) 

            #Test by pairs
            for i_comb in xrange(len(ordered_SSs_variations[ndxSSA])): #all_possible_combinations_of_solutions:
                if( looped_results_matrix[ndxSSA,ndxSSB] >= maxNumLoops_per_sspair):#HACKtest
                    print "!!!HACKBREAK, this combination is DONE as maxNumLoops_per_sspair defines"
                    break

                if (ordered_SSs_variations[ndxSSA][i_comb][2] != 0): #If variation of the n-term on chain 1 ignore and go to next :
                    continue

                this_comb_loop_num=0
                for j_comb in xrange(len(ordered_SSs_variations[ndxSSB])):
                    if( looped_results_matrix[ndxSSA,ndxSSB] >= maxNumLoops_per_sspair):#HACKtest
                        print "!!!HACKBREAK, this combination is DONE as maxNumLoops_per_sspair defines"
                        break

                    #Temporary containers
                    this_test_pose_cp=ordered_SSs_variations[ndxSSA][i_comb][0].clone()
                    this_ispans_cp=[ordered_SSs_variations[ndxSSA][i_comb][1]]
                    this_n_extra_sizes_cp=[ordered_SSs_variations[ndxSSA][i_comb][2]]
                    this_c_extra_sizes_cp=[ordered_SSs_variations[ndxSSA][i_comb][3]]
                    this_test_pose_ss=[ordered_SSs_variations[ndxSSA][i_comb][4]]
                    if (ordered_SSs_variations[ndxSSB][j_comb][3] == 0):
                        this_test_pose_cp.append_pose_by_jump(ordered_SSs_variations[ndxSSB][j_comb][0], 1)
                        this_ispans_cp.append(ordered_SSs_variations[ndxSSB][j_comb][1])
                        this_n_extra_sizes_cp.append(ordered_SSs_variations[ndxSSB][j_comb][2])
                        this_c_extra_sizes_cp.append(ordered_SSs_variations[ndxSSB][j_comb][3])
                        this_test_pose_ss.append(ordered_SSs_variations[ndxSSB][j_comb][4])
                    else:
                        continue

                    #Variate by size is no made out of here before passing the poses to the code !:)
                    print "\t(%d of %d) Playing with size combination for ss_%02d <-> %02d: "%(((i_comb*len(ordered_SSs_variations[ndxSSB]))+j_comb),
                                                                                               (len(ordered_SSs_variations[ndxSSA])*len(ordered_SSs_variations[ndxSSB])),
                                                                                               ndxSSA+1,
                                                                                               ndxSSB+1)
                    print "\t\tN-comb:", this_n_extra_sizes_cp, "C-comb:", this_c_extra_sizes_cp
                    print "SSs of poses: ", this_test_pose_ss
                    result_looped_poses=dlooper_alg.find_closure_loops( in_pose_for_loops=this_test_pose_cp,
                                                                        pose_by_chain_ss=this_test_pose_ss,
                                                                        allowVariateSSsizeBy=0,
                                                                        only_as_input_conectivity_order_mode=True,
                                                                        b_use_global_context=b_use_global_context,
                                                                        global_context_pose_hash=global_context_pose_hash,
                                                                        num_global_context_clashes_initial=num_global_context_clashes_initial,
                                                                        max_num_loop_solutions=maxNumLoops_per_sspair_variation)

                    print result_looped_poses

                    tmp_sol_counter=0
                    for key in result_looped_poses:
                        if (len(result_looped_poses[key]) > 0):
                            tmp_pass_sol_counter=0
                            for jresult in result_looped_poses[key]:
                                tmp_pass_sol_counter+=1
                                #Store Results
                                #This might need to be checked, I am not convinced it is bugfree, might be discarding good results
                                [do_clash_withSScore,
                                 num_clashes_withSScore]=dlooper_alg.check_pm_clashes_hash_vs_pose(comparision_minSS_context_pose_hash, 
                                                                                                   jresult[1],
                                                                                                    pick_mode=pyrosetta.rosetta.core.pose.PoseCoordPickMode_N_CA_C_CB)
                                ##print "  TESTclash:", do_clash_withSScore
                                if len(num_clashes_withSScore)<=max_clashes_perLoop_withMinimalSSscaffold: #Is 1 better than 0?
                                    print "\t\tYAY! This has num loop-solutions, %d of %d"%(tmp_pass_sol_counter,
                                                                                            len(result_looped_poses[key]))
                                    looped_results_matrix[ndxSSA,ndxSSB]+=1
                                    looped_result_poses_container[ndxSSA][ndxSSB].append([jresult[1].clone(),
                                                                                          this_n_extra_sizes_cp,
                                                                                          this_c_extra_sizes_cp])
                                    this_comb_loop_num+=1
                                    tmp_sol_counter+=1
                                    sol_counter+=1

                                    #Debug pose-out
                                    if b_debug_out:
                                        tmp_filename="%s/debug/%s_testPairs_ss%02d_ss%02d_sol%06d_varN%s_varC%s_loop%03d.pdb"%(out_file_dir,
                                                                                                                            unique_run_id,
                                                                                                                            ndxSSA+1,
                                                                                                                            ndxSSB+1,
                                                                                                                            sol_counter,
                                                                                                                            this_n_extra_sizes_cp,
                                                                                                                            this_c_extra_sizes_cp,
                                                                                                                            this_comb_loop_num)
                                        tmp_filename=tmp_filename.replace(" ","_")
                                        ###print "\t\tYAY!, looped... OUT:", os.path.basename(tmp_filename)
                                        jresult[1].dump_pdb(tmp_filename)
                            print "\t\tSurviving after SS-core clash check : %d"%(tmp_sol_counter)
                    print "\t Done Playing with size combination for ss_%02d <-> %02d: )"%(ndxSSA+1,
                                                                                           ndxSSB+1)
                #assert()         

    print " Connecting by all pairs Finished!!! Yay!!! :)"
    print " The connection matrix/num-loops looks like: \n ", looped_results_matrix
    print "ALL DONE"


    #Cluster reconnection solutions per pair of SSs in the set

    #Debug Hacks:
    #clustered_loops_res_threshold=0.7
    print "Loop-msd-dasm clustering with resolution %0.2f", clustered_loops_res_threshold
    print "Searching for the best loops and clustering results, max=%d"%(keep_only_top_loops_num_max_during_clustering)
    print "IN connectivity Matrix:\n ", looped_results_matrix

    #Some containers
    shortest_looped_results_matrix=np.zeros((num_ideal_ss, num_ideal_ss), int)
    shortest_looped_results_matrix.fill(0)
    shortest_looped_result_poses_container=[]
    for i in xrange(num_ideal_ss):
        shortest_looped_result_poses_container.append([])
        for j in xrange(num_ideal_ss):
            shortest_looped_result_poses_container[i].append([])



    for irow in  xrange(len(looped_result_poses_container)):
        for jcol in xrange(len(looped_result_poses_container[irow])):
            if (len(looped_result_poses_container[irow][jcol]) > 0):
                print "Working on loop selection for comb: SS%02d<->SS%02d"%(irow+1,
                                                                             jcol+1)

                tmp_pose_container=[]
                for kndx in xrange(len(looped_result_poses_container[irow][jcol])):
                    tmp_pose_container.append(looped_result_poses_container[irow][jcol][kndx][0])

                print "Clustering loops by PIKAA msd with resolution threshold of:", clustered_loops_res_threshold
                [keep_ndxs_clustered,
                 keep_clusters_size]=k_center_clusterer_heterogeneusPoses(in_poses=tmp_pose_container,
                                                                         align_method="msd_label",
                                                                        convergence_dist=clustered_loops_res_threshold,
                                                                        pdb_info_label="DLOOPER")

                print "\tCluster center's NDXs: (num=%d):"%(len(keep_ndxs_clustered)), keep_ndxs_clustered 
                print "\tCluster's size:", keep_clusters_size
                print "\tSelecting by threshold population of : %d"%threshold_min_matches_for_clustered_loops
                clusters_better_than_threshold_ndxs=np.where(keep_clusters_size>=threshold_min_matches_for_clustered_loops)[0]
                if (len(clusters_better_than_threshold_ndxs)<1) and b_make_single_membered_cluster_if_not_friends:
                    print "!!!This connection has no solutions within the population threshold of %d"%threshold_min_matches_for_clustered_loops
                    print "!!!Since you requested b_make_single_membered_cluster_if_not_friends=True, I'll reduce the threhold until something passes"
                    tmp_threshold_min_matches_for_clustered_loops=threshold_min_matches_for_clustered_loops
                    while (len(clusters_better_than_threshold_ndxs)<1):
                        tmp_threshold_min_matches_for_clustered_loops-=1
                        print "!!!Trying threshold population of: %d"%(tmp_threshold_min_matches_for_clustered_loops)
                        clusters_better_than_threshold_ndxs=np.where(keep_clusters_size>=tmp_threshold_min_matches_for_clustered_loops)[0]
                        assert(tmp_threshold_min_matches_for_clustered_loops>0)
                keep_ndxs_clustered=keep_ndxs_clustered[clusters_better_than_threshold_ndxs]
                keep_clusters_size=keep_clusters_size[clusters_better_than_threshold_ndxs]
                print "\tPrunned cluster's centers NDXs: (num=%d):"%(len(keep_ndxs_clustered)), keep_ndxs_clustered 
                print "\tPrunend cluster's size:", keep_clusters_size
                best_keep_clustered_ndxs=[]
                print "\tKeeping a maximum of the best %d loop's clusters (as defined in keep_only_top_loops_num_max_during_clustering)"%keep_only_top_loops_num_max_during_clustering
                sorted_clusters_ndxs=np.argsort(-keep_clusters_size)[:min(len(keep_clusters_size),keep_only_top_loops_num_max_during_clustering)]
                for kcenter_ndx in sorted_clusters_ndxs:
                    best_keep_clustered_ndxs.append(kcenter_ndx)
                print "\tBest cluster's centers NDXs: (num=%d):"%(len(best_keep_clustered_ndxs)), keep_ndxs_clustered[best_keep_clustered_ndxs] 
                print "\tBest cluster's size:", keep_clusters_size[best_keep_clustered_ndxs]
                for k in best_keep_clustered_ndxs:
                    kndx=keep_ndxs_clustered[k]
                    #Append the pose to the results
                    shortest_looped_result_poses_container[irow][jcol].append(looped_result_poses_container[irow][jcol][kndx][:])
                    #Append the size of the cluster to it
                    shortest_looped_result_poses_container[irow][jcol][-1].append(keep_clusters_size[k])
                    #Append a binary count to the connectivity matrix
                    shortest_looped_results_matrix[irow][jcol]+=1
                    #HACK:  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    #ToDo: figure out loop len somehow
                    tmp_loop_len=-1
                    tmp_seq_len=-1
                    #Debug out loop
                    if b_debug_out:
                        #print "Best/shortest kndx:", irow,jcol, kndx, best_loop_len
                        tmp_filename="%s/debug/%s_testClusteredPairs_ss%02d_ss%02d_pop%d_sol%04d_ll%d_sl%d.pdb"%(out_file_dir,
                                                                                                                unique_run_id,
                                                                                                                irow+1,
                                                                                                                jcol+1,
                                                                                                                keep_clusters_size[k],
                                                                                                                kndx,
                                                                                                                tmp_loop_len,
                                                                                                                tmp_seq_len)
                        print "\t\tYAY!, best loop pair result... OUT:", tmp_filename
                        looped_result_poses_container[irow][jcol][kndx][0].dump_pdb(tmp_filename)
                        ##END test-dump PDBS

                    #Add something smart to filter by len or so
                print "Loop selection DONE for comb SS%02d<->SS%02d"%(irow+1, 
                                                                      jcol+1)


    print "New connectivity Matrix:\n ", shortest_looped_results_matrix
    print "All done"


    #Figure out the possible ways to re-connect this thing
    #only_as_input_conectivity_order_mode=True

    print "Searching connectivity solutions!"
    print "Right now the connectivity matrix passed here looks like: "
    print shortest_looped_results_matrix.astype(int)

    possible_ways_to_conect_chains=[]
    #b_find_incomplete_elements_reconnections=True
    incomplete_elements_reconnections_min_size=3
    if (not b_find_incomplete_elements_reconnections):
        num_ideal_ss_connected=len(shortest_looped_results_matrix)
        if only_as_input_conectivity_order_mode:
            print "Connectivity order has been set to AS IN THE INPUT connectivity mode, I will not try any permutations"
            is_possible_b=True
            for indx in xrange(num_ideal_ss_connected-1):
                if not shortest_looped_results_matrix[indx][indx+1]:
                    is_possible_b=False
                    break
            if is_possible_b:
                possible_ways_to_conect_chains.append(np.asarray(range(num_ideal_ss)))
        else:
            print "Connectivity order has been set to SEARCH connectivity mode, I will try all the permutations (complexity is exponential)"
            if (len(force_chain_pairs)>0):
                print "Forcing pairs:", force_chain_pairs
            for iperm in itertools.permutations(range(num_ideal_ss_connected),num_ideal_ss_connected):
                is_possible_b=True
                for jndx in xrange(len(iperm)-1):
                    if not shortest_looped_results_matrix[iperm[jndx]][iperm[jndx+1]]:
                        is_possible_b=False
                        break
                if is_possible_b:
                    is_good_c=True
                    for jforcedpair in force_chain_pairs:
                        for kndx in xrange(len(iperm)-1):
                            if (((iperm[jndx]-1 in jforcedpair) and 
                                 (iperm[jndx+1]-1 not in jforcedpair)) or ((iperm[jndx]-1 not in jforcedpair) and 
                                                                         (iperm[jndx+1]-1 in jforcedpair))):
                                is_good_c=False
                                break
                    if is_good_c:
                        possible_ways_to_conect_chains.append(np.asarray(iperm,int))
                    #else:
                    #    print "Bad:", np.asarray(iperm,int)

    elif (b_find_incomplete_elements_reconnections):
        min_elements_size=max(incomplete_elements_reconnections_min_size, 2)
        max_elements_size=len(shortest_looped_results_matrix) #min(incomplete_elements_reconnections_max_size, len(shortest_looped_results_matrix)) 
        print "I will allow incomplete reconnections as min as: "
        print "I will allow incomplete reconnections as max as: "
        for num_ideal_ss_connected in range(min_elements_size, 
                                            max_elements_size+1):
            print "TEST", num_ideal_ss_connected
            if only_as_input_conectivity_order_mode:
                print "Connectivity order has been set to AS IN THE INPUT connectivity mode, I will not try any permutations"
            else:
                print "Connectivity order has been set to SEARCH connectivity mode, I will try all the permutations (complexity is exponential)"
                print list(itertools.permutations(range(len(shortest_looped_results_matrix)), num_ideal_ss_connected))
                for iperm in itertools.permutations(range(len(shortest_looped_results_matrix)), num_ideal_ss_connected):
                    is_possible_b=True
                    for jndx in xrange(len(iperm)-1):
                        if not shortest_looped_results_matrix[iperm[jndx]][iperm[jndx+1]]:
                            is_possible_b=False
                            break
                    if is_possible_b:
                        #print "Verifying solution contains:", b_find_incomplete_elements_reconnections_must_contain_elements
                        is_allowed=True
                        for kelement in b_find_incomplete_elements_reconnections_must_contain_elements:
                            if not kelement in iperm:
                                is_allowed=False
                                break
                        if is_allowed:
                            possible_ways_to_conect_chains.append(np.asarray(iperm,int))#"""
    else:
        print "What are you doing dude???"
        assert(0==1)


    if (len(possible_ways_to_conect_chains) <= 0):
        print "I couldn't find any potential solution to reconnect all the stuff, sorry for that. I'll sadly have to end right now"
        assert(0==1)
    else:
        print "Let's see. These are the possible paths to reconnect this stuff (total= %d):"%(len(possible_ways_to_conect_chains))
        for ipath in possible_ways_to_conect_chains:
            print "**SSs->", ipath+1
        ##possible_ways_to_conect_chains


    #Reconnect it in plausible ways :)_
    #ToDo: 
    #1. Add the pathway of connection to the list of results, so it can be output in the filename
    #2. Add option to allow 2-/3-/4-/n-pieces solutions, should be as trivial as a counter of max-gaps and append_pose_by_jump 
    #   (BUT be careful at the N-ter C-ter size optimization step)
    #3. Maybe: make this much faster

    max_bb_clashes_during_pose_assembly=0

    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())
    avg_degree_filter=pyrosetta.rosetta.protocols.protein_interface_design.filters.AverageDegreeFilter()#?
    avg_degree_filter.task_factory(tf)
    avg_degree_filter.distance_threshold(8.0) #SergeySays


    print "Now I'll try to reconnect your SS-thing, hold on! Each path has a lenght of: ", len(possible_ways_to_conect_chains[0])
    solution_poses=[]
    for ipath in xrange(len(possible_ways_to_conect_chains)):
        paths=possible_ways_to_conect_chains[ipath]
        curr_poses=[]
        is_first=True
        ###print paths
        print "-Trying to reconnect path(#%d of %d):"%(ipath+1,
                                                       len(possible_ways_to_conect_chains)), "SS-seq is:", paths+1
        for jndx in xrange(len(paths)-1):
            ###print jndx
            previous_size_change=0
            if is_first:
                print "Working on first path-connection:", paths[jndx]+1,"<->", paths[jndx+1]+1
                tmp_curr_poses=[]
                for isol in xrange(len(shortest_looped_result_poses_container[paths[jndx]][paths[jndx+1]])):
                    #shortest_looped_result_poses_container  0= pose, 
                    #shortest_looped_result_poses_container  1= CtermInfo                                    
                    #shortest_looped_result_poses_container  2= NtermInfo
                    #shortest_looped_result_poses_container  3= Cluster_size
                    tmp_curr_poses.append( [shortest_looped_result_poses_container[paths[jndx]][paths[jndx+1]][isol][0].clone(),
                                        shortest_looped_result_poses_container[paths[jndx]][paths[jndx+1]][isol][1][1], #C-terminii SSn-1 size change
                                        shortest_looped_result_poses_container[paths[jndx]][paths[jndx+1]][isol][2][0], #N-terminii SSn size change
                                        shortest_looped_result_poses_container[paths[jndx]][paths[jndx+1]][isol][3],
                                        0.0, #This is a placeholder for putting stuff
                                        paths])
                                        
                is_first=False

            else:
                tmp_curr_poses=[]
                print "Working on path-connection:", paths[jndx]+1,"<->", paths[jndx+1]+1
                print "--Num of loops available for this connection: %d-- "% len(shortest_looped_result_poses_container[paths[jndx]][paths[jndx+1]])
                #ToDo1:check that things (N-C connection) are not too far at this point
                for ipose in xrange(len(curr_poses)):

                    for isol in xrange(len(shortest_looped_result_poses_container[paths[jndx]][paths[jndx+1]])):
                        ####
                        ####
                        #Copy from first pose up to where the new loop ends
                        tmp_previous_pose=pyrosetta.pose_from_sequence(fake_aminoacid_letter_for_design_loop)
                        tmp_previous_pose.copy_segment(1, curr_poses[ipose][0], 1, 1)
                        #Copy PDB-info
                        curr_pose_pdbinfo=curr_poses[ipose][0].pdb_info()
                        for jlabel in curr_pose_pdbinfo.get_reslabels(1):
                            tmp_previous_pose.pdb_info().add_reslabel(tmp_previous_pose.total_residue(), jlabel)
                        ori_secSS_size=contiguous_SS_spans[paths[jndx]][1]-contiguous_SS_spans[paths[jndx]][0]+1

                        new_end=curr_poses[ipose][0].total_residue()-ori_secSS_size
                        #Shift everything to +1 residue
                        residue_shift_ndx=0
                        if (curr_poses[ipose][1] < max_bb_clashes_during_pose_assembly): #If we removed residues from the SSn-1 N-side
                            residue_shift_ndx=-curr_poses[ipose][1]
                            new_end+=residue_shift_ndx

                        #Because the loop affects up to +2
                        new_end+=2
 
                        if(new_end > curr_poses[ipose][0].total_residue()):
                            print "W%* is wrong here??? ", residue_shift_ndx, new_end, curr_poses[ipose][0].total_residue()
                            break
                        ##print "AKSJHDF", new_end, curr_poses[ipose][0].total_residue(), ori_secSS_size

                        for ires in range( 1, new_end):
                            tmp_previous_pose.append_residue_by_bond( curr_poses[ipose][0].residue(ires+1), False);
                            #Copy PDB-info
                            for jlabel in curr_pose_pdbinfo.get_reslabels(ires+1):
                                tmp_previous_pose.pdb_info().add_reslabel(tmp_previous_pose.total_residue(), jlabel)





                        #Get a pointer to the pose that we are going to add 
                        this_pose_chunk=shortest_looped_result_poses_container[paths[jndx]][paths[jndx+1]][isol][0]
                        #Shift-ndx: The shift should be equal to the lenght of the "previous chain in the path"
                        #Shift-ndx: and minus what the previous chain grow in the C termini [1][1] and N-termini
                        #Copy from third+shift up to the last residue
                        new_chain_start_ndx=4+residue_shift_ndx
                        assert (new_chain_start_ndx<=this_pose_chunk.total_residue())
                        tmp_this_pose_minus_shift=pyrosetta.pose_from_sequence(fake_aminoacid_letter_for_design_loop)
                        tmp_this_pose_minus_shift.copy_segment(1, 
                                                               this_pose_chunk, 
                                                               1, 
                                                               new_chain_start_ndx)

                        for ires in range( new_chain_start_ndx, this_pose_chunk.total_residue()):
                            tmp_this_pose_minus_shift.append_residue_by_bond( this_pose_chunk.residue(ires+1), False);

                        #LR clash check, very important!!!
                        tmp_previous_pose_bb_hash=pyrosetta.rosetta.core.pose.xyzStripeHashPose(tmp_previous_pose,
                                                                                                pyrosetta.rosetta.core.pose.PoseCoordPickMode_N_CA_C,
                                                                                                radius=2.0)
                        b_clashing,tmp_clash_list=dlooper_alg.check_pm_clashes_hash_vs_pose(tmp_previous_pose_bb_hash, 
                                                                                            tmp_this_pose_minus_shift,
                                                                                             pick_mode=pyrosetta.rosetta.core.pose.PoseCoordPickMode_N_CA_C)
                        if (len(tmp_clash_list)>0): #This means it is clashing at the BB level
                            print "This loop/combination will clash!!!Skipping it"
                            continue


                        #Clone the previous to append the result and store it
                        tmp_result_pose=tmp_previous_pose.clone()
                        #Copy all the way from the second-res of second pose starting + the necessary shift
                        this_pose_chunk_pdbinfo=this_pose_chunk.pdb_info()
                        new_chain_start_ndx=2+residue_shift_ndx
                        for ires in range( new_chain_start_ndx, this_pose_chunk.total_residue()):
                            tmp_result_pose.append_residue_by_bond( this_pose_chunk.residue(ires+1), False);
                            #Copy PDB-info
                            for jlabel in this_pose_chunk_pdbinfo.get_reslabels(ires+1):
                                tmp_result_pose.pdb_info().add_reslabel(tmp_result_pose.total_residue(), jlabel)

                        #This check could happen waaay before when building the loops???
                        #ToDo: Move it up the way!
                        [is_good_bond,
                         nc_dist]=check_pose_bonds_are_correct(tmp_result_pose,
                                                             NC_lowBound=1.2,
                                                             NC_highBound=1.6) #1.45
                        if (not is_good_bond):
                            if (nc_dist > 3.2): #Arbitrary too-long distance to be even a long BAAAD bond
                                print "THIS N-C dist SEEMS IMPOSIBLE!!! WTF??? Something could be wrong in the code!!!"
                                print "BAD POSE AFTER NO-M!N: ", ipose, isol, nc_dist

                                debug_pose_name="%s/debug/%s_DposeFailNC_p%02d-ip%02d_is%04d_%04d_nc%0.2f"%(out_file_dir,
                                                                                                   unique_run_id,
                                                                                                 paths[jndx]+1,
                                                                                                 paths[jndx+1]+1,
                                                                                                 ipose, 
                                                                                                 isol,
                                                                                                 nc_dist)
                                print "Dump:", debug_pose_name
                                tmp_result_pose.dump_pdb(debug_pose_name+".pdb")
                                tmp_previous_pose.dump_pdb(debug_pose_name+".pdb")
                                #Result pose debug
                                debug_pose_name="%s/debug/ERROR_%s_DposeTMPr_p%02d-ip%02d_is%04d_%04d"%(out_file_dir,
                                                                                         unique_run_id,
                                                                                         paths[jndx]+1,
                                                                                         paths[jndx+1]+1,
                                                                                         ipose, 
                                                                                         isol)

                                #print "Will stop RIGHT NOW"
                                #assert(0==1)
                            print "BAD POSE AFTER NO-M!N: ", ipose, isol, nc_dist, "Silently ignoring and continuing"

                            continue

                        if True: #Add some sort of rama, score &| energy check here
                            #closure_results_dic_trimmed  0= loop_len,  -> 2
                            #closure_results_dic_trimmed  1= pose,      -> 0                                              
                            #closure_results_dic_trimmed  2= population, -> 3
                            #closure_results_dic_trimmed  3= nc_points_rmsd, -> None
                            tmp_curr_poses.append( [tmp_result_pose.clone(),
                                                   shortest_looped_result_poses_container[paths[jndx]][paths[jndx+1]][isol][1][1], 
                                                   shortest_looped_result_poses_container[paths[jndx]][paths[jndx+1]][isol][2][0],
                                                   shortest_looped_result_poses_container[paths[jndx]][paths[jndx+1]][isol][3]+curr_poses[ipose][3], #Number of loops in this cluster (accumulative)
                                                   0.0, #This will store the score and others
                                                   curr_poses[ipose][5]] )


            #Make this an array so indexing is happy
            tmp_curr_poses=np.asarray(tmp_curr_poses)

            if (len(tmp_curr_poses)==0):
                print " > OhOh, we alrready failed here! TERMINATING this solution dreams NOW. "
                curr_poses=[] #empty this so at the End doesnt' print anything weird
                break



            #DASMDASM
            #TODO: SORT USING METHOD A: SHORTEST LOOP (maybe WITH BEST AGNOSTIC SCORE?)


            #NOTE: Check this, the first loop  used to get a lot of diversity but later it doesn't happen that much.
            #Cluster the results: (does this makes sense at all)? NOTE: This may fix the previous NOTE
            #ToDo: clustering might return the most "pure" poses as centers (i.e. those with less labels, meaning less loops have been rebuilt on them!)
            tmp_curr_poses_for_clustering=[]
            for jresult in tmp_curr_poses:
                tmp_curr_poses_for_clustering.append(jresult[0].clone())
            print "Clustering loops by PIKAA msd with resolution threshold of:", clustered_loops_res_threshold
            [keep_ndxs_clustered,
             keep_clusters_size]=k_center_clusterer_heterogeneusPoses(in_poses=tmp_curr_poses_for_clustering,
                                                                     align_method="msd_label",
                                                                     convergence_dist=clustered_loops_res_threshold,
                                                                     pdb_info_label="DLOOPER")
            ##print "Keep clusters,size: ", keep_ndxs_clustered, keep_clusters_size
            #Keep only the cluster centers
            print "Keeping only the cluster centers up to a resolution of: %0.2f"% clustered_loops_res_threshold
            num_poses_before_clustering=len(tmp_curr_poses)
            tmp_curr_poses=tmp_curr_poses[keep_ndxs_clustered]
            print "OK, we reduced the redundancy by %d poses (before: %d, after: %d)"%((num_poses_before_clustering-len(tmp_curr_poses)), 
                                                               num_poses_before_clustering,
                                                               len(tmp_curr_poses))

            #UNTESTED METHOD #A minumum loop lenght
            tmp_curr_poses=np.asarray(tmp_curr_poses)
            ##tmp_pop_max=tmp_curr_poses[:,3].max()
            for krow in range(len(tmp_curr_poses)):
                test_pose=tmp_curr_poses[krow][0] #Just a pointer
                tmp_dssp_instance = pyrosetta.rosetta.core.scoring.dssp.Dssp(test_pose)
                tmp_dssp_sequence = np.ascontiguousarray(list(tmp_dssp_instance.get_dssp_reduced_IG_as_L_secstruct()),str)
                tmp_score=((tmp_dssp_sequence=='L').sum())
                tmp_curr_poses[krow][4]=tmp_score #This was: jresult[4]=tmp_score, works because is an np.array of pointers, but was not clean.
            print "Keeping only %d poses with Shortest loops (of %d possible poses)"%(max_num_loops_per_hop_when_building_topology_filtA,
                                                                                       len(tmp_curr_poses))
            #Sort by field[4] that contains the result
            curr_poses=sorted( tmp_curr_poses, 
                                key=lambda x: (x[4]) )[:min(len(tmp_curr_poses), max_num_loops_per_hop_when_building_topology_filtA)]

            print "Loop size (field[4]) sorted poses:\n"
            for krow in curr_poses:
                print krow

            #UNTESTED METHOD #6: avg_degree_filter
            tmp_curr_poses=np.asarray(tmp_curr_poses)
            ##tmp_pop_max=tmp_curr_poses[:,3].max()
            for krow in range(len(tmp_curr_poses)):
                test_pose=tmp_curr_poses[krow][0] #Just a pointer
                tmp_score=avg_degree_filter.compute(test_pose) #/test_pose.total_residue()
                tmp_curr_poses[krow][4]=tmp_score  #This was: jresult[4]=tmp_score, works because is an np.array of pointers, but was not clean.
            print "Keeping only %d poses with Best Avg Degree (of %d possible poses)"%(max_num_loops_per_hop_when_building_topology_filtB,
                                                                                       len(tmp_curr_poses))
            #Sort by field[4] that contains the result
            curr_poses=sorted( tmp_curr_poses, 
                                key=lambda x: (-x[4]) )[:min(len(tmp_curr_poses), max_num_loops_per_hop_when_building_topology_filtB)]

            print "AvgDegree (field[4]) sorted poses:\n"
            for krow in curr_poses:
                print krow

            print "  When walked by SS-positions:", paths[jndx]+1,"-", paths[jndx+1]+1, "of path:", paths+1, " and doing some mess. The number of results after trimming by size, score, etc! : ",  len(curr_poses)
            if (len(curr_poses)==0):
                print " > OhOh, we alrready failed here! TERMINATING this solution dreams NOW. "
                print "------------------PATHS with at least one solution until now: %d"% len(solution_poses)
                break


        if (len(curr_poses)>0):
            solution_poses.append(curr_poses)
            print "------------------PATHS with at least one solution until now: %d"% len(solution_poses)


    total_solutions=0
    for isol in solution_poses:
        total_solutions+=len(isol)

    if   (total_solutions < 10):
        print "Total number of results generated is: ",  total_solutions
    elif (total_solutions  < 100):            
        print "WOW! Total number of results generated is: ",  total_solutions
    elif (total_solutions  < 1000):            
        print "OMG! Total number of results generated is: ",  total_solutions
    else:            
        print "HA! The total number of results generated is... well... A lot: ",  total_solutions



    #Optimize N- C- lenghts!! (i.e. grow or make shorter)
    #ToDo: 
    ##1. MAKE THIS MUCH FASTER, possibly by changing one end at a time, opposed to optimizing both by combinations
    ##1. Add a clash check! (i.e. a growth structure might clash with a loop). Maybe some-sort of FArep could be good enough
    ### or uses hashes Like:
    ### global_context_pose_hash=pyrosetta.rosetta.core.pose.xyzStripeHashPose(global_context_pose,
    ###                                                                pyrosetta.rosetta.core.pose.PoseCoordPickMode_N_CA_C_CB, #maybe .BB instead of .HVY
    ###                                                                radius=2.0)   
    ### has_GC_clashes,GC_clash_map=dlooper_alg.check_pm_clashes_hash_vs_pose(global_context_pose_hash, 
    ###                                                                    inpose,
    ###                                                                    pick_mode=pyrosetta.rosetta.core.pose.PoseCoordPickMode_N_CA_C_CB)
    ### num_global_context_clashes_initial=len(list(GC_clash_map))

    print "Optimizing the N- C- size by maximizing the core+limbo volume"
    for isol in xrange(len(solution_poses)):
        print "\nOptimizing solutions for topology %d of %d (num=%d)"%(isol+1, len(solution_poses), len(solution_poses[isol]))
        for jndx in xrange(len(solution_poses[isol])):
            print "Optimizing solution %d_%d"%(isol,jndx)
            if (len(solution_poses[isol][jndx])==6):
                solution_poses[isol][jndx]=list(solution_poses[isol][jndx])+[None]
            else:
                print "\n\nWARNING. THIS IS EXTREMELY WEIRD! I am expecting this array (size %d now) to have size 6 at this point.\n\n"%len(solution_poses[isol][jndx])
                if (solution_poses[isol][jndx][6]!=None):
                    print "Something is in the solution slot... I am assuming it is already optimized, I'll skip to the next"

            tmp_scored_poses_array=[]
            count=0
            for kcount,ksol in enumerate(ordered_SSs_variations[solution_poses[isol][jndx][5][0]]): #First SS idealized helixes array
                ###for lsol in ordered_SSs_variations[solution_poses[isol][jndx][5][-1]]: #Last SS idealized helixes array
                #print kcount, count
                count+=1
                if( (ksol[3]!=0) ):
                    continue

                print("Profiling N-side (num_core_res,min_lenght) variation:", ksol[2:4])

                result_pose=pyrosetta.pose_from_sequence("")
                #N-term variation #add some sort of clash check here
                for lres in xrange(1,ksol[2]+1):#N-terminal size
                    result_pose.append_residue_by_bond(ksol[0].residue(lres))
                #Middle variation
                for lres in xrange(max(1,((ksol[2])*-1)+1),
                                   solution_poses[isol][jndx][0].total_residue()+1):
                    result_pose.append_residue_by_bond(solution_poses[isol][jndx][0].residue(lres))

                tmp_res_per_layer=get_residues_per_layer(in_pose=result_pose,
                                                        test_aa_large=test_aa_large,
                                                        test_aa_muta_contact=test_aa_muta_contact,
                                                        b_use_global_context=False, #b_use_global_context,
                                                        in_global_context_pose=None, #global_context_pose,
                                                        layer_sasa_difinitions=layer_sasa_difinitions,
                                                        test_aa_large_ref_sasa=test_aa_large_ref_sasa,
                                                        layer_contact_distance_difinitions=layer_contact_distance_difinitions,
                                                        b_debug=False)

                tmp_scored_poses_array.append([result_pose.clone(),
                                               len(tmp_res_per_layer["core"]),#+len(tmp_res_per_layer["limbo"]), #num
                                               result_pose.total_residue(),
                                               ksol[2], #var in N-term
                                               kcount
                                              ])


            #Sort and Keep best solution:
            #Choose best variation until now
            tmp_ksol_ndx= sorted( tmp_scored_poses_array, key=lambda x: (-x[1],x[2]) )[0][4]
            ksol=ordered_SSs_variations[solution_poses[isol][jndx][5][0]][tmp_ksol_ndx] #Last baest for N-term

            tmp_scored_poses_array_b=[]
            for lsol in ordered_SSs_variations[solution_poses[isol][jndx][5][-1]]: #Last SS idealized helixes array 
                if( (lsol[2]!=0)):
                    continue

                print("Profiling N-sideFIX C-side (num_core_res,min_lenght) variation:", ksol[2:4],lsol[2:4])
                result_pose=pyrosetta.pose_from_sequence("")
                #N-term variation #add some sort of clash check here
                for lres in xrange(1,ksol[2]+1):#N-terminal size
                    result_pose.append_residue_by_bond(ksol[0].residue(lres))
                    #Copy PDBinfolabels!!!
                    for mlabel in ksol[0].pdb_info().get_reslabels(lres):
                        result_pose.pdb_info().add_reslabel(result_pose.total_residue(), mlabel)
                #Middle variation
                for lres in xrange(max(1,((ksol[2])*-1)+1),
                                   min(solution_poses[isol][jndx][0].total_residue()+1,solution_poses[isol][jndx][0].total_residue()+lsol[3]+1)):
                    result_pose.append_residue_by_bond(solution_poses[isol][jndx][0].residue(lres))
                    #Copy PDBinfolabels!!!
                    for mlabel in solution_poses[isol][jndx][0].pdb_info().get_reslabels(lres):
                        #print mlabel
                        result_pose.pdb_info().add_reslabel(result_pose.total_residue(), mlabel)

                #C-variation #add some sort of clash check here
                if (lsol[3]>0):
                    for lres in xrange(lsol[0].total_residue()-lsol[3]+1,lsol[0].total_residue()+1):#C-terminal size
                        result_pose.append_residue_by_bond(lsol[0].residue(lres))
                        #Copy PDBinfolabels!!!
                        for mlabel in lsol[0].pdb_info().get_reslabels(lres):
                            result_pose.pdb_info().add_reslabel(result_pose.total_residue(), mlabel)

                tmp_res_per_layer=get_residues_per_layer(in_pose=result_pose,
                                                        test_aa_large=test_aa_large,
                                                        test_aa_muta_contact=test_aa_muta_contact,
                                                        b_use_global_context=False, #b_use_global_context,
                                                        in_global_context_pose=None, #global_context_pose,
                                                        layer_sasa_difinitions=layer_sasa_difinitions,
                                                        test_aa_large_ref_sasa=test_aa_large_ref_sasa,
                                                        layer_contact_distance_difinitions=layer_contact_distance_difinitions,
                                                        b_debug=False)

                tmp_scored_poses_array_b.append([result_pose.clone(),
                                               len(tmp_res_per_layer["core"]),#+len(tmp_res_per_layer["limbo"]), #num
                                               result_pose.total_residue(),
                                               ksol[2], #var in N-term
                                               lsol[3] #var in C-term
                                              ])

            #Sort and Keep best solution:
            solution_poses[isol][jndx][6]=sorted( tmp_scored_poses_array_b, key=lambda x: (-x[1],x[2]) )[0]

            #Debug out
            print "Best solution:", solution_poses[isol][jndx][6]

    print "All N- C- Optimizations DONE"


    total_solutions=0
    for isol in solution_poses:
        total_solutions+=len(isol)

    if   (total_solutions < 10):
        print "Total number of results generated is: ",  total_solutions
    elif (total_solutions  < 100):            
        print "WOW! Total number of results generated is: ",  total_solutions
    elif (total_solutions  < 1000):            
        print "OMG! Total number of results generated is: ",  total_solutions
    else:            
        print "HA! The total number of results generated is... well... A lot: ",  total_solutions


    #DebugOutPDBs
    print "TMP (non-profiled) sol num:", len(solution_poses)
    print "Generating outputs just for Fun!"
    for isol in xrange(len(solution_poses)):
        #Results are sorted already :)
        for jndx in xrange(min(len(solution_poses[isol]),max_num_results_per_topology)):
            result_pose=None
            outpdbname="%s/debug/%s_testFinalSol_%s"%(out_file_dir,
                                                   unique_run_id,
                                                    in_pose_basename[0:-4])
            if b_use_global_context:
                outpdbname+="_c%s"%in_context_pose_basename[0:-4]
                result_pose=solution_poses[isol][jndx][6][0].clone() 
                result_pose.append_pose_by_jump(global_context_pose.clone(), 1)
            else:
                result_pose=solution_poses[isol][jndx][6][0].clone()

            print "Debug out:", outpdbname+".pdb"
            result_pose.dump_scored_pdb(outpdbname+".pdb",
                                        scorefxn_vanilla_loops )

    print "Done"



    #Layer detection and assignment code (a.k.a. Rosetta is horrible and I will not use that function again if possible)
    #ToDo: Make sure that it includes the HOTSPOT a.a., as well as at least Ala/Gly if no-layer match happens
    print "Doing my version of LayerDesign, this will likely save you a great hedache and make your designs work. "
    print "Hold on! This will take some time because I am not cheap"


    solution_poses_with_layers=[]
    for isol in xrange(len(solution_poses)): #isol should read itop actually :)
        solution_poses_with_layers.append([]) 
        for jtop in xrange(min(len(solution_poses[isol]),max_num_results_per_topology)): #jtop should read jsol ^, sometimes you wonder what you were thinking at the time DASM...
            print "Working in generating layers for Topology %02d-%02d"%(isol+1,jtop+1)
            print "--Calculating layers with definitions", layer_sasa_difinitions, "(nothing is your usual subway sandwich!)"
            tmp_mut_pose_largeAA=make_polyAA_pose(solution_poses[isol][jtop][6][0],
                                                        test_aa_large)
            tmp_pose_size=tmp_mut_pose_largeAA.total_residue()
            [total_SASA, 
            per_res_sasa]=calc_pose_sasa_total_per_residue_and_atom(tmp_mut_pose_largeAA,
                                                                            probe_radii=1.7)
            by_layer_dic={}
            by_layer_dic["core"]=[]
            by_layer_dic["limbo"]=[]
            by_layer_dic["surface"]=[]
            by_layer_dic["interface"]=[]
            if b_use_global_context:
                tmp_mut_pose_largeAA_wContext=tmp_mut_pose_largeAA.clone()
                tmp_mut_pose_largeAA_wContext.append_pose_by_jump(global_context_pose.clone(), 1)
                [total_wContext_SASA, 
                per_res_sasa_wContext]=calc_pose_sasa_total_per_residue_and_atom(tmp_mut_pose_largeAA_wContext,
                                                                            probe_radii=1.7)
                tmp_delta_sasa=np.asarray(per_res_sasa[:tmp_pose_size])-np.asarray(per_res_sasa_wContext[:tmp_pose_size])
                interface_ndxs=np.where(tmp_delta_sasa>=layer_sasa_difinitions["interface"])[0]
                for lres in interface_ndxs:
                    by_layer_dic["interface"].append(lres)

            per_res_sasa=np.asarray(per_res_sasa)/test_aa_large_ref_sasa*100.0

            print "--Using pixie dust to detect core-neigbors"
            tmp_test_sasa=[]
            for lres in xrange(len(per_res_sasa)):
                    if (per_res_sasa[lres] <= layer_sasa_difinitions["core"]):
                        by_layer_dic["core"].append(lres)
                    elif (per_res_sasa[lres] <= layer_sasa_difinitions["limbo"]):
                        by_layer_dic["limbo"].append(lres)

            #WHAT???? Why we are not using the function???
            #To Do: fix calling the lib-function
            tmp_mut_pose_smallAA=make_polyAA_pose(solution_poses[isol][jtop][6][0],
                                                        test_aa_muta_contact) 
            #tmp_mut_pose_smallAA.dump_pdb("%s/test_layers/testF2.pdb"%out_file_dir)
            by_layer_dic["core"]=np.asarray(by_layer_dic["core"])
            for lres in xrange(len(per_res_sasa)):
                if (  (lres not in by_layer_dic["core"]) and 
                      (lres not in by_layer_dic["limbo"]) and
                      (lres not in by_layer_dic["interface"]) ):
                    min_dist=res_nonBB_min_distance_to( target_pose=tmp_mut_pose_smallAA,
                                                target_resA=(lres+1),
                                                target_resBs=(by_layer_dic["core"]+1))
                    if (len(min_dist)>0):
                        min_dist=np.asarray(min_dist.values()).min()
                        if (min_dist <= layer_contact_distance_difinitions["limbo"]):
                            if (per_res_sasa[lres] <layer_sasa_difinitions["surf"]):
                                by_layer_dic["limbo"].append(lres)
                            else:
                                by_layer_dic["surface"].append(lres)

                        else: #Otherwise is a surface residue
                            by_layer_dic["surface"].append(lres)
                    else: #Otherwise is a surface residue
                        by_layer_dic["surface"].append(lres)

            #out pymol selection strings:
            """
            for klayer in by_layer_dic.keys():
                tmp_line_out="select sasa_%s, resi "%klayer
                for lres in by_layer_dic[klayer]:
                    tmp_line_out=tmp_line_out+"%d+"%(lres+1)
                print tmp_line_out[:-1]+";"
            """

            print "--Now merging this into usable resfiles!"
            #Finally define the residues that can be in the layers
            tmp_test_pose=solution_poses[isol][jtop][6][0].clone()
            can_mutate_positions=[]
            for klayer in by_layer_dic.keys():
                for lres in by_layer_dic[klayer]:
                    
                    tmp_reslabel_aas="PIKAA_"
                    if (tmp_test_pose.pdb_info().res_haslabel(lres+1, pdb_info_label_keyword)):
                        tmp_reslabel_aas="PIKAA_%s"%tmp_test_pose.residue(lres+1).name1()
                    #Keep track of non-hotspot positions
                    else:
                        can_mutate_positions.append(lres+1)
                    
                    tmp_thisres_reslabels=tmp_test_pose.pdb_info().get_reslabels(lres+1)
                    has_previous_label=False
                    for jlabel in tmp_thisres_reslabels:
                        def_aa=[]
                        if (jlabel[0:6]=="PIKAA_"):
                            has_previous_label=True
                            #print tmp_test_pose.pdb_info().get_reslabels(lres+1)
                            def_aa=list(jlabel[6:])
                            common_aa=list(set(def_aa).intersection(layer_aa_dic[klayer]))
                            if (len(common_aa) <= 0):
                                print "WARNING: I don't know what to do with a.a. %d, since the layers doesn't match the PIKKA predefined, Hence I'll not touch it, sorry!"%(lres+1)
                                tmp_reslabel_aas=jlabel
                            else:
                                for jaa in common_aa:
                                    tmp_reslabel_aas=tmp_reslabel_aas+jaa
                    if not has_previous_label:
                        for jaa in layer_aa_dic[klayer]:
                            tmp_reslabel_aas=tmp_reslabel_aas+jaa
                    new_labels=[]
                    for jlabel in tmp_thisres_reslabels:
                        if (jlabel[0:5]!="PIKAA"):
                            new_labels.append(jlabel)
                    new_labels.append(klayer)
                    #new_labels.append("PIKAA")
                    new_labels.append(tmp_reslabel_aas)
                    #Write labels to the out PDB
                    tmp_test_pose.pdb_info().clear_reslabel(lres+1)
                    for mlabel in new_labels:
                        tmp_test_pose.pdb_info().add_reslabel(lres+1, mlabel)

            #Separate PIKAA into multiple PIKAA_labels for ROSETTASCRIPTS
            for lres in xrange(1,tmp_test_pose.total_residue()+1):
                tmp_thisres_reslabels=tmp_test_pose.pdb_info().get_reslabels(lres)
                #print(tmp_thisres_reslabels)
                final_labels=[]
                this_allowed_aas=[]
                for jlabel in tmp_thisres_reslabels:
                    if (jlabel[0:6]=="PIKAA_"):
                        for kaa_name in jlabel[6:]:
                            this_allowed_aas.append(kaa_name)
                    else:
                        final_labels.append(jlabel) #append right the way the non-PIKAA labels

                #Uggly Hack for adding the inverse labels, it should be 
                for kaa in layer_aa_dic["all"]: #This are all the possible aminoacids!!!
                    if kaa not in this_allowed_aas: #If not in the identified ones, then just forbidit
                        final_labels.append("notPIKAA_"+kaa)
                    else:
                        final_labels.append("PIKAA_"+kaa)

                #Replace labels with new labels
                tmp_test_pose.pdb_info().clear_reslabel(lres)
                for jlabel in final_labels:
                    tmp_test_pose.pdb_info().add_reslabel(lres, jlabel)

            if b_generate_final_output_pose_mutating_all_non_hotspots_to_single_absurd_aa:
                tmp_test_pose=mutate_residues_without_neighbors_consideration( inpose=tmp_test_pose , 
                                                                               mutant_positions=can_mutate_positions , 
                                                                               mutant_aas=(aa_for_mutating_all_non_hotspots_in_absurd_output*len(can_mutate_positions)) )


            #print tmp_test_pose_with_out_aa
            solution_poses_with_layers[-1].append([tmp_test_pose.clone(),
                                               solution_poses[isol][jtop][1],
                                               solution_poses[isol][jtop][2],
                                               solution_poses[isol][jtop][3],
                                               solution_poses[isol][jtop][4]])
            print "-Done and to the next"
    print "ACTUALLY...ALL DONE!"


    #Finally out PDBs for design
    if (len(solution_poses_with_layers) > 0):
        print "OK, so you got solutions for this problem! (Not waranted to work in the real world indeed, I can't do that because you might sue me if unhappy)"
        print "\nNow, I'll write the poses to the Disk (this might be snail sloooowww in the DIGs, not my fault indeed!)\n"
        print "The output path is: \n\t ", out_file_dir 
        print "NOTE: I'll output a maximum of %d structures per topology (as specified in max_num_results_per_topolog)"%max_num_results_per_topology
        for isol in xrange(len(solution_poses_with_layers)):
                #Results are sorted already :)
                for jndx in xrange(min(len(solution_poses_with_layers[isol]),max_num_results_per_topology)):
                    result_pose=None
                    outpdbname="%s/%s_%s_dlooperf_s%04d%04d"%(out_file_dir,
                                                                unique_run_id,
                                                                in_pose_basename[0:-4],
                                                                isol,
                                                                jndx,)
                    print "Generating output to:", outpdbname
                    if b_use_global_context:
                        outpdbname+="_c%s"%in_context_pose_basename[0:-4]
                        result_pose=solution_poses_with_layers[isol][jndx][0].clone()
                        result_pose.append_pose_by_jump(global_context_pose.clone(), 1)
                    else:
                        result_pose=solution_poses_with_layers[isol][jndx][0].clone()
                    #ToDo1:add a constrained minimization here??? Need to do Val if you want to do this.
                    if b_generate_final_output_pose_mutating_all_non_hotspots_to_single_absurd_aa:
                        result_pose.dump_pdb( "%s_poly%s.pdb"%(outpdbname,
                                                               test_aa_for_output))
                    else:
                        result_pose.dump_pdb( "%s.pdb"%(outpdbname))

        print "OK, I am done. Have a happy session of residue design, packing and camping!"
        print "\nEnd of loop reconnection, run design now using the HOTSPOT labels restrictions? See ya soon. DASM.\n"
    else:
        print "We have failed so badly reconnecting your SSthing. Mhhh... Are you sure that you gave me something reasonable as input?"


    print "That was all... Bye"
    #END END END . . .


if __name__ == '__main__':
     main()

