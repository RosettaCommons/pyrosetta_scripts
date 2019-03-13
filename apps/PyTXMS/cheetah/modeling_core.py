#!/usr/bin/env python

from rosetta import *
from pyrosetta import *
from rosetta.protocols.rigid import *
from PyDocking import PyDocking
from rosettaxlv import rosettaxlv
from rosettaxl import rosettaxl
from string import digits
from math import sqrt, pow, exp
import glob
import os
import argparse

def modeling_core(input_pdb, partners, top_XL_file, cut_off, num_of_models, num_of_top_filters):

    init()
    pose = Pose()
    pose_from_file(pose, input_pdb)
    score_t_list = []
    score_normal_list = []
    pose_list = []
    out_XLs_list = []
    top_score = 0
    top_pose = Pose()
    
    model_dir = "LR_models"
    if not os.path.exists(model_dir):
    	os.makedirs(model_dir)
    job_output = "LR_models/lr_model"
    
    ## This counter will stop the calculation if after 10 rounds (can be user specified) it can't 
    ## find a new pose justifying more XLs. Then we report the pose with MAX number of XLs. 
    break_counter = 1
    # num_of_models = 10 #for test
    # num_of_top_filters = 5 #for test

    #top_XL_file_name = 'top_model_XLs.txt'
    #top_model_XLs_file = open(top_XL_file_name, 'w')

    ## Global docking (low resolution)
    while break_counter <= num_of_models:
        print "Modeling is started ...\n"
        print "Round number ", break_counter
        XL_below_cutoff = []
        temp_pose = Pose()
        temp_pose.assign(pose)
        dock_pose = PyDocking(temp_pose, partners, 3.0, 8.0, 1, job_output, True, False)
        
        if dock_pose is not None:
            
            score_t, score_normal, XL_below_cutoff = rosettaxlv(dock_pose, top_XL_file, cut_off)

            ## keeping top structures according to the "num_of_top_filters"
            if len(score_t_list) > num_of_top_filters:
                min_score = min(score_t_list)
                index_min = score_t_list.index(min_score)
                if score_t > min_score:
                	score_t_list[index_min] = score_t
                	score_normal_list[index_min] = score_normal
                	pose_list[index_min] = dock_pose
            
            else:
                score_t_list.append(score_t)
                score_normal_list.append(score_normal)
                pose_list.append(dock_pose)
            
            break_counter += 1
                

    print "Modeling is finished!\n"
    # print "score of top models: ", score_t_list, score_normal_list
    ## storing all the top models
    for num_pos, struct in enumerate(pose_list):
        pdb_name = job_output+"_"+str(num_pos)
    	out_name = job_output+"_"+str(num_pos)+"_XLs.txt"
        out_XLs_list.append(out_name)
    	scorefxn_low = create_score_function('interchain_cen')
    	jd = PyJobDistributor(pdb_name, 1, scorefxn_low)
        struct.pdb_info().name(pdb_name + '_' + str( num_pos ) + '_fa')
    	jd.output_decoy(struct)
        rosettaxl(struct, partners, cut_off, out_name)

    return out_XLs_list
