#!/usr/bin/env python

# from __future__ import print_function

from rosetta import *
from pyrosetta import *
from rosetta.protocols.rigid import *

from string import digits
from math import sqrt, pow, exp

# init(extra_options = "-constant_seed")

## rosettaxlv parse the input pdb (according to the partners involved)
## to obtain number of XLs it fulfills.
## inputs: pdb, list of XLs, partners, dist cut-off

def rosetta_eu_dist(pose, a, b):
	dist = sqrt(pow(( pose.residue( b ).atom( "CB" ).xyz().x - pose.residue( a ).atom( "CB" ).xyz().x ) , 2) \
		+ pow(( pose.residue( b ).atom( "CB" ).xyz().y - pose.residue( a ).atom( "CB" ).xyz().y ) , 2) \
		+ pow(( pose.residue( b ).atom( "CB" ).xyz().z - pose.residue( a ).atom( "CB" ).xyz().z ) , 2))
	return dist


## 0. initializing

def rosettaxlv(pose, top_XL_file, cut_off):

    ## 1. Reading and storing the pdb in a pose
    # pose = Pose()
    # pose_from_file(pose, pdb)
    sequence = pose.sequence()
    pdb_info = pose.pdb_info()


    # print (str( pose.residue( 20 ).type().base_name() ))
    # print (str( pose.residue( 100 ).type().base_name() ))
    print "\nThe protein sequence is: \n", sequence, "\n"


    ## 2. Reading the xl_file
    with open(top_XL_file,'r') as xlfile:
        top_XL = xlfile.read().splitlines() # each rows as one element of the list

    ## output file ready to be written
    # out_sc_file = open('xl_based_score.txt', 'w')
    # out_sc_file.write("PDB, list_xl, score_t, n_dist_score\n")
    # out_sc_file.write(pdb+", ")

    output_xl_number = 0
    normal_dist_score = 0.0
    good_XL_list = []
    for num_xl,xl in enumerate(top_XL):
        print num_xl+1, xl

        K_pos_P1 = 0
        K_pos_P2 = 0
	eulidean_dist = 10000.0

	xl_without_digit = xl.translate(None, digits)
	peptide1 = xl_without_digit.split('--')[0].replace("-","").replace(".","")[:-2]
	peptide2 = xl_without_digit.split('--')[1].replace("-","")[:-3]

	K_pos_P1 = peptide1.find('K') + 1
	K_pos_P2 = peptide2.find('K') + 1

	multiple_occ_list1 = [x for x in xrange(len(sequence)) if sequence.find(peptide1, x) == x]
	multiple_occ_list2 = [y for y in xrange(len(sequence)) if sequence.find(peptide2, y) == y]

	## finding minimum distance if multiple occurance happened
	seq_pos_p1_k = 0
	seq_pos_p2_k = 0
	tmp_dist = 10000.0
	for pos1 in multiple_occ_list1:
	    for pos2 in multiple_occ_list2:
		if tmp_dist > rosetta_eu_dist(pose, pos1+K_pos_P1, pos2+K_pos_P2):
		    tmp_dist = rosetta_eu_dist(pose, pos1+K_pos_P1, pos2+K_pos_P2)
		    seq_pos_p1_k = pos1+K_pos_P1
		    seq_pos_p2_k = pos2+K_pos_P2
				 
	print "aa positions on the sequence: ", seq_pos_p1_k, seq_pos_p2_k

	if ((peptide1 in sequence) and (peptide2 in sequence)):
	    eulidean_dist = rosetta_eu_dist(pose, seq_pos_p1_k, seq_pos_p2_k)

	    print "Euclidean distance is:  ", eulidean_dist, "\n"
		
	else:
	    print "The XL is not found on the protein sequence. Check each peptide to be valid!\n"
	    # out_sc_file.write("XL is not found\n")

		
	if eulidean_dist <= cut_off:
            good_XL_list.append(xl)
	    output_xl_number += 1
	    normal_dist_score += (10 * (1/(exp ( pow( (eulidean_dist-17.5),2 )/30 )))) #Normal Distribution Formula Ae^-(x-M)/2q^2
	    # out_sc_file.write(xl+", ")
	    print "found a valid XL!", eulidean_dist, "\n"
			
    # out_sc_file.write(str(output_xl_number)+", "+str(normal_dist_score))
    return output_xl_number, normal_dist_score, good_XL_list



