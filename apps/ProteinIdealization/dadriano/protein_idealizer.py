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

####
#ToDo add diferent n-mer len

#Some general python imports
#import json
#import pickle
import pandas 
import itertools
import numpy as np
import os as os
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from collections import OrderedDict
	
	
def licence():
 print "   ***THIS IS AN INTERNAL DEVELOPER VERSION. DO NOT REDISTRIBUTE***\n\
   DASM Protein Idealizer (not an official name)\n\
   This is a novel method that uses Statistical Clustered Information from the vall to idealize SS and build de-novo loops\n\
   Author: Daniel-Adriano Silva\n\
   Ver: internal release v0.1b\n\
   Date:  Nov/6/2014\n\
   Intended use: Fold-It, general profiling of protein's design quality\n\
   Citation: To the date this is still unpublished work, if you inted to citate contact the authors.\n\
   Contact: Daniel-Adriano Silva <dadriano@gmail.com> or  David Baker <dabaker@u.washington.edu>\n\
     DASM Protein Idealizer (not an official name)\n\
     Version: internal release v0.1b, Nov/6/2014\n\
     This is a novel method that uses Statistical Clustered Information from the vall to idealize SS and build de-novo loops\n\
     Authors: Daniel-Adriano Silva, Enrique Marcos and David Baker. 2014. Rosetta Commons & The Baker lab.\n\
     Questions: dadriano@gmail.com\n\
     Documentation: Does not exist yet\n\
     Copyright: The Baker Lab & Rosetta Commons (2014)\n\
   ***THIS IS AN INTERNAL DEVELOPER VERSION. DO NOT REDISTRIBUTE***"
	
	
#plot in 3D
def plot_3D(XYZmat ):
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	x=XYZmat.T[0]
	y=XYZmat.T[1]
	z=XYZmat.T[2]
	ax.scatter(x, y, -z, zdir='z', c= 'red')
	ax.set_xlabel("X")
	ax.set_ylabel("Y")
	ax.set_zlabel("Z")
	#ax.set_xlim(-5,5)
	#ax.set_ylim(-5,5)
	#ax.set_zlim(-5,5)
	ax.view_init(21, -21)
	
	
#Take a pose convert it to a np array (look: "N","CA","C" ordering to match as Alex specification in the database(s))
def pose_endp_to_array(pose, ndxs, atoms=["CA","N","C"]):
	ep_coors = (np.zeros(((len(atoms)*2),3),float))
	#Extract the coordinates
	tmp_index=0
	for resi in (ndxs):
		for atom in (atoms):
			for dim in range(0,3): 
				ep_coors[tmp_index][dim] = pose.residue(resi).xyz(atom)[dim]
			tmp_index+=1
	return ep_coors
	
	
def pose_res_to_array(pose, ndxs, atoms=["CA","N","C"]):
	coors = np.zeros(( (len(atoms)*(ndxs[1]-ndxs[0]+1)),3),float)
	#Extract the coordinates
	tmp_index=0
	for resi in range (ndxs[0], ndxs[1]+1):
		for atom in (atoms):
			for dim in range(0,3): 
				coors[tmp_index][dim] = pose.residue(resi).xyz(atom)[dim]
			tmp_index+=1
	return coors
	
	
def rmsd_2_np_arrays(crds1, crds2):
	"""Returns RMSD between 2 sets of [nx3] numpy array"""
	#D assert(crds1.shape[1] == 3)
	#D assert(crds1.shape == crds2.shape)
	##Corrected to account for removal of the COM
	COM1 = np.sum(crds1,axis=0) / crds1.shape[0]
	COM2 = np.sum(crds2,axis=0) / crds2.shape[0]
	crds1-=COM1
	crds2-=COM2
	n_vec = np.shape(crds1)[0]
	correlation_matrix = np.dot(np.transpose(crds1), crds2)
	v, s, w_tr = np.linalg.svd(correlation_matrix)
	is_reflection = (np.linalg.det(v) * np.linalg.det(w_tr)) < 0.0
	if is_reflection:
			s[-1] = - s[-1]
			v[:,-1] = -v[:,-1]
	E0 = sum(sum(crds1 * crds1)) + sum(sum(crds2 * crds2))
	rmsd_sq = (E0 - 2.0*sum(s)) / float(n_vec)
	rmsd_sq = max([rmsd_sq, 0.0])
	return np.sqrt(rmsd_sq)
	
def rmsd_2_np_arrays_wRotationMatrix(crds1, crds2):
	"""Returns RMSD between 2 sets of [nx3] numpy array"""
	#D assert(crds1.shape[1] == 3)
	#D assert(crds1.shape == crds2.shape)
	##Corrected to account for removal of the COM
	COM1 = np.sum(crds1,axis=0) / crds1.shape[0]
	COM2 = np.sum(crds2,axis=0) / crds2.shape[0]
	crds1-=COM1
	crds2-=COM2
	n_vec = np.shape(crds1)[0]
	correlation_matrix = np.dot(np.transpose(crds1), crds2)
	v, s, w_tr = np.linalg.svd(correlation_matrix)
	is_reflection = (np.linalg.det(v) * np.linalg.det(w_tr)) < 0.0
	if is_reflection:
			s[-1] = - s[-1]
			v[:,-1] = -v[:,-1]
	E0 = sum(sum(crds1 * crds1)) + sum(sum(crds2 * crds2))
	rmsd_sq = (E0 - 2.0*sum(s)) / float(n_vec)
	rmsd_sq = max([rmsd_sq, 0.0])
	#Calculate the rotation matrix
	rMtx=np.dot(v, w_tr)
	#Calculate the translation Vector
	tVec=COM1-(np.dot(COM2, np.linalg.inv(rMtx)))
	###Are you kidding me??? Is this the correct way to build arrays inside 
	#of Rosetta core? No indexing?? Who was the creator of this???
	###rMtx_xyzM=rosetta.numeric.xyzMatrix_double()
	###rMtx_xyzM.xx=rMtx[0,0]
	###rMtx_xyzM.xy=rMtx[0,1]
	###rMtx_xyzM.xz=rMtx[0,2]
	###rMtx_xyzM.yx=rMtx[1,0]
	###rMtx_xyzM.yy=rMtx[1,1]
	###rMtx_xyzM.yz=rMtx[1,2]
	###rMtx_xyzM.zx=rMtx[2,0]
	###rMtx_xyzM.zy=rMtx[2,1]
	###rMtx_xyzM.zz=rMtx[2,2]
	###tVec_xyzV=rosetta.numeric.xyzVector_double()
	###tVec_xyzV.x=tVec[0]
	###tVec_xyzV.y=tVec[1]
	###tVec_xyzV.z=tVec[2]
	return np.sqrt(rmsd_sq), rMtx, tVec #rMtx_xyzM, tVec_xyzV
	
	
def rmsd_bb(pose1, pose2):
	numRes=pose1.total_residue()
	#D assert (numRes == pose2.total_residue())
	coorA=np.zeros(((4*numRes),3), float)
	coorB=np.zeros(((4*numRes),3), float)
	counter=0
	for res in range (1, (numRes+1)):
		for dim in range(0,3): 
			coorA[counter,dim]=(pose1.residue(res).xyz("CA")[dim])
			coorB[counter,dim]=(pose2.residue(res).xyz("CA")[dim])
		counter+=1
		for dim in range(0,3): 
			coorA[counter,dim]=(pose1.residue(res).xyz("C")[dim])
			coorB[counter,dim]=(pose2.residue(res).xyz("C")[dim])
		counter+=1
		for dim in range(0,3): 
			coorA[counter,dim]=(pose1.residue(res).xyz("O")[dim])
			coorB[counter,dim]=(pose2.residue(res).xyz("O")[dim])
		counter+=1
		for dim in range(0,3): 
			coorA[counter,dim]=(pose1.residue(res).xyz("N")[dim])
			coorB[counter,dim]=(pose2.residue(res).xyz("N")[dim])
		counter+=1
	rmsdVal = rmsd_2_np_arrays(coorA, coorB)
	return rmsdVal
	
	

	
def msd_bb(pose1, pose2):
	numRes=pose1.total_residue()
	#D assert (numRes == pose2.total_residue())
	coorA=np.zeros(((4*numRes),3), float)
	coorB=np.zeros(((4*numRes),3), float)
	counter=0
	for res in range (1, (numRes+1)):
		for dim in range(0,3): 
			coorA[counter,dim]=(pose1.residue(res).xyz("CA")[dim])
			coorB[counter,dim]=(pose2.residue(res).xyz("CA")[dim])
		counter+=1
		for dim in range(0,3): 
			coorA[counter,dim]=(pose1.residue(res).xyz("C")[dim])
			coorB[counter,dim]=(pose2.residue(res).xyz("C")[dim])
		counter+=1
		for dim in range(0,3): 
			coorA[counter,dim]=(pose1.residue(res).xyz("O")[dim])
			coorB[counter,dim]=(pose2.residue(res).xyz("O")[dim])
		counter+=1
		for dim in range(0,3): 
			coorA[counter,dim]=(pose1.residue(res).xyz("N")[dim])
			coorB[counter,dim]=(pose2.residue(res).xyz("N")[dim])
		counter+=1
	msdVal=0.0
	for i in range(len(coorA)):
		tmpVal_buff= np.linalg.norm(coorA[i] - coorB[i])
		msdVal+=(tmpVal_buff*tmpVal_buff)
	msdVal=msdVal/len(coorA)
	msdVal=np.sqrt(msdVal)
	return msdVal
	
	
def msd_bb_by_ndxs(pose1, init_res1, end_res1, pose2, init_res2, end_res2):
	numRes=(end_res1-init_res1+1)
	coorA=np.zeros(((4*numRes),3), float)
	coorB=np.zeros(((4*numRes),3), float)
	
	counter=0
	for res in range (init_res1, (end_res1+1)):
		for dim in range(0,3): 
			coorA[counter,dim]=(pose1.residue(res).xyz("CA")[dim])
		counter+=1
		for dim in range(0,3): 
			coorA[counter,dim]=(pose1.residue(res).xyz("C")[dim])
		counter+=1
		for dim in range(0,3): 
			coorA[counter,dim]=(pose1.residue(res).xyz("O")[dim])
		counter+=1
		for dim in range(0,3): 
			coorA[counter,dim]=(pose1.residue(res).xyz("N")[dim])
		counter+=1
	counter=0
	for res in range (init_res2, (end_res2+1)):
		for dim in range(0,3): 
			coorB[counter,dim]=(pose2.residue(res).xyz("CA")[dim])
		counter+=1
		for dim in range(0,3): 
			coorB[counter,dim]=(pose2.residue(res).xyz("C")[dim])
		counter+=1
		for dim in range(0,3): 
			coorB[counter,dim]=(pose2.residue(res).xyz("O")[dim])
		counter+=1
		for dim in range(0,3): 
			coorB[counter,dim]=(pose2.residue(res).xyz("N")[dim])
		counter+=1
		
	#Calculate the MSD
	msdVal=0.0
	for i in range(len(coorA)):
		tmpVal_buff= np.linalg.norm(coorA[i] - coorB[i])
		msdVal+=(tmpVal_buff*tmpVal_buff)
	msdVal=msdVal/len(coorA)
	msdVal=np.sqrt(msdVal)
	return msdVal
	
	
def rmsd_Ca(pose1, pose2):
	numRes=pose1.total_residue()
	#D assert (numRes == pose2.total_residue())
	coorA=np.zeros(((1*numRes),3), float)
	coorB=np.zeros(((1*numRes),3), float)
	counter=0
	for res in range (1, (numRes+1)):
		for dim in range(0,3): 
			coorA[counter,dim]=(pose1.residue(res).xyz("CA")[dim])
			coorB[counter,dim]=(pose2.residue(res).xyz("CA")[dim])
		counter+=1
	rmsdVal = rmsd_2_np_arrays(coorA, coorB)
	return rmsdVal
	
	
def distance_between_two_Ca(pose1, resnum1, pose2, resnum2):
	dist=np.zeros(3,float)
	tmpvec=pose1.residue(resnum1).xyz("CA") - pose2.residue(resnum2).xyz("CA")
	dist[0] = tmpvec[0]
	dist[1] = tmpvec[1]
	dist[2] = tmpvec[2]
	return np.linalg.norm(dist)
	
	
def generate_fragment_combination_list_by_TCM( start_ndxs, 
											  this_depth, 
											  max_depth, 
											  TCM_ndx_p, 
											  TCM_ndx_p_len,
											  combinations_buffer, 
											  combinations_array ):
	if(this_depth >= max_depth):
	   # if(len(combinations_buffer) == max_depth):
		combinations_array.append(combinations_buffer[:])
		return
	else:
		for i in start_ndxs:
			if(TCM_ndx_p_len[i] > 0):
				curr_combination=combinations_buffer[:]
				curr_combination.append(i)
				generate_fragment_combination_list_by_TCM( TCM_ndx_p[i], (this_depth+1), 
													  max_depth, TCM_ndx_p, TCM_ndx_p_len, curr_combination, combinations_array )
	return
	
	

#Determine the average probability for each position/a.a. and print a PSSM 
def create_fragment_based_pssm(self, aaProbabilies, outfilename="default.pssm"):
	aaProbabilities_conjunct=np.zeros((len(aaProbabilies),20),float)
	pert_pseudocount=0.000001
	for i in range( len(aaProbabilies)):
		#print aaProbabilies[i]
		for j in range(len(aaProbabilies[i])):
			aaProbabilities_conjunct[i]+=(aaProbabilies[i][j]+pert_pseudocount)
		aaProbabilities_conjunct[i]=(aaProbabilities_conjunct[i]/len(aaProbabilies[i]))
	
	
	
	"""aaA_1 aaR_1 aaN_1 aaD_1 aaC_1 aaE_1 aaQ_1 aaG_1 aaH_1 aaI_1 aaL_1 aaK_1 aaM_1 aaF_1 aaP_1 aaS_1 aaT_1 aaW_1 aaY_1 aaV_1"""
	
	"""
	FavorSequenceProfile Mover doc says:
	You can set how to scale the given values with the "scaling" settings. 
	The default value of "prob" does a per-residue Boltzmann-weighted probability 
	based on the profile score (the unweighted scores for all 20 amino acid identities 
	at any given position sum to -1.0). A setting of "global" does a global linear fixed-zero 
	rescaling such that all (pre-weighted) values fall in the range of -1.0 to 1.0. 
	A setting of "none" does no adjustment of values.
	
	NOTE: Rosetta documentation is WRONG, it seems the input probability is either:
	a) from 0.0 to 1.0
	or, most probably:
	b) from -something to +something, where is? the weighted Boltzman probability of 0.0 to 1.0.
	For example: log(P/rPaa,2)
	"""
	"""
	#Background probabilities of each a.a.
	#BLAST_LetterProb Robinson_prob (Thanks Sergey)
	'A',	78.05
	'R',	51.29
	'N',	44.87
	'D',	53.64
	'C',	19.25
	'E',	62.95
	'Q',	42.64
	'G',	73.77
	'H',	21.99
	'I',	51.42
	'L',	90.19
	'K',	57.44
	'M',	22.43
	'F',	38.56
	'P',	52.03
	'S',	71.2
	'T',	58.41
	'W',	13.3
	'Y',	32.16
	'V',	64.41
	"""     
	
	for aa_position in aaProbabilities_conjunct:
		aa_position[0]=np.log2((aa_position[0]/0.07805))
		aa_position[1]=np.log2((aa_position[1]/0.05129))
		aa_position[2]=np.log2((aa_position[2]/0.04487))
		aa_position[3]=np.log2((aa_position[3]/0.05364))
		aa_position[4]=np.log2((aa_position[4]/0.01925))
		aa_position[5]=np.log2((aa_position[5]/0.06295))
		aa_position[6]=np.log2((aa_position[6]/0.04264))
		aa_position[7]=np.log2((aa_position[7]/0.07377))
		aa_position[8]=np.log2((aa_position[8]/0.02199))
		aa_position[9]=np.log2((aa_position[9]/0.05142))
		aa_position[10]=np.log2((aa_position[10]/0.09019))
		aa_position[11]=np.log2((aa_position[11]/0.05744))
		aa_position[12]=np.log2((aa_position[12]/0.02243))
		aa_position[13]=np.log2((aa_position[13]/0.03856))
		aa_position[14]=np.log2((aa_position[14]/0.05203))
		aa_position[15]=np.log2((aa_position[15]/0.07120))
		aa_position[16]=np.log2((aa_position[16]/0.05841))
		aa_position[17]=np.log2((aa_position[17]/0.01330))
		aa_position[18]=np.log2((aa_position[18]/0.03216))
		aa_position[19]=np.log2((aa_position[19]/0.06441))
	
		
	aa_profile=""
	for aas_prob in aaProbabilities_conjunct:
		if (aas_prob.argmax() == 0):
			aa_profile+="A"
		elif(aas_prob.argmax() == 1):
			aa_profile+="R"
		elif(aas_prob.argmax() == 2):
			aa_profile+="N"
		elif(aas_prob.argmax() == 3):
			aa_profile+="D"
		elif(aas_prob.argmax() == 4):
			aa_profile+="C"
		elif(aas_prob.argmax() == 6):
			aa_profile+="Q"
		elif(aas_prob.argmax() == 5):
			aa_profile+="E"
		elif(aas_prob.argmax() == 7):
			aa_profile+="G"
		elif(aas_prob.argmax() == 8):
			aa_profile+="H"
		elif(aas_prob.argmax() == 9):
			aa_profile+="I"
		elif(aas_prob.argmax() == 10):
			aa_profile+="L"
		elif(aas_prob.argmax() == 11):
			aa_profile+="K"
		elif(aas_prob.argmax() == 12):
			aa_profile+="M"
		elif(aas_prob.argmax() == 13):
			aa_profile+="F"
		elif(aas_prob.argmax() == 14):
			aa_profile+="P"
		elif(aas_prob.argmax() == 15):
			aa_profile+="S"
		elif(aas_prob.argmax() == 16):
			aa_profile+="T"
		elif(aas_prob.argmax() == 17):
			aa_profile+="W"
		elif(aas_prob.argmax() == 18):
			aa_profile+="Y"
		elif(aas_prob.argmax() == 19):
			aa_profile+="V"
		else:
			print "ERROR!!! That aa does not exist?!"
			assert(0==1)
			
		
	#Be careful E Q are reversed in my database (Double check!)
	lineout_aaPSSM="PSSM from Fragment Protein Idealizer\nPSSM from Fragment Protein Idealizer\n"
	lineout_aaPSSM+="          A     R     N     D     C     Q     E     \
	G     H     I     L     K     M     F     P     S     T     W     Y     V\n"
	for i in range(len(aaProbabilities_conjunct)):
		iRowMax=aaProbabilities_conjunct[i].max()
		lineout_aaPSSM+= "%4d %s  " % ((i+1), aa_profile[i])
		lineout_aaPSSM+= " %0.2f" % (aaProbabilities_conjunct[i][0])
		lineout_aaPSSM+= " %0.2f" % (aaProbabilities_conjunct[i][1])
		lineout_aaPSSM+= " %0.2f" % (aaProbabilities_conjunct[i][2])
		lineout_aaPSSM+= " %0.2f" % (aaProbabilities_conjunct[i][3])
		lineout_aaPSSM+= " %0.2f" % (aaProbabilities_conjunct[i][4])
		lineout_aaPSSM+= " %0.2f" % (aaProbabilities_conjunct[i][6])
		lineout_aaPSSM+= " %0.2f" % (aaProbabilities_conjunct[i][5])
		lineout_aaPSSM+= " %0.2f" % (aaProbabilities_conjunct[i][7])
		lineout_aaPSSM+= " %0.2f" % (aaProbabilities_conjunct[i][8])
		lineout_aaPSSM+= " %0.2f" % (aaProbabilities_conjunct[i][9])
		lineout_aaPSSM+= " %0.2f" % (aaProbabilities_conjunct[i][10])
		lineout_aaPSSM+= " %0.2f" % (aaProbabilities_conjunct[i][11])
		lineout_aaPSSM+= " %0.2f" % (aaProbabilities_conjunct[i][12])
		lineout_aaPSSM+= " %0.2f" % (aaProbabilities_conjunct[i][13])
		lineout_aaPSSM+= " %0.2f" % (aaProbabilities_conjunct[i][14])
		lineout_aaPSSM+= " %0.2f" % (aaProbabilities_conjunct[i][15])
		lineout_aaPSSM+= " %0.2f" % (aaProbabilities_conjunct[i][16])
		lineout_aaPSSM+= " %0.2f" % (aaProbabilities_conjunct[i][17])
		lineout_aaPSSM+= " %0.2f" % (aaProbabilities_conjunct[i][18])
		lineout_aaPSSM+= " %0.2f" % (aaProbabilities_conjunct[i][19])
		lineout_aaPSSM+= "\n"
	f = open(outfilename,'w')
	f.write(lineout_aaPSSM)
	f.close()
		
		
def convert_protein_sequence_to_indexes( protein_seq ):
		sequence_indexes=np.zeros(len(protein_seq),int)
		position=-1
		"""aaA_1 aaR_1 aaN_1 aaD_1 aaC_1 aaE_1 aaQ_1 aaG_1 aaH_1 aaI_1 aaL_1 aaK_1 aaM_1 aaF_1 aaP_1 aaS_1 aaT_1 aaW_1 aaY_1 aaV_1"""
		for aa in protein_seq:
			position+=1
			if (aa == 'A'):
				sequence_indexes[position]=0
			elif(aa == "R"):
				sequence_indexes[position]=1
			elif(aa == "N"):
				sequence_indexes[position]=2
			elif(aa == "D"):
				sequence_indexes[position]=3
			elif(aa == "C"):
				sequence_indexes[position]=4
			elif(aa == "E"):
				sequence_indexes[position]=5
			elif(aa == "Q"):
				sequence_indexes[position]=6
			elif(aa == "G"):
				sequence_indexes[position]=7
			elif(aa == "H"):
				sequence_indexes[position]=8
			elif(aa == "I"):
				sequence_indexes[position]=9
			elif(aa == "L"):
				sequence_indexes[position]=10
			elif(aa == "K"):
				sequence_indexes[position]=11
			elif(aa == "M"):
				sequence_indexes[position]=12
			elif(aa == "F"):
				sequence_indexes[position]=13
			elif(aa == "P"):
				sequence_indexes[position]=14
			elif(aa == "S"):
				sequence_indexes[position]=15
			elif(aa == "T"):
				sequence_indexes[position]=16
			elif(aa == "W"):
				sequence_indexes[position]=17
			elif(aa == "Y"):
				sequence_indexes[position]=18
			elif(aa == "V"):
				sequence_indexes[position]=19
			else:
				print "ERROR!!! That aa does not exist?!"
				assert(0==1)
				
		return sequence_indexes
		
		
#Compare a sequence against a profile
def compare_seq_2polar_profile(seq, 
								profile_dic):
	for i in range(len(seq)):
		if (seq[i] not in profile_dic[i]):
			return False
	return True
	
		
class database():
	#Specific imports from interface_fragment_matching libraries
	import interface_fragment_matching.fragment_fitting.clustering as fragment_matching
	##from interface_fragment_matching.fragment_fitting import FragmentDatabase as FragmentDatabase
	##from interface_fragment_matching.fragment_fitting.clustering.cluster_manager import FragmentClusterManager
	###from interface_fragment_matching.fragment_fitting.clustering.cluster_manager import FragmentHierarchicalClusterManager
	from interface_fragment_matching.fragment_fitting.rmsd_calc import rmsd_calc as FragmentRmsdCalculator
	#Class private variables
	##1. Clustering databases:
	num_threads_global=None
	rmsd_calculator=None
	cluster_lr_centers_store_data=None
	cluster_lr_centers_coordinates=None
	cluster_lr_centers_store_specs=None
	assignments_lr_store_data=None
	cluster_hr_centers_store_data=None
	cluster_hr_centers_coordinates=None
	cluster_hr_centers_store_specs=None
	assignments_hr_store_data=None
	clusters_lr_population=None
	clusters_hr_population=None
	clusters_lr_TCM=None
	clusters_lr_TCM_ndx_p=None
	clusters_lr_TCM_ndx_len_p=None
	correspondence_dic_LR_HR = None
	correspondence_dic_HR_LR = None
	aa_seq_dic=None
	aa_seq_2ndx_dic=None
	aa_assignments_ndx_dic=None
	
	def __init__(self):
		#Instantiate Asford fast RMSD calculation Code  implementation based on my original fast RMSD Code based in QCP
		self.num_threads_global=10
		self.rmsd_calculator = self.FragmentRmsdCalculator(self.num_threads_global)
		print "Number of RMSD-OpenMP-threads: ", self.rmsd_calculator.num_threads
		
		
	def parse_dadriano_clustered_fragment_table(self, table_name, pert_pseudocount=0.000001):
		def get_fragment_entry_dtype(num_fragment_atoms):
			return np.dtype([("id", int), ("avg_distance", float), ("threshold_distance", float), 
							 ("coordinates", float, ((num_fragment_atoms), 3)), 
								("bb", float, (num_fragment_atoms/4, 3)), ("ss", float, (num_fragment_atoms/4, 3)),
								("aa", float, (num_fragment_atoms/4, 20)), ("aa_w", float, (num_fragment_atoms/4, 20)) ])
		fragment_table = pandas.read_csv(table_name, index_col=False,  delim_whitespace=True)
		num_coordinates = max([int(c.split("_")[1]) for c in fragment_table.columns if c.startswith("X")])
		fragment_entries = np.zeros(len(fragment_table), get_fragment_entry_dtype(num_coordinates))
		fragment_entries["id"] = fragment_table.Clust
		fragment_entries["avg_distance"] = fragment_table.AvgRad
		fragment_entries["threshold_distance"] = fragment_table.MaxRad
		for i in xrange(1, num_coordinates + 1):
			fragment_entries["coordinates"][:,i - 1,0] = fragment_table["X_%s" % i]
			fragment_entries["coordinates"][:,i - 1,1] = fragment_table["Y_%s" % i]
			fragment_entries["coordinates"][:,i - 1,2] = fragment_table["Z_%s" % i]
			
		"""aaA_1 aaR_1 aaN_1 aaD_1 aaC_1 aaE_1 aaQ_1 aaG_1 aaH_1 aaI_1 aaL_1 aaK_1 aaM_1 aaF_1 aaP_1 aaS_1 aaT_1 aaW_1 aaY_1 aaV_1"""
		for i in xrange(1, (num_coordinates/4)+1):
			fragment_entries["bb"][:,i - 1,0 ] = fragment_table["Phi_%s" % i]
			fragment_entries["bb"][:,i - 1,1 ] = fragment_table["Psi_%s" % i]
			fragment_entries["bb"][:,i - 1,2 ] = fragment_table["Ome_%s" % i]
			
			fragment_entries["ss"][:,i - 1,0 ] = fragment_table["ssH_%s" % i] + fragment_table["ssG_%s" % i] + fragment_table["ssI_%s" % i]
			fragment_entries["ss"][:,i - 1,1 ] = fragment_table["ssB_%s" % i] + fragment_table["ssE_%s" % i] 
			fragment_entries["ss"][:,i - 1,2 ] = fragment_table["ssT_%s" % i] + fragment_table["ssS_%s" % i] + fragment_table["ssL_%s" % i]
			
			fragment_entries["aa"][:,i - 1,0 ] = fragment_table["aaA_%s" % i]
			fragment_entries["aa"][:,i - 1,1 ] = fragment_table["aaR_%s" % i]
			fragment_entries["aa"][:,i - 1,2 ] = fragment_table["aaN_%s" % i]
			fragment_entries["aa"][:,i - 1,3 ] = fragment_table["aaD_%s" % i]
			fragment_entries["aa"][:,i - 1,4 ] = fragment_table["aaC_%s" % i]
			fragment_entries["aa"][:,i - 1,5 ] = fragment_table["aaE_%s" % i]
			fragment_entries["aa"][:,i - 1,6 ] = fragment_table["aaQ_%s" % i]
			fragment_entries["aa"][:,i - 1,7 ] = fragment_table["aaG_%s" % i]
			fragment_entries["aa"][:,i - 1,8 ] = fragment_table["aaH_%s" % i]
			fragment_entries["aa"][:,i - 1,9 ] = fragment_table["aaI_%s" % i]
			fragment_entries["aa"][:,i - 1,10 ] = fragment_table["aaL_%s" % i]
			fragment_entries["aa"][:,i - 1,11 ] = fragment_table["aaK_%s" % i]
			fragment_entries["aa"][:,i - 1,12 ] = fragment_table["aaM_%s" % i]
			fragment_entries["aa"][:,i - 1,13 ] = fragment_table["aaF_%s" % i]
			fragment_entries["aa"][:,i - 1,14 ] = fragment_table["aaP_%s" % i]
			fragment_entries["aa"][:,i - 1,15 ] = fragment_table["aaS_%s" % i]
			fragment_entries["aa"][:,i - 1,16 ] = fragment_table["aaT_%s" % i]
			fragment_entries["aa"][:,i - 1,17 ] = fragment_table["aaW_%s" % i]
			fragment_entries["aa"][:,i - 1,18 ] = fragment_table["aaY_%s" % i]
			fragment_entries["aa"][:,i - 1,19 ] = fragment_table["aaV_%s" % i]
			"""
			fragment_entries["aa_w"][:,i - 1,0 ] = np.log2((fragment_table["aaA_%s" % i]/0.07805)+pert_pseudocount)
			fragment_entries["aa_w"][:,i - 1,1 ] = np.log2((fragment_table["aaR_%s" % i]/0.05129)+pert_pseudocount)
			fragment_entries["aa_w"][:,i - 1,2 ] = np.log2((fragment_table["aaN_%s" % i]/0.04487)+pert_pseudocount)
			fragment_entries["aa_w"][:,i - 1,3 ] = np.log2((fragment_table["aaD_%s" % i]/0.05364)+pert_pseudocount)
			fragment_entries["aa_w"][:,i - 1,4 ] = np.log2((fragment_table["aaC_%s" % i]/0.01925)+pert_pseudocount)
			fragment_entries["aa_w"][:,i - 1,5 ] = np.log2((fragment_table["aaE_%s" % i]/0.06295)+pert_pseudocount)
			fragment_entries["aa_w"][:,i - 1,6 ] = np.log2((fragment_table["aaQ_%s" % i]/0.04264)+pert_pseudocount)
			fragment_entries["aa_w"][:,i - 1,7 ] = np.log2((fragment_table["aaG_%s" % i]/0.07377)+pert_pseudocount)
			fragment_entries["aa_w"][:,i - 1,8 ] = np.log2((fragment_table["aaH_%s" % i]/0.02199)+pert_pseudocount)
			fragment_entries["aa_w"][:,i - 1,9 ] = np.log2((fragment_table["aaI_%s" % i]/0.05142)+pert_pseudocount)
			fragment_entries["aa_w"][:,i - 1,10 ] = np.log2((fragment_table["aaL_%s" % i]/0.09019)+pert_pseudocount)
			fragment_entries["aa_w"][:,i - 1,11 ] = np.log2((fragment_table["aaK_%s" % i]/0.05744)+pert_pseudocount)
			fragment_entries["aa_w"][:,i - 1,12 ] = np.log2((fragment_table["aaM_%s" % i]/0.02243)+pert_pseudocount)
			fragment_entries["aa_w"][:,i - 1,13 ] = np.log2((fragment_table["aaF_%s" % i]/0.03856)+pert_pseudocount)
			fragment_entries["aa_w"][:,i - 1,14 ] = np.log2((fragment_table["aaP_%s" % i]/0.05203)+pert_pseudocount)
			fragment_entries["aa_w"][:,i - 1,15 ] = np.log2((fragment_table["aaS_%s" % i]/0.07120)+pert_pseudocount)
			fragment_entries["aa_w"][:,i - 1,16 ] = np.log2((fragment_table["aaT_%s" % i]/0.05841)+pert_pseudocount)
			fragment_entries["aa_w"][:,i - 1,17 ] = np.log2((fragment_table["aaW_%s" % i]/0.01330)+pert_pseudocount)
			fragment_entries["aa_w"][:,i - 1,18 ] = np.log2((fragment_table["aaY_%s" % i]/0.03216)+pert_pseudocount)
			fragment_entries["aa_w"][:,i - 1,19 ] = np.log2((fragment_table["aaV_%s" % i]/0.06441)+pert_pseudocount)
			"""
		fragment_spec = self.fragment_matching.cluster_manager.FragmentSpecification(4, ("CA", "C", "O", "N"))
		return fragment_spec, fragment_spec.fragment_data_to_location_and_atom(fragment_entries)
		
		
	def parse_dadriano_clustered_fragment_assignments(self, table_name):
		def get_assignment_entry_dtype():
			return np.dtype([("name", str, 200), ("span", int, 2) , ("assignment", int), ("distance", float), ("aaSequence", str, 4)])
		assignment_table = pandas.read_csv(table_name, index_col=False,  delim_whitespace=True)
		assignment_entries = np.zeros(len(assignment_table), get_assignment_entry_dtype())
		assignment_entries["name"] = assignment_table.Structure
		assignment_entries["span"][:,0] = assignment_table.SpanI
		assignment_entries["span"][:,1] = assignment_table.SpanE
		assignment_entries["assignment"] = assignment_table.Assignment
		assignment_entries["distance"] = assignment_table.Distance
		assignment_entries["aaSequence"] = assignment_table.aaSequence
		return assignment_entries
		
		
		
	def read_clustering_databases(self, clustering_database_path="", lres="0.5", hres="0.2", tcm_jump_step=3):
		print "Reading clustering databases and calculating on-the-fly information, be patient"
		#Read clustering (original: "dadriano version")
		print "Reading vall cluster databases (LR and HR)"
		#LR clustering 
		lr_stats_in_path="%s/kcenters_stats.dat_r%s"%(clustering_database_path, lres)
		print "Reading LR clustering from: ", lr_stats_in_path
		self.cluster_lr_centers_store_specs, self.cluster_lr_centers_store_data = self.parse_dadriano_clustered_fragment_table(lr_stats_in_path)
		self.cluster_lr_centers_coordinates = self.cluster_lr_centers_store_data["coordinates"]
		self.cluster_lr_centers_coordinates = self.cluster_lr_centers_coordinates.copy().view(
									dtype=float).reshape((self.cluster_lr_centers_coordinates.shape[0], -1, 3))
		#HR clustering
		hr_stats_in_path="%s/kcenters_stats.dat_r%s"%(clustering_database_path, hres)
		print "Reading LR clustering from: ", hr_stats_in_path
		self.cluster_hr_centers_store_specs, self.cluster_hr_centers_store_data = self.parse_dadriano_clustered_fragment_table(hr_stats_in_path)
		self.cluster_hr_centers_coordinates = self.cluster_hr_centers_store_data["coordinates"]
		self.cluster_hr_centers_coordinates = self.cluster_hr_centers_coordinates.copy().view(
									dtype=float).reshape((self.cluster_hr_centers_coordinates.shape[0], -1, 3))
		#Read assignments (dadriano version)
		print "Reading assignments (LR and HR)"
		#LR assignments
		lr_assign_in_path="%s/kcenters_assignments.dat_r%s"%(clustering_database_path, lres)
		print "Reading LR vall assignments from: ", lr_assign_in_path
		self.assignments_lr_store_data = self.parse_dadriano_clustered_fragment_assignments(lr_assign_in_path)
		#HR assignments
		hr_assign_in_path="%s/kcenters_assignments.dat_r%s"%(clustering_database_path, hres)
		print "Reading LR vall assignments from: ", hr_assign_in_path
		self.assignments_hr_store_data = self.parse_dadriano_clustered_fragment_assignments(hr_assign_in_path)
		print "Done Reading Databases"
		
		#Calculate clustering stats
		print "Calculating clustering stats on-the-fly"
		#LR
		num_clusters_lr=len(self.cluster_lr_centers_coordinates)
		#LR calculate the population of each cluster
		self.clusters_lr_population = np.zeros(num_clusters_lr, int) 
		for i in self.assignments_lr_store_data["assignment"]:
			#Assignmetns start in 1, so -1 it!
			self.clusters_lr_population[i-1]+=1
		#LR calculate the TCM
		self.clusters_lr_TCM = np.zeros((num_clusters_lr,num_clusters_lr), int)
		print "Fragments in LR dataset: ", len(self.assignments_lr_store_data["assignment"])
		print "Clusters from LR dataset: ", num_clusters_lr
		#Get contiguous indexes
		tmpIndex=0
		assignments_lr_contiguousIndex= np.zeros(len(self.assignments_lr_store_data["assignment"]),int)
		assignments_lr_contiguousIndex[0]=tmpIndex
		for i in range(0,len(self.assignments_lr_store_data["assignment"])-1):
			if ( (self.assignments_lr_store_data["span"][i,0]) == (self.assignments_lr_store_data["span"][i+1,0]-1) ):
				assignments_lr_contiguousIndex[i+1]=tmpIndex
			else:
				tmpIndex+=1
				assignments_lr_contiguousIndex[i+1]=tmpIndex
		print "Number of contiguous LR fragments: ", assignments_lr_contiguousIndex.max()+1
		#Calculate TCM with jump step=arg.jump step
		for i in range(0,len(self.assignments_lr_store_data["assignment"])-tcm_jump_step):
			if ( assignments_lr_contiguousIndex[i] == assignments_lr_contiguousIndex[i+tcm_jump_step]):
				self.clusters_lr_TCM[(self.assignments_lr_store_data["assignment"][i]-1)][(self.assignments_lr_store_data["assignment"][i+tcm_jump_step]-1)]+=1
		print "LR TCM created: ", self.clusters_lr_TCM.shape, ", with jump step: ", tcm_jump_step
		#HR
		num_clusters_hr=len(self.cluster_hr_centers_coordinates)
		#HR calculate the population of each cluster
		self.clusters_hr_population = np.zeros(num_clusters_hr, int) 
		for i in self.assignments_hr_store_data["assignment"]:
			self.clusters_hr_population[i-1]+=1
		print "Fragments in HR dataset: ", len(self.assignments_hr_store_data["assignment"])
		print "Clusters from HR dataset: ", num_clusters_hr
		print "Done calculating stats"
		
		#Pre calculate tcm bitmasks and counts at trimmed cut-offs
		#Calculate trimmed cut-offs indexes in the TCM
		print "Precalculating clustered databases trimmed bitmasks"
		tmp_clusters_lr_TCM_ndx_p={}
		tmp_clusters_lr_TCM_ndx_len_p={}
		tmp_clusters_lr_TCM_ndx_p[1]=[]
		tmp_clusters_lr_TCM_ndx_p[10]=[]
		tmp_clusters_lr_TCM_ndx_p[20]=[]
		tmp_clusters_lr_TCM_ndx_p[30]=[]
		tmp_clusters_lr_TCM_ndx_p[40]=[]
		tmp_clusters_lr_TCM_ndx_p[50]=[]
		tmp_clusters_lr_TCM_ndx_p[60]=[]
		tmp_clusters_lr_TCM_ndx_p[70]=[]
		tmp_clusters_lr_TCM_ndx_p[80]=[]
		tmp_clusters_lr_TCM_ndx_p[90]=[]
		tmp_clusters_lr_TCM_ndx_p[100]=[]
		tmp_clusters_lr_TCM_ndx_p[200]=[]
		tmp_clusters_lr_TCM_ndx_len_p[1]=[]
		tmp_clusters_lr_TCM_ndx_len_p[10]=[]
		tmp_clusters_lr_TCM_ndx_len_p[20]=[]
		tmp_clusters_lr_TCM_ndx_len_p[30]=[]
		tmp_clusters_lr_TCM_ndx_len_p[40]=[]
		tmp_clusters_lr_TCM_ndx_len_p[50]=[]
		tmp_clusters_lr_TCM_ndx_len_p[60]=[]
		tmp_clusters_lr_TCM_ndx_len_p[70]=[]
		tmp_clusters_lr_TCM_ndx_len_p[80]=[]
		tmp_clusters_lr_TCM_ndx_len_p[90]=[]
		tmp_clusters_lr_TCM_ndx_len_p[100]=[]
		tmp_clusters_lr_TCM_ndx_len_p[200]=[]
		for i in range(len(self.clusters_lr_TCM)):
			tmp_clusters_lr_TCM_ndx_p[1].append(np.where(self.clusters_lr_TCM[i] >= 1)[0])
			tmp_clusters_lr_TCM_ndx_p[10].append(np.where(self.clusters_lr_TCM[i] >= 10)[0])
			tmp_clusters_lr_TCM_ndx_p[20].append(np.where(self.clusters_lr_TCM[i] >= 20)[0])
			tmp_clusters_lr_TCM_ndx_p[30].append(np.where(self.clusters_lr_TCM[i] >= 30)[0])
			tmp_clusters_lr_TCM_ndx_p[40].append(np.where(self.clusters_lr_TCM[i] >= 40)[0])
			tmp_clusters_lr_TCM_ndx_p[50].append(np.where(self.clusters_lr_TCM[i] >= 50)[0])
			tmp_clusters_lr_TCM_ndx_p[60].append(np.where(self.clusters_lr_TCM[i] >= 60)[0])
			tmp_clusters_lr_TCM_ndx_p[70].append(np.where(self.clusters_lr_TCM[i] >= 70)[0])
			tmp_clusters_lr_TCM_ndx_p[80].append(np.where(self.clusters_lr_TCM[i] >= 80)[0])
			tmp_clusters_lr_TCM_ndx_p[90].append(np.where(self.clusters_lr_TCM[i] >= 90)[0])
			tmp_clusters_lr_TCM_ndx_p[100].append(np.where(self.clusters_lr_TCM[i] >= 100)[0])
			tmp_clusters_lr_TCM_ndx_p[200].append(np.where(self.clusters_lr_TCM[i] >= 200)[0])
			#append the len of the last added element
			tmp_clusters_lr_TCM_ndx_len_p[1].append(len(tmp_clusters_lr_TCM_ndx_p[1][-1]))
			tmp_clusters_lr_TCM_ndx_len_p[10].append(len(tmp_clusters_lr_TCM_ndx_p[10][-1]))
			tmp_clusters_lr_TCM_ndx_len_p[20].append(len(tmp_clusters_lr_TCM_ndx_p[20][-1]))
			tmp_clusters_lr_TCM_ndx_len_p[30].append(len(tmp_clusters_lr_TCM_ndx_p[30][-1]))
			tmp_clusters_lr_TCM_ndx_len_p[40].append(len(tmp_clusters_lr_TCM_ndx_p[40][-1]))
			tmp_clusters_lr_TCM_ndx_len_p[50].append(len(tmp_clusters_lr_TCM_ndx_p[50][-1]))
			tmp_clusters_lr_TCM_ndx_len_p[60].append(len(tmp_clusters_lr_TCM_ndx_p[60][-1]))
			tmp_clusters_lr_TCM_ndx_len_p[70].append(len(tmp_clusters_lr_TCM_ndx_p[70][-1]))
			tmp_clusters_lr_TCM_ndx_len_p[80].append(len(tmp_clusters_lr_TCM_ndx_p[80][-1]))
			tmp_clusters_lr_TCM_ndx_len_p[90].append(len(tmp_clusters_lr_TCM_ndx_p[90][-1]))
			tmp_clusters_lr_TCM_ndx_len_p[100].append(len(tmp_clusters_lr_TCM_ndx_p[100][-1]))
			tmp_clusters_lr_TCM_ndx_len_p[200].append(len(tmp_clusters_lr_TCM_ndx_p[200][-1]))
		#Convert to np arrays
		self.clusters_lr_TCM_ndx_p={}
		self.clusters_lr_TCM_ndx_len_p={}
		self.clusters_lr_TCM_ndx_p[1]      =np.asarray(tmp_clusters_lr_TCM_ndx_p[1])
		self.clusters_lr_TCM_ndx_p[10]     =np.asarray(tmp_clusters_lr_TCM_ndx_p[10])
		self.clusters_lr_TCM_ndx_p[20]     =np.asarray(tmp_clusters_lr_TCM_ndx_p[20])
		self.clusters_lr_TCM_ndx_p[30]     =np.asarray(tmp_clusters_lr_TCM_ndx_p[30])
		self.clusters_lr_TCM_ndx_p[40]     =np.asarray(tmp_clusters_lr_TCM_ndx_p[40])
		self.clusters_lr_TCM_ndx_p[50]     =np.asarray(tmp_clusters_lr_TCM_ndx_p[50])
		self.clusters_lr_TCM_ndx_p[60]     =np.asarray(tmp_clusters_lr_TCM_ndx_p[60])
		self.clusters_lr_TCM_ndx_p[70]     =np.asarray(tmp_clusters_lr_TCM_ndx_p[70])
		self.clusters_lr_TCM_ndx_p[80]     =np.asarray(tmp_clusters_lr_TCM_ndx_p[80])
		self.clusters_lr_TCM_ndx_p[90]     =np.asarray(tmp_clusters_lr_TCM_ndx_p[90])
		self.clusters_lr_TCM_ndx_p[100]    =np.asarray(tmp_clusters_lr_TCM_ndx_p[100])
		self.clusters_lr_TCM_ndx_p[200]    =np.asarray(tmp_clusters_lr_TCM_ndx_p[200])
		self.clusters_lr_TCM_ndx_len_p[1]  =np.asarray(tmp_clusters_lr_TCM_ndx_len_p[1])
		self.clusters_lr_TCM_ndx_len_p[10] =np.asarray(tmp_clusters_lr_TCM_ndx_len_p[10])
		self.clusters_lr_TCM_ndx_len_p[20] =np.asarray(tmp_clusters_lr_TCM_ndx_len_p[20])
		self.clusters_lr_TCM_ndx_len_p[30] =np.asarray(tmp_clusters_lr_TCM_ndx_len_p[30])
		self.clusters_lr_TCM_ndx_len_p[40] =np.asarray(tmp_clusters_lr_TCM_ndx_len_p[40])
		self.clusters_lr_TCM_ndx_len_p[50] =np.asarray(tmp_clusters_lr_TCM_ndx_len_p[50])
		self.clusters_lr_TCM_ndx_len_p[60] =np.asarray(tmp_clusters_lr_TCM_ndx_len_p[60])
		self.clusters_lr_TCM_ndx_len_p[70] =np.asarray(tmp_clusters_lr_TCM_ndx_len_p[70])
		self.clusters_lr_TCM_ndx_len_p[80] =np.asarray(tmp_clusters_lr_TCM_ndx_len_p[80])
		self.clusters_lr_TCM_ndx_len_p[90] =np.asarray(tmp_clusters_lr_TCM_ndx_len_p[90])
		self.clusters_lr_TCM_ndx_len_p[100]=np.asarray(tmp_clusters_lr_TCM_ndx_len_p[100])
		self.clusters_lr_TCM_ndx_len_p[200]=np.asarray(tmp_clusters_lr_TCM_ndx_len_p[200])
		print "Done calculating databases trimmed bitmasks"
		
		#Generate the correspondence table between the HR and LR clusterings and viceversa
		print "Creating LR to HR (and viceversa) correspondence dictionaries"
		correspondence_table_LR_HR = []
		correspondence_table_HR_LR = []
		for i in range(len(self.assignments_lr_store_data["assignment"])):
				correspondence_table_LR_HR.append( 
							(self.assignments_lr_store_data["assignment"][i]-1,
							 self.assignments_hr_store_data["assignment"][i]-1) )
				correspondence_table_HR_LR.append( 
							(self.assignments_hr_store_data["assignment"][i]-1,
							 self.assignments_lr_store_data["assignment"][i]-1) )
		#Remove duplicates
		correspondence_table_LR_HR = list(OrderedDict.fromkeys(correspondence_table_LR_HR))
		correspondence_table_HR_LR = list(OrderedDict.fromkeys(correspondence_table_HR_LR))
		
		#Convert the correspondence table to dictionary :
		print "Creating LR to HR (and viceversa) dictionary"
		self.correspondence_dic_LR_HR={}
		for i in correspondence_table_LR_HR:
			c=i[0]
			itm=i[1]
			if not c in self.correspondence_dic_LR_HR:
				self.correspondence_dic_LR_HR[c] = [itm]
			else:
				self.correspondence_dic_LR_HR[c].append(itm)
		#Convert to dictionary :  
		self.correspondence_dic_HR_LR={}
		for i in correspondence_table_HR_LR:
			c=i[0]
			itm=i[1]
			if not c in self.correspondence_dic_HR_LR:
				self.correspondence_dic_HR_LR[c] = [itm]
			else:
				self.correspondence_dic_HR_LR[c].append(itm)
		print "Done with correspondence dictionaries"
		
		#Convert assignments to a dictionary of the observed sequences:  
		print "Creating seq count/ndx/assignment dictionaries"
		self.aa_seq_dic={}
		for i in self.assignments_lr_store_data["aaSequence"]:
			if not i in self.aa_seq_dic:
				self.aa_seq_dic[i] = [1]
			else:
				self.aa_seq_dic[i][0]+=1
				
		self.aa_seq_2ndx_dic={}
		counter=0
		for i in self.assignments_lr_store_data["aaSequence"]:
			if not i in self.aa_seq_2ndx_dic:
				self.aa_seq_2ndx_dic[i] = [counter]
			else:
				self.aa_seq_2ndx_dic[i].append(counter)
			counter+=1
			
		#Convert assignments to a dictionary of the observed sequences:  
		self.aa_assignments_ndx_dic={}
		counter=0
		#Remember to -1 because the assignments start on 1
		for i in (self.assignments_lr_store_data["assignment"]-1):
			if not i in self.aa_assignments_ndx_dic:
				self.aa_assignments_ndx_dic[i] = [counter]
			else:
				self.aa_assignments_ndx_dic[i].append(counter)
			counter+=1
		print "Done with seq dictionaries"
		
		print "Finally done reading Clustering databases!!!"
		
class proteinLayers():
	##2a. Design related layers:
	CoreAAdic=None
	HYDROPHIBICaa=None
	NEUTRALaa=None
	CoreAAdic_extended=None
	CoreAAdic=None
	BoundaryAAdic=None
	SurfAAdic=None
	SepecialAAdic=None
	
	def __init__(self):
		print "Creating Defaul profile a.a. dictionaries"
		self.create_aa_profile_dictionaries()
		self.check_aa_layer_dictionaries()
		print "Done"
		
	def create_aa_profile_dictionaries(self, POLARaa={'R','N','D','C','E','Q','H','K','S','T','Y'},
										HYDROPHIBICaa={'I','L','M','F','P','W','V'},
										NEUTRALaa={'A','G'},
										CoreAAdic_extended={'A','F','I','L','M','P','V','W','Y','D','N','S','T'},
										CoreAAdic={'A','F','I','L','M','P','V','W'},
										BoundaryAAdic={'A','D','E','F','G','I','K','L','M','N','P','Q','R','S','T','V','W','Y'},
										SurfAAdic={'D','E','G','H','K','N','P','Q','R','S','T'},
										SepecialAAdic={'C','U'}):
		#By polarity
		self.CoreAAdic			=CoreAAdic
		self.HYDROPHIBICaa		=HYDROPHIBICaa
		self.NEUTRALaa			=NEUTRALaa
		#By localization
		self.CoreAAdic_extended	=CoreAAdic_extended
		self.CoreAAdic			=CoreAAdic
		self.BoundaryAAdic		=BoundaryAAdic
		self.SurfAAdic			=SurfAAdic
		self.SepecialAAdic		=SepecialAAdic
		
	def check_aa_layer_dictionaries(self):
		CheckDic=[]
		for i in self.CoreAAdic_extended: CheckDic.append(i)
		for i in self.CoreAAdic: CheckDic.append(i)
		for i in self.SurfAAdic: CheckDic.append(i)
		for i in self.SepecialAAdic: CheckDic.append(i)
		for i in self.BoundaryAAdic: CheckDic.append(i)
		CheckDic=np.unique(np.asarray(CheckDic))
		print "Num a.a. in layer dictionary:", len(CheckDic),"; aas:", CheckDic
		if (len(CheckDic) == 21):
			print "All seems correct, all 21 a.a. are considered in the layer dictionary"
		else:
			print "***WARNING! Only", len(CheckDic), "a.a. are present in the layer dictionary instead of 21."
			print "***If this is not intended something must be wrong in your setup. You have been warned!"
		
class idealizer():
	#Some rosetta imports
	import rosetta as rosetta
	from rosetta.core.optimization import CartesianMinimizer as CartMin
	from rosetta.core.optimization import AtomTreeMinimizer as mzr
	from rosetta.core.scoring.constraints import HarmonicFunc as HarmFunc
	from rosetta.core.scoring.constraints import CoordinateConstraint as CoorCstr
	from rosetta.protocols.relax import FastRelax as FastRel
	from rosetta.protocols.simple_moves import PackRotamersMover as PackRot
	from rosetta.protocols.idealize import IdealizeMover as IdealBondMover
	##from rosetta.core import sequence as rosetta_sequence
	##from rosetta.protocols.relax import ClassicRelax as ClassRel
	##from rosetta.protocols.simple_moves import FavorSequenceProfile as FSP 
	
	##Global variables
	fake_aminoacid_letter_for_design=None
	#sfxs
	scorefxn_tal=None
	scorefxn_tal_low_farep=None
	scorefxn_tal_low_farep_sstrand=None
	scorefxn_soft_rep=None
	scorefxn_remodel_cen=None
	
	def __init__(self):
		print "Setting Fake aminoacid letter for design to Default"
		self.fake_aminoacid_letter_for_design='A'
		print "Default aminoacid for scaffold design: ", self.fake_aminoacid_letter_for_design
		self.init_rosetta("")
		self.create_energy_functions()
		
		
	def create_energy_functions(self):
		print "Creating energy functions"
		scorefxn_tal_name="talaris2013"
		self.scorefxn_tal = self.rosetta.core.scoring.ScoreFunctionFactory.create_score_function(scorefxn_tal_name)
		
		scorefxn_tal_low_farep_name="talaris2013"
		self.scorefxn_tal_low_farep = self.rosetta.core.scoring.ScoreFunctionFactory.create_score_function(scorefxn_tal_low_farep_name)
		self.scorefxn_tal_low_farep.set_weight(self.rosetta.core.scoring.fa_rep, 0.20)
		
		scorefxn_tal_low_farep_name_sstrand="talaris2013"
		self.scorefxn_tal_low_farep_sstrand = self.rosetta.core.scoring.ScoreFunctionFactory.create_score_function(scorefxn_tal_low_farep_name_sstrand)
		self.scorefxn_tal_low_farep_sstrand.set_weight(self.rosetta.core.scoring.ss_pair, 1.0)
		self.scorefxn_tal_low_farep_sstrand.set_weight(self.rosetta.core.scoring.hbond_bb_sc, 1.0)
		self.scorefxn_tal_low_farep_sstrand.set_weight(self.rosetta.core.scoring.hbond_sr_bb, 1.0)
		self.scorefxn_tal_low_farep_sstrand.set_weight(self.rosetta.core.scoring.hbond_lr_bb, 1.0)
		
		scorefxn_soft_rep_name="soft_rep"
		self.scorefxn_soft_rep = self.rosetta.core.scoring.ScoreFunctionFactory.create_score_function(scorefxn_soft_rep_name)
		
		scorefxn_remodel_cen_name="remodel_cen"
		self.scorefxn_remodel_cen = self.rosetta.core.scoring.ScoreFunctionFactory.create_score_function(scorefxn_remodel_cen_name)
		print "Done creating energy functions"
	
	
	def init_rosetta(self, args=" -mute all "):
		##Init Rosetta and Rosetta database(s)
		print "Initiating Rosetta, be patient"
		self.rosetta.init(args)
		print "Rosseta loading Done"
		
	def read_input_pdbs_stage1(self, 
						general_input_pdb="", 
						general_context_pdb_name="",
						general_out_dir="./",
						stage1_tmp_dir="stage1_tmp"):
		print "Setting global input_output variables"
		print "Working dir root: ", general_out_dir
		#Input PDB
		general_input_pdb=("%s/%s"%(general_out_dir,general_input_pdb)).strip()
		print "Input PDB: ", general_input_pdb
		#Input context pose
		general_context_pdb_name=("%s/stage1_context.pdb"%general_out_dir).strip()
		print "Context PDB: ", general_context_pdb_name
		b_make_plots = True
		print "Plot generation set to: ", b_make_plots
		out_pssm_filename = ("%s/pssm_design.pssm"%general_out_dir).strip()
		#Read input pose
		print "Reading input"
		#Input PDB
		inputPose_original=self.rosetta.pose_from_pdb(general_input_pdb)
		print "Done, thanks :) I have readed", general_input_pdb, ", a pose with: ", inputPose_original.total_residue(), "residues"
		
		context_pose=None
		context_pose_centroid=None
		b_has_context_pose=False
		if (os.path.isfile(general_context_pdb_name)):
			b_has_context_pose=True
			context_pose=self.rosetta.pose_from_pdb(general_context_pdb_name)
			#Convert it to centroid
			context_pose_centroid=context_pose.clone()
			to_centroid_mover = self.rosetta.protocols.simple_moves.SwitchResidueTypeSetMover ( "centroid" )
			to_centroid_mover.apply(context_pose_centroid)
			print "Ohh, I also readed a context pose: ", general_context_pdb_name, "a pose with: \
					", context_pose_centroid.total_residue(), "residues"
		out_pssm_filename_rosetta = self.rosetta.utility.file.FileName(out_pssm_filename)
		print "DONE reading inputs"
		
		#Now Idealize bond distances which is needed for the next step
		print "Idealizing input pose for idealization bond distances"
		inputPose_idealized_bonds=inputPose_original.clone()
		self.idealize_bond_distances(inputPose_idealized_bonds)
		
		tmp_out_name="%s/test_startP.pdb" % stage1_tmp_dir
		print "Dumped initial structure to: ", tmp_out_name
		inputPose_original.dump_pdb( tmp_out_name )
		
		rmsd_ori_new=self.align_by_bb(inputPose_idealized_bonds, inputPose_original)
		if (rmsd_ori_new > 0.3):
			print "Warning, after idealizing bond distances the RMSD is: ", rmsd_ori_new
		print "Done idealizing bond distances, RMSD to input: ", rmsd_ori_new
		
		return inputPose_original, inputPose_idealized_bonds, context_pose, context_pose_centroid, b_has_context_pose
		
		
	def split_pose_byto_SScore(self, 
								inputPose_original, 
								clustersDB=None, 
								general_out_dir="./",
								stage1_tmp_dir="stage1_tmp"):
		#Determine if the pose is multi chain
		b_read_by_chain=False
		poses_by_chain=inputPose_original.split_by_chain()
		
		if (len(poses_by_chain) > 1):
			print "Using separated chains as SS defnition, if you don't want this use a single chain as input"
			b_read_by_chain=True
			#Silly print
			##print len(poses_by_chain)
			for pose_i in poses_by_chain:
				#print pose_i.total_residue()
				#print pose_i.sequence() 
				if (pose_i.total_residue() < 4):
					print "Can't use less than 4 aminoacids in a given SS, sorry. Aborting!"
					assert (0==1)
		else:
			print "Using automatic SS detection mode"
			
		#Create angle representation of our pose (NOTE: the idealize protocol should be used before on the input pose)
		pose_array=[]
		if (not b_read_by_chain):
			print "Copying Idealized bond lenghts (and angles)"
			tmpline=""
			for i in range(0, inputPose_original.total_residue()):
				tmpline+=self.fake_aminoacid_letter_for_design
				
			inputPose=self.rosetta.pose_from_sequence(tmpline)
			
			#Apply the original angles
			for i in range (0, inputPose_original.total_residue()):
				inputPose.set_phi(i+1, inputPose_original.phi(i+1))
				inputPose.set_psi(i+1, inputPose_original.psi(i+1))
				inputPose.set_omega(i+1, inputPose_original.omega(i+1)) 
			#Align
			self.align_by_bb(inputPose, inputPose_original)
			print "Done"
			
			#User tunable regions
			##if (not b_read_by_chain):
			
			avoid_protein_regions=[]
			force_keep_protein_regions=[]
			
			#Read files with user's regions constraints
			avoidRegions_fileName= ("%s/avoid_regions.dat" % general_out_dir )
			if (os.path.isfile(avoidRegions_fileName)):
				avoid_protein_regions = np.loadtxt(avoidRegions_fileName, dtype=int)
				print "Forced to avoid structural regions: ", avoid_protein_regions
			
			keepRegions_fileName= ("%s/keep_regions.dat" % general_out_dir )
			if (os.path.isfile(keepRegions_fileName)):
				force_keep_protein_regions = np.loadtxt(keepRegions_fileName, dtype=int)
				print "Forced to keep structural regions: ", force_keep_protein_regions
				
			#Determine SS probability
			##if (not b_read_by_chain):
			pose_ss_probability=np.zeros((inputPose.total_residue(),3),float) 
			print "Determining SS probability using fragments "
			##rmsd_calculator = rmsd_calc(num_threads_global)
			for testres in range(1, inputPose.total_residue()-2): #inputPose.total_residue()-2):
					#Compute the span for a 4mer
					testresndx=[testres,testres+3]
					#Calculate the RMS
					rmsd_result=clustersDB.rmsd_calculator.get_broadcast_coordinate_rmsd(
								np.expand_dims( pose_res_to_array(inputPose, testresndx, atoms=["CA","C","O","N"]), 0), 
								clustersDB.cluster_lr_centers_coordinates.copy() 
								).ravel()
					#The SS is determined by the probability of the SS by the overlaping fragments
					counter=0
					for ss_prob in clustersDB.cluster_lr_centers_store_data["ss"][rmsd_result.argmin()]:
						pose_ss_probability[testres-1+counter]+=ss_prob
						counter+=1
			print "Done determining SS probability"
			
			##ToDO: SS probability of B-sheets is complicated, they seem very much like loops. 
			## -- therefore remove this hack, by training? 
			##if (not b_read_by_chain):
			#Merge DSSP prob with cluster prob
			pose_ss_probability_w = np.copy(pose_ss_probability)
			tmp_dssp = self.rosetta.core.scoring.dssp.Dssp(inputPose_original)
			tmp_dssp_ss =tmp_dssp.get_dssp_secstruct()
			
			ss_E_assignment_prob_low_limit=0.8
			for i in range(len(tmp_dssp_ss)):
				if ( tmp_dssp_ss[i] == 'H' ):
					if(2 == pose_ss_probability_w[i].argmax()):
						print "We are likely wrong for assignment 'L' of position: ", i+1, "using DSSP 'H'"
						pose_ss_probability_w[i]=[1,0,0]
					elif(1 == pose_ss_probability_w[i].argmax()):
						print "We are likely wrong for assignment 'E' of position: ", i+1, "using DSSP 'H'"
						pose_ss_probability_w[i]=[1,0,0]
					elif(0 == pose_ss_probability_w[i].argmax()):
						pose_ss_probability_w[i]=[1,0,0]
				elif ( tmp_dssp_ss[i] == 'E' ):
					if(1 != pose_ss_probability_w[i].argmax()):
						##if(pose_ss_probability_w[i][1] >= ss_E_assignment_prob_low_limit):
							print "We are likely wrong for assignment of position: ", i+1, "using DSSPs 'E'"
							pose_ss_probability_w[i]=[0,1,0]
						##else:
							#Carefull this could be wrong, anyway is a hack, fix
						##    print "DSSP is likely wrong for assignment of position: ", i+1, "using ours 'L'"
						##    pose_ss_probability_w[i]=[0,0,1]
					else:
						pose_ss_probability_w[i]=[0,1,0]
				elif ( tmp_dssp_ss[i] == 'L' ):
					if(2 != pose_ss_probability_w[i].argmax()):
						if(1 == pose_ss_probability_w[i].argmax()):
							print "DSSP is likely wrong for assignment 'L' of position: ", i+1, "using ours 'E'"
							pose_ss_probability_w[i]=[0,1,0]
						elif(0 == pose_ss_probability_w[i].argmax()):
							print "DSSP is likely wrong for assignment 'H' of position: ", i+1, "using ours 'H'"
							pose_ss_probability_w[i]=[1,0,0]
			
			#for i in range(len(pose_ss_probability_w)):
			#    print i+1, pose_ss_probability_w[i]
			
			##if (not b_read_by_chain):
			print "Current assignments (0==H 1==E 2==L)"
			tmp_ss_assignments_array = np.zeros(len(pose_ss_probability_w), int)
			for i in range(len(pose_ss_probability_w)):
				tmp_ss_assignments_array[i] = pose_ss_probability_w[i].argmax()
				#print i+1, tmp_ss_assignments_array[i]
			print tmp_ss_assignments_array
				
			##if (not b_read_by_chain):
			print "Enforcing user defined regions (0==H 1==E 2==L)"
			#Hack regions to avoid (i.e. convert them to loops)
			for region in avoid_protein_regions:
				print "Avoid", region
				for j in range(region[0], region[1]+1):
					tmp_ss_assignments_array[j-1]=2
					pose_ss_probability_w[j-1]=[0,0,1]
			for region in force_keep_protein_regions:
				print "Force", region
				for j in range(region[0], region[1]+1):
					tmp_ss_assignments_array[j-1]=3
					pose_ss_probability_w[j-1]=[1,1,0]
					
			print tmp_ss_assignments_array, len(tmp_ss_assignments_array)
			#for i in range(len(tmp_ss_assignments_array)):
			#    print i, tmp_ss_assignments_array[i]
			
			##if (not b_read_by_chain):
			#Calculate SS assignments (0=H, 1=E, 2=L) and contiguous SS fragments
			print "Determining contiguous SS fragments"
			#g_counter=0
			#l_counter=1
			#current_topology=-1
			
			low_cont_aa_lim_for_ss=3
			
			tmp_ss_assignments_array = np.zeros(len(pose_ss_probability_w), int)
			lastSS_res=-1
			contigous_fragments=[]
			for i in range(len(pose_ss_probability_w)):
				tmp_ss_assignments_array[i] = pose_ss_probability_w[i].argmax()
				if (tmp_ss_assignments_array[i] == 2):
					if (i-(low_cont_aa_lim_for_ss-1) > lastSS_res+1):   
						#The extra +1 is to be compatible with rosetta numbers
						contigous_fragments.append([lastSS_res+1+1, i])
					lastSS_res=i
			if ((len(pose_ss_probability_w)-1- lastSS_res) > (low_cont_aa_lim_for_ss-1)):
				contigous_fragments.append([lastSS_res+2, len(pose_ss_probability_w)-1])
				
			print "Found %d contiguous SS Fragments: " %(len(contigous_fragments)), contigous_fragments
			
			#Merge those separated by only one residue
			"""
			if False:
				print "Merging those SS separated by only one residue..."
				needs_merge=True
				while(needs_merge):
					needs_merge=False
					contigous_fragments_merged=[]
					i = 0
					while(i<len(contigous_fragments)-1):
						if ((contigous_fragments[i][1]+2) == contigous_fragments[i+1][0]):
							#print  [contigous_fragments[i][0],contigous_fragments[i+1][1]]
							contigous_fragments_merged.append([contigous_fragments[i][0],contigous_fragments[i+1][1]])
							needs_merge=True
							i+=1
						else:
							#print contigous_fragments[i], contigous_fragments[i+1]
							contigous_fragments_merged.append(contigous_fragments[i][:])
						i+=1
					if(not needs_merge):
						contigous_fragments_merged.append(contigous_fragments[-1][:])
					contigous_fragments=contigous_fragments_merged[:]
				#print contigous_fragments
			print "\n Now (after merging) I have found %d contiguous SS fragment spans : "%(len(contigous_fragments)), contigous_fragments
			"""
			
			#Copy contiguous SS elements (H and E) to separate poses
			##if (not b_read_by_chain):
			print "Copying  %d contigous SS elements to separate poses" % ( len(contigous_fragments) )
			pose_array=[]
			for segmentndx in contigous_fragments:
				tmpline=""
				for i in range(segmentndx[0], segmentndx[1]+1):
					tmpline+=self.fake_aminoacid_letter_for_design
				tmpPose=self.rosetta.pose_from_sequence(tmpline)
				tmpPose.copy_segment((segmentndx[1]-segmentndx[0]+1),inputPose,1,segmentndx[0])
				pose_array.append(tmpPose.clone())
				##pose_array[len(pose_array)-1].dump_pdb("%s/testSS_before_ideal_%02d.pdb"%(general_out_dir,len(pose_array)-1))
		#If in separate chains:
		elif (b_read_by_chain):
			print "Copying  %d contigous SS elements to separate poses" % ( len(poses_by_chain) )
			pose_array=[]
			
			for i in range(1, (len(poses_by_chain)+1)):
				print poses_by_chain[i].total_residue()
				print poses_by_chain[i].sequence() 
				print "Idealizing bond lenghts (and angles) of chain", i
				tmpline=""
				for j in range(0, poses_by_chain[i].total_residue()):
					tmpline+=self.fake_aminoacid_letter_for_design    
				inputPose_tmp=self.rosetta.pose_from_sequence(tmpline)
				#Apply the original angles
				for j in range (0, poses_by_chain[i].total_residue()):
					inputPose_tmp.set_phi(j+1, poses_by_chain[i].phi(j+1))
					inputPose_tmp.set_psi(j+1, poses_by_chain[i].psi(j+1))
					inputPose_tmp.set_omega(j+1, poses_by_chain[i].omega(j+1)) 
				rmsd_after= self.align_by_bb(inputPose_tmp, poses_by_chain[i])
				print "RMSD after",rmsd_after
				if(rmsd_after > 0.25):
					print "Error, the RMSD is too large~~~!!!"
					assert(0==1)
				pose_array.append(inputPose_tmp.clone())
				##pose_array[len(pose_array)-1].dump_pdb("%s/testSS_before_ideal_%02d.pdb"%(general_out_dir,len(pose_array)-1))
		else:
			print "Error, I shouldnt be here"
		
		#Output disjointed SSs that will be used as targets for idealization
		for i in range(len(pose_array)):
			tmp_out_name="%s/testSS_before_ideal_%02d.pdb"%(stage1_tmp_dir, i)
			pose_array[i].dump_pdb(tmp_out_name)
			print "Dumped disjointed SS input target structure to: ", tmp_out_name
			print "Done"
			
		return pose_array
		
	def idealize_ss_array(self,
		pose_array, 
		context_pose_centroid,
		b_has_context_pose,
		clustersDB=None,
		ss_fragment_population_threshold=30,  #100 is default "quite flexible maybe", 10 is allow almost all, 200 is for very ideal structs
		fragment_designability_ratio=0.5,     #0.8 is default "be VERY very carefull when modifing this value", < will afect your ability to design
		rmsd_tolerance=0.8,                   #0.8 is default, < is OK (more strict)
		rmsd_sliding_window_size=8,           #8 is default The shorter, the nearer to input, the larger  the more idealized struct will come out
		stage1_tmp_dir="./"):
		
		#Loop them all!
		print "Idealizing SS elements array"
		idealized_pose_array=[]
		##for i in range(4,5):
		for i in range(len(pose_array)):
			print "Idealizing SS fragment %d of %d" % ( (i+1), len(pose_array) )
			#Generate a context pose for optimization of non-local interactions:
			b_is_first=True
			tmp_context_pose = self.rosetta.pose_from_sequence("")
			if b_has_context_pose:
				print "Note: Clashes and scores will be calculated in the presence of the context protein/object"
				tmp_context_pose = context_pose_centroid.clone()
				b_is_first=False
			for j in range(len(idealized_pose_array)):
				if (j != i):
					if(b_is_first):
						tmp_context_pose = idealized_pose_array[j].clone()
						b_is_first=False
					else:
						tmp_context_pose.append_pose_by_jump(pose_array[j],1)
			for j in range(len(idealized_pose_array), len(pose_array)):
				if (j != i):
					if(b_is_first):
						tmp_context_pose = pose_array[j].clone()
						b_is_first=False
					else:
						tmp_context_pose.append_pose_by_jump(pose_array[j],1)
			##debug out pdb
			##tmp_context_pose.dump_pdb("%s/test_SScontext_%02d.pdb"%(general_out_dir,i))
			
			#Call idealizer
			##debug out pdb
			##pose_array[i].dump_pdb("%s/test_before_ideal_%02d.pdb"%(general_out_dir,i))
			##tmp_context_pose.dump_pdb("%s/test_context_%02d.pdb"%(general_out_dir,i))
			#Idealize SS, (12a.a. rmsd spans works fine, use less for close matching)
			result_pose = self.idealize_ss( inputPose=pose_array[i], 
										contextPoses=pose_array, 
										init_res=1, 
										end_res=pose_array[i].total_residue(), 
										context_pose=tmp_context_pose, 
										scorefxn=self.scorefxn_remodel_cen,
										clustersDB=clustersDB,
										fragment_population_threshold=ss_fragment_population_threshold, 
										rmsd_test_window=rmsd_sliding_window_size,
										rmsd_acceptance_tolerance_ratio=rmsd_tolerance, 
										fragment_population_acceptance_ratio=fragment_designability_ratio,
										debug=True)
			idealized_pose_array.append(result_pose.clone())
			#tmp file output
			outfilename="%s/test_idealSS_%02d.pdb"%(stage1_tmp_dir,i)
			result_pose.dump_pdb(outfilename)
			print "Dumped  idealized SS structure to: ", outfilename
			
			print "Done idealizing span ", 1, "-", result_pose.total_residue(), " of segment: ", i+1
		print "Done idealizing the SS elements array"
		return idealized_pose_array
		
	def idealize_bond_distances(self, pose_for_idea):
		#Fist clear the constraint set
		pose_for_idea.remove_constraints()
		##scorefxn_rel = rosetta.core.scoring.ScoreFunctionFactory.create_score_function("talaris2013")
		IdeBondMov=self.IdealBondMover()
		IdeBondMov.chainbreaks(True)
		IdeBondMov.report_CA_rmsd(True)
		IdeBondMov.apply(pose_for_idea)
		pose_for_idea.remove_constraints()
		return (IdeBondMov.get_last_move_status() == self.rosetta.MS_SUCCESS)
	
		###%%timeit
	#TODO: Add fast clash checking after a "good angle" test pass
	def idealize_ss(self, 
					inputPose, 
					contextPoses, 
					init_res, 
					end_res, 
					context_pose, 
					scorefxn=None,
					clustersDB=None,
					fragment_population_threshold = 100, 
					rmsd_test_window=4, 
					rmsd_acceptance_tolerance_ratio=1.0, 
					fragment_population_acceptance_ratio=0.8, 
					debug=False):
		print "Idealizing using RMSD restrain a pose of #res = ", inputPose.total_residue(), ", span: ", init_res, end_res
		
		inputPose =inputPose.clone()
		targetPose =inputPose.clone()
		refPose =inputPose.clone()
		contPose = context_pose.clone()
		#Create a Hash of the context
		ss_context_pose_hash = self.rosetta.core.pose.xyzStripeHashPose(contPose, self.rosetta.core.pose.PoseCoordPickMode.N_CA_C_CB)
		
		
		#The test window cannot excedd the pose's size
		if (inputPose.total_residue() < (rmsd_test_window*2)):
			#This rounds to the lower number
			rmsd_test_window = int(inputPose.total_residue()/2)

		
		#Debugoutput PDB
		#refPose.dump_pdb("%s/test_reference.pdb"%general_out_dir)
		firstLoopRun = True
		lastBestCluster=-1
		sequence_jump_size=3
		for testres in range(max(1,init_res), min(end_res, refPose.total_residue()), sequence_jump_size): 
		#for testres in range(1, 10): 
			#Compute the span for a 4mer
			if (testres+3 <= refPose.total_residue()):
				testresndx=[testres,testres+3]
			else:
				testresndx=[testres, refPose.total_residue()]
			if debug:
				print "Idealizing residues in span: ", testresndx
			rmsd_result=np.zeros((1),float)
			rmsd_threshold=0.6
			
			number_of_matches_treshold = 3
			good_fragments_ndx=[]
			if firstLoopRun:
				rmsd_result=clustersDB.rmsd_calculator.get_broadcast_coordinate_rmsd(
						np.expand_dims( pose_res_to_array(refPose, testresndx, atoms=["CA","C","O","N"]), 0), 
						clustersDB.cluster_lr_centers_coordinates.copy() 
						).ravel()
				
				#Calculate those clusters that match this fragment
				while ( len(good_fragments_ndx) < number_of_matches_treshold  ):       
					good_fragments_ndx = np.where (rmsd_result < rmsd_threshold)[0]
					good_fragments_ndx = good_fragments_ndx[np.where(
													clustersDB.clusters_lr_population[good_fragments_ndx] > fragment_population_threshold)]
					rmsd_threshold += 0.01
				rmsd_threshold-=0.01
				if debug:
					print "Auto rmsd threshold: ", rmsd_threshold, good_fragments_ndx
			else:
				lasBestPositions = clustersDB.correspondence_dic_HR_LR[lastBestCluster]
				#ToDo: Is index [1] correct?
				good_fragments_ndx = np.where (clustersDB.clusters_lr_TCM[lasBestPositions] > 0)[1]
				
				good_fragments_ndx = good_fragments_ndx[np.where(
													clustersDB.clusters_lr_population[good_fragments_ndx] > fragment_population_threshold)]
				if debug:
					print "Using TCM with last best HR cluster: ", lastBestCluster, "Corresponding to #LR cluster(s):", len(lasBestPositions)
					print "Following the LR TCM, #clusters: ", len(good_fragments_ndx) 
					print "After prunning by population #clusters: ", len(good_fragments_ndx)
		
			
			#Use HR clustering for replacement, careful, this could introduce a bug
			hr_good_fragments_ndx=[]
			max_num_fragments_threshold = 5000
			while ( (len(hr_good_fragments_ndx) < number_of_matches_treshold) or ( 
					len(hr_good_fragments_ndx) > max_num_fragments_threshold ) ):
				hr_good_fragments_ndx=[]
				for i in good_fragments_ndx:
					for j in clustersDB.correspondence_dic_LR_HR[i]:
						hr_good_fragments_ndx.append(j)
				hr_good_fragments_ndx = np.asarray(hr_good_fragments_ndx)
				#Previously mentioned bug may happen here
				hr_good_fragments_ndx = hr_good_fragments_ndx[np.where(
										clustersDB.clusters_hr_population[hr_good_fragments_ndx] > fragment_population_threshold)]
				
				#blahh... stupid mistake before, dont waste computer work Daniel :)
				#brief: only test unique fragments
				hr_good_fragments_ndx=np.unique(hr_good_fragments_ndx)
			
				if (len(hr_good_fragments_ndx) < number_of_matches_treshold):
					fragment_population_threshold -= 1
				elif (len(hr_good_fragments_ndx) > max_num_fragments_threshold ):
					hr_good_fragments_ndx = hr_good_fragments_ndx[np.argsort(
												clustersDB.clusters_hr_population[hr_good_fragments_ndx])[-(max_num_fragments_threshold+1):-1]]
					fragment_population_threshold=min(clustersDB.clusters_hr_population[hr_good_fragments_ndx])
			if debug:
				print "Now going to high resolution sampling using #clusters: %d, with a population threshold of: %d"%(
																	len(hr_good_fragments_ndx), 
																	fragment_population_threshold)
			
			#HACK magic number, remind to remove it
			min_rmsd_1=9999
			min_rmsd_2=9999
			#HACK magic number, remind to remove it
			best_score=9999
			best_pop_score=0
			best_hr_cluster=-1
			
			fragOverlaps4Check=[]
			if(testresndx[1] > 4):
				for indx in range(testresndx[0]-3, testresndx[1]-3):
					fragOverlaps4Check.append([indx, indx+3]) 
			#Debug    
			##print len(hr_good_fragments_ndx), hr_good_fragments_ndx
			##hr_good_fragments_ndx=np.unique(hr_good_fragments_ndx)
			##print len(hr_good_fragments_ndx), hr_good_fragments_ndx
			#End debug
			
			
			bestTargetPose=self.rosetta.pose_from_sequence("")
			best_rmsd=-9999
			
			for cluster_hr_number in hr_good_fragments_ndx:
				counter=0
				for angle in clustersDB.cluster_hr_centers_store_data["bb"] [cluster_hr_number]:
					if( (testresndx[0]+counter) > testresndx[1]):
						break
					targetPose.set_phi(testresndx[0]+counter, angle[0])
					targetPose.set_psi(testresndx[0]+counter, angle[1])
					##if( (testresndx[0]+counter) < testresndx[1]):
					targetPose.set_omega(testresndx[0]+counter, angle[2])
					counter+=1
				   
				#rmsd = rmsd_bb(targetPose_4rmsd, refPose_4rmsd)
				#Local RMSD
				rmsd1 = self.align_bb_by_fragment(targetPose, inputPose, max(1, testresndx[0]-rmsd_test_window), 
											min(targetPose.total_residue(), testresndx[1]+rmsd_test_window))
				#Global RMSD
				rmsd2 = self.align_bb_by_fragment(targetPose, inputPose, 1, inputPose.total_residue())
				#rmsd = self.align_by_bb(targetPose_4rmsd, refPose_4rmsd)
				#If the RMSD improves
				if ( (rmsd1 <= (min_rmsd_1*rmsd_acceptance_tolerance_ratio) ) 
					and (rmsd2 <= (min_rmsd_2*rmsd_acceptance_tolerance_ratio)) ):
					##and If there is no clash
					##if (not check_bb_clashes_hash_vs_pose(ss_context_pose_hash, targetPose_4rmsd)):
						##and if the score improves (or does not unninprove)
							
						tmp_pose_for_score = contPose.clone()
						scorefxn(tmp_pose_for_score)
						tmp_pose_for_score.append_pose_by_jump( targetPose, 1 )
						this_score=scorefxn( tmp_pose_for_score )
						#HACK: +2.5 is about 1aa in rosetta units
						if ( this_score <= ( best_score+2.5 )):
							#Check fragment overlaps
							isGoodOverlap=True
							for fragOverlapsNdx in fragOverlaps4Check:
								min_rmsd_result=clustersDB.rmsd_calculator.get_broadcast_coordinate_rmsd(
										np.expand_dims( pose_res_to_array(targetPose, 
														fragOverlapsNdx, atoms=["CA","C","O","N"]), 0
													   ), 
										clustersDB.cluster_lr_centers_coordinates.copy() 
									).ravel()
								min_rmsd_result_ndx=min_rmsd_result.argmin()
								if (min_rmsd_result[min_rmsd_result_ndx] > 
									clustersDB.cluster_lr_centers_store_data["threshold_distance"][min_rmsd_result_ndx] * 1.0 ):
									##print "UPPPPSSS, this overlap is wrong: ", targetPose.total_residue(), fragOverlapsNdx, \
									##    min_rmsd_result[min_rmsd_result_ndx], \
									##    self.cluster_lr_centers_store_data["threshold_distance"][min_rmsd_result_ndx], \
									##    " :: ", cluster_hr_number, rmsd, this_score, "Skipping"
									isGoodOverlap=False
									break
							if(isGoodOverlap):
								#Check the minimum population of fragments
								#ToDo: check this function
								##score_population = pose_fragment_population_within(targetPose, testresndx[0], 
								##                        testresndx[1], cut_off=0.51).min()
								score_population = self.pose_fragment_population_within(targetPose, 1, 
														testresndx[1], cut_off=0.51, clustersDB=clustersDB).min()
								
								if ( score_population >= (best_pop_score*fragment_population_acceptance_ratio)):  
									min_rmsd_1 = rmsd1
									min_rmsd_2 = rmsd2
									best_pop_score=score_population
									best_hr_cluster = cluster_hr_number
									best_score = this_score
									bestTargetPose=targetPose.clone()
									best_rmsd=rmsd2
								##targetPose.dump_pdb("%s/testbuild_target_%d.pdb"%(general_out_dir, cluster_hr_number))
								##inputPose.dump_pdb("%s/testbuild_ref_%d.pdb"%(general_out_dir, cluster_hr_number))
							
					##else:
					##    targetPose_4rmsd.dump_pdb("%s/testbuild.pdb"%general_out_dir)
					##    assert(0==1)
					##if (min_rmsd <= 0.3):
					##    print "Target RMSD reached with:", min_rmsd, ". Continuing" 
					##    break
				
			#If there is any solution
			if (best_hr_cluster > -1):
				if debug:
					print ("Best cluster :", best_hr_cluster, "RMSDs: ", min_rmsd_1, min_rmsd_2, "Cen_score: ", best_score, 
					   "Population: ", best_pop_score)
				lastBestCluster = best_hr_cluster
			else:
				#DIE because there is no solution!
				if debug:
					print "I couldn't find a solution for this fragment"
				assert(0==1)
			firstLoopRun = False 
			
			#if (min_rmsd > rmsd_threshold ):
				#print "Failed min RMSD target to original pose: ", rmsd_threshold
				#assert (True == False)
			
			
			#Apply the angles of the best cluster to the targetPose
			targetPose=bestTargetPose.clone()
		"""
		#Aligns the resulting idealized structure aligned to the original query coordinates
		rmsdVal=self.align_by_bb(targetPose, inputPose)
		"""
		
		print "Final total-RMSD after idealization= ", best_rmsd
		
		#Return the result
		return(targetPose)
		
		
	def build_loops(self,
					idealized_pose_array,
					context_pose_centroid,
					b_has_context_pose,
					target_num_loops=100, #Stop condition
					sufficient_num_loops=100, #Suffience condition to succed
					sufficient_num_loops_for_each_len=40, #Suffience condition to succed per loop len
					min_loop_size=0, #Min desired loop size
					max_loop_size=10, #Max desired loop size
					max_num_frag_matches = 50,  #N-side seeds number of neighbors
					max_resolution_frag_matches = 0.6,  #N-side seeds max matching res
					clustersDB=None,
					out_frag_tmp_dir="possible_loops_tmp"):
		#ToDo: Implement various SS delta trimed tests (need to have independent paths for each)
		##%run
		##%%timeit
		##Exhaustive sampler of SS loop connectors
		b_debug=False
		
		#Don't change
		min_loop_size+=2
		max_loop_size+=2
		#

		#RMSD thresholds
		rmsd_treshold_lr= 1.0   #0.8 good, 1.0 OK
		tip_dist_treshold_lr=rmsd_treshold_lr*1.5

		#TCMs: used the same cutoff on both!! This is only for the beggining
		this_sampling_TCM_ndx = clustersDB.clusters_lr_TCM_ndx_p[10]      
		this_sampling_TCM_len = clustersDB.clusters_lr_TCM_ndx_len_p[10]  

		#Allowed degrees that omega is permited to deviate from planarity
		allowed_omega_deviation=30  #30 OK
		#max allowed fragment overlap rmsd from center
		max_allowed_fragment_overlap_rmsd_deviation=0.6  #0.7 seems OK
		population_score_max_rmsd=0.8   #score populations up to a distance of, 0.8 seems OK
		allowed_A_fragment_delta=0  #num a.a.
		allowed_B_fragment_delta=0  #num a.a.

		#Calculate the lenght of the context SS
		total_SS_residues_minus_ends=0
		for ideal_ss_pose in idealized_pose_array:
			total_SS_residues_minus_ends += ideal_ss_pose.total_residue()-2
		##print total_SS_residues_minus_ends
		#Create a new pose to hold the context
		tmpline=""
		for i in range(total_SS_residues_minus_ends):
			tmpline+=self.fake_aminoacid_letter_for_design
		all_SS_pose = self.rosetta.pose_from_sequence(tmpline)
		#Copy the fragments
		tmp_ndx=1
		for ideal_ss_pose in idealized_pose_array:
			#Hack: copy the fragment but not the first or last residue
			#ToDo: Remove it.
			all_SS_pose.copy_segment((ideal_ss_pose.total_residue()-2),ideal_ss_pose,tmp_ndx,2)
			tmp_ndx+=(ideal_ss_pose.total_residue()-2)
			
		#all_SS_pose_wContext=all_SS_pose.clone()
		if b_has_context_pose:
			print "Note: Clashes and scores for loops will be calculated in the presence of the context protein/object"
			all_SS_pose.append_pose_by_jump(context_pose_centroid, 1)
			##all_SS_pose.dump_pdb("%s/testSSandContext1.pdb"%general_out_dir)
		#Hash the ss context
		context_pose_hash = self.rosetta.core.pose.xyzStripeHashPose(
			all_SS_pose,
			self.rosetta.core.pose.PoseCoordPickMode.N_CA_C_CB)  #was rosetta.core.pose.PoseCoordPickMode.BB)

		good_loop_poses = []
		for ideal_ss_pose_i in range( len(idealized_pose_array)-1):
		###for ideal_ss_pose_i in range( 2, 3 ):
			##rmsd_calculator = rmsd_calc(num_threads_global)
			print "Idealizing connection # %d of %d" % ( (ideal_ss_pose_i+1), (len(idealized_pose_array)-1) )
			
			queryPose1 = idealized_pose_array[ideal_ss_pose_i]
			queryPose2 = idealized_pose_array[ideal_ss_pose_i+1]
			
			#Generate a context pose for optimization of non-local interactions:
			b_is_first=True
			tmp_context_pose_for_minimization = self.rosetta.pose_from_sequence("")
			for j in range(len(idealized_pose_array)):
				if( (j != ideal_ss_pose_i) and  (j != (ideal_ss_pose_i+1)) ):
					if(b_is_first):
						tmp_context_pose_for_minimization = idealized_pose_array[j].clone()
						b_is_first=False
					else:
						tmp_context_pose_for_minimization.append_pose_by_jump(idealized_pose_array[j],1)
			
			#Quickly estimate the minimum necessary number of a.a. to close the loop
			"""BMC Structural Biology 2013, 13(Suppl 1):S5  doi:10.1186/1472-6807-13-S1-S5
			The expected length was calculated as 3.8A x(n + 1), 
			where n is the number of the amino acids on the loop and 3.8A
			is the average distance between two amino acids."""
			fragment_np_dist = distance_between_two_Ca(queryPose1, queryPose1.total_residue(), queryPose2, 1)
			#the -1 is just in case
			num_amino_guess=int(fragment_np_dist/3.8)-1
			print "Rough guess: In order to connect these segments I'll need a minimum # of aminoacids = " ,  num_amino_guess
			
			#Calculate the last nearest fragments to queryPose1 C-term
			#ToDO: URGENT!!! NEED TO CHANGE THIS DOWN TO ACOUNT FOR allowed_A_fragment_delta and allowed_B_fragment_delta
			rmsd_result=clustersDB.rmsd_calculator.get_broadcast_coordinate_rmsd(
							np.expand_dims( pose_res_to_array(queryPose1, ([queryPose1.total_residue()-(4-1), 
											queryPose1.total_residue()]), atoms=["CA","C","O","N"]), 0), 
											clustersDB.cluster_lr_centers_coordinates.copy() ).ravel()
			fragment_match_cut_off=rmsd_result.min()
			print fragment_match_cut_off
			p1_C_fragments_ndx = []
			while ( (len(p1_C_fragments_ndx)) < max_num_frag_matches 
				   and (fragment_match_cut_off <= max_resolution_frag_matches) ):
				p1_C_fragments_ndx = np.where(rmsd_result <= fragment_match_cut_off)[0]
				fragment_match_cut_off+=0.01
			print "F1-C resolution, fragments: ", (fragment_match_cut_off-0.01), p1_C_fragments_ndx
			
			#make_list_of_post_transition_fragments
			post_fragments_ndx_lr=[]
			for indx in p1_C_fragments_ndx:
				for jndx in this_sampling_TCM_ndx[indx]:
					if( this_sampling_TCM_len[jndx] > 0):
						post_fragments_ndx_lr.append(jndx) 
			##post_fragments_ndx_lr = list(itertools.chain.from_iterable(post_fragments_ndx_lr))
			##post_fragments_ndx_lr = list(OrderedDict.fromkeys(post_fragments_ndx_lr))
			print "Number of post transition fragments: ", len(post_fragments_ndx_lr)
			
			min_insertion_len = min_loop_size  #( max( num_amino_guess, min_loop_size ) )
			max_insertion_len =( max( max_loop_size, (min_insertion_len) ) )
			##b_found_solution=False
			lr_good_poses=[]
			#Initialize the lr_combination_array
			fragment_lr_combinations_array=[]
			for ifragmentndx in post_fragments_ndx_lr:
				fragment_lr_combinations_array.append([ifragmentndx])
			
			for window_len in range(min_insertion_len, max_insertion_len+1):
				this_window_len_loop_num_solutions=0
				target_num_aa=window_len
				##print target_num_aa ##(+2 because roseta makes stupid things with the BB angles of the first and last residue)
				##break
				tmpline=""
				for i in range(0, target_num_aa+2):
					tmpline+=self.fake_aminoacid_letter_for_design
				tmp_loop_pose = self.rosetta.pose_from_sequence(tmpline)
				#make another copy with the correct len (CAUTION can break if the loop len is < 2 ??)
				tmpline=""
				##Be carefull with the next for, it can break it all if it's size == 0 !!! (can happen!)
				for i in range(0, target_num_aa):
					tmpline+=self.fake_aminoacid_letter_for_design
				tmp_loop_pose_wcorrect_len = self.rosetta.pose_from_sequence(tmpline)
				#make another copy with the correct len-2 (CAUTION this will break if the loop len is =< 2 ??)
				tmpline=""
				##Be carefull with the next for, it can break it all if it's size == 0 !!! (can happen!)
				for i in range(0, target_num_aa-2):
					tmpline+=self.fake_aminoacid_letter_for_design
				tmp_loop_pose_minus_ends = self.rosetta.pose_from_sequence(tmpline)
				##break
				#Calculate the combinations
				needed_depth=int(np.ceil((window_len-2)/2.0))
				#if b_debug:
				print "Testing de-novo loops of len: ", window_len-2
				 
				
				#New combinatios code. ToDo: IMprove writting style
				#It uses tooo much memory
				tmp_current_TCM_deepth=len(fragment_lr_combinations_array[0])
				print "Current TCM Deepth Sampling Calcultated vs Needed: ", tmp_current_TCM_deepth, needed_depth

				##Experimental code to ramp the TCM cut-off
				if (needed_depth == 1):
					print "Switching TCM to level 1"
					this_sampling_TCM_ndx = clustersDB.clusters_lr_TCM_ndx_p[10]      
					this_sampling_TCM_len = clustersDB.clusters_lr_TCM_ndx_len_p[10]  
				elif(needed_depth == 2):
					print "Switching TCM to level 2"
					this_sampling_TCM_ndx = clustersDB.clusters_lr_TCM_ndx_p[20]      
					this_sampling_TCM_len = clustersDB.clusters_lr_TCM_ndx_len_p[20]  
				elif(needed_depth == 3):
					print "Switching TCM to level 3"
					this_sampling_TCM_ndx = clustersDB.clusters_lr_TCM_ndx_p[40]      
					this_sampling_TCM_len = clustersDB.clusters_lr_TCM_ndx_len_p[40]  
				elif(needed_depth >= 4):
					print "Switching TCM to level 4"
					this_sampling_TCM_ndx = clustersDB.clusters_lr_TCM_ndx_p[80]      
					this_sampling_TCM_len = clustersDB.clusters_lr_TCM_ndx_len_p[80]  
				## END Experimental code
				
				if (tmp_current_TCM_deepth < needed_depth):
					curr_combinations_array=fragment_lr_combinations_array[:]
					fragment_lr_combinations_array=[]
					for combinations_buffer in curr_combinations_array:
						if(this_sampling_TCM_len[combinations_buffer[-1]] > 0):  #(>0):
							start_ndxs=this_sampling_TCM_ndx[combinations_buffer[-1]]
							generate_fragment_combination_list_by_TCM( start_ndxs, 
														  tmp_current_TCM_deepth, 
														  needed_depth, 
														  this_sampling_TCM_ndx, 
														  this_sampling_TCM_len,
														  combinations_buffer, 
														  fragment_lr_combinations_array )
					if b_debug:
						print "Sample FPaths(3): ", fragment_lr_combinations_array[0:3]
						
				#print   fragment_lr_combinations_array
				
				##if b_debug:
				print "Testing # LR combinations by walking the TCM: ", len(fragment_lr_combinations_array), ":)"
				#New code END
				
				#ToDo It should be possible to test all the combinations of the same depth at once.
				#Loop to test the Fragment combinations in LR
				
				##possible_lr_fragment_paths=[]
				#lr_good_poses=[]
				
				#Hack : In case that we want to break early do a Random shuffle of the testing order list, 
				#because near->by_indexes give place to near-like solutions
				np.random.shuffle(fragment_lr_combinations_array)
				
				## Do LimitNumOfSolutionsToTest?
				##if (len(fragment_lr_combinations_array) > 1e5):
				##    tmp_new_paths_num=int(len(fragment_lr_combinations_array)*0.1)
				##    print "Limiting the number of solutions to test to 10%:"
				##    fragment_lr_combinations_array=fragment_lr_combinations_array[:tmp_new_paths_num]
				##    print "New size of LR combinations: ", len(fragment_lr_combinations_array), ":))!"

				#ToDo: Make this much faster
				tmp_counter=-1
				for lr_fragment_path in fragment_lr_combinations_array[:]:
					tmp_counter+=1
					curr_res_ndx=2
					tmp_rmsd= 9999.99
					for lr_fragment in lr_fragment_path:
						for angle in clustersDB.cluster_lr_centers_store_data["bb"] [lr_fragment]:
							if(curr_res_ndx > (target_num_aa+1)):
								break
							tmp_loop_pose.set_phi(curr_res_ndx, angle[0])
							tmp_loop_pose.set_psi(curr_res_ndx, angle[1])
							tmp_loop_pose.set_omega(curr_res_ndx, angle[2])
							curr_res_ndx+=1
						curr_res_ndx-=1
						
					#Experimental
					#Try to check +/- 2 residues in each extreme at the same time:
					tmp_ndx1a=2
					tmp_ndx1b=(target_num_aa+1)
					tmp_ndx2=0
					tmp_ndx3=0
					best_rmsd=999999
					loop_np_distance=distance_between_two_Ca(tmp_loop_pose, 2, tmp_loop_pose, tmp_loop_pose.total_residue()-1)
					queryPoseASize=queryPose1.total_residue()
					queryPoseBSize=queryPose2.total_residue()
					
					#Here starts the "loop" checking "the p loops" ::))
					
					for delta_resA in range(allowed_A_fragment_delta+1):
						if ( (queryPoseASize-delta_resA) < 3 ):
							#The A framgnet is too small
							break
						for delta_resB in range(allowed_B_fragment_delta+1):
							if ( (queryPoseBSize-delta_resB) < 3 ):
								#The B framgnet is too small
								break
							tmp_ndx2=queryPoseASize-delta_resA
							tmp_ndx3=1+delta_resB
							
							fragment_np_dist = distance_between_two_Ca(queryPose1, tmp_ndx2, queryPose2, tmp_ndx3)
							if (abs(fragment_np_dist-loop_np_distance) > tip_dist_treshold_lr):
								continue
							tmp_rmsd, tmp_max_pair_dist = self.align_bb_1pose_to_2poses_by_ndxs_tips(tmp_loop_pose, 
																								[tmp_ndx1a], 
																								[tmp_ndx1b], 
																								queryPose1, 
																								[tmp_ndx2], 
																								queryPose2, 
																								[tmp_ndx3])
									
							#End Experimental
			
						
							if ( (tmp_rmsd < rmsd_treshold_lr) and (tmp_max_pair_dist < tip_dist_treshold_lr) ):
								#print "Possible match LR rmsd/lr-path: ", tmp_rmsd, lr_fragment_path, tmp_counter
								#Check for clashes (Hacked to remove first and last residue)
								##tmp_loop_pose_wcorrect_len.copy_segment((tmp_loop_pose.total_residue()-2),tmp_loop_pose,1,2)
							
								#Test clashes
								##hasBBclashes=False
								if ( (target_num_aa-2) > 0):
									tmp_loop_pose_minus_ends.copy_segment((tmp_loop_pose.total_residue()-4),tmp_loop_pose,1,3)
									hasBBclashes = self.check_bb_clashes_hash_vs_pose(context_pose_hash, tmp_loop_pose_minus_ends)
									if hasBBclashes:
										continue
								#If the cart constraint is to high you'll be in trouble when checking OMEGA
								#Works best with sofrep, why? Ask Frank one of this days.
								test_mini_pose, score_merge, score_cen = self.connect_3_poses_byNDX_and_minize_w_constraints(
																				queryPose1,  
																				tmp_loop_pose, 
																				queryPose2, 
																				tmp_ndx2, 
																				tmp_ndx3, 
																				tmp_context_pose_for_minimization,  
																				sfx=self.scorefxn_soft_rep,
																				harmonic_constraint_streght=0.6, 
																				cart_constraint_weight=0.8)
								
								#Also test clashes after minimization
								##hasBBclashes=False
								if ( (target_num_aa-2) > 0):
									tmp_loop_pose_minus_ends.copy_segment((tmp_loop_pose.total_residue()-4),tmp_loop_pose,1,3)
									hasBBclashes = self.check_bb_clashes_hash_vs_pose(context_pose_hash, tmp_loop_pose_minus_ends)
									if hasBBclashes:
										continue
								
								#Experimental code (Omega angle, we have to be a little permisive because this will undergo further minimization),
								#Usual omega values are abs(-160 to 180 and 0 to 20), we will allow  150 to 180 and 0 to 30
								#ToDo: Speed up this part
								#ToDo, current is a Hack: remove, maybe we need code to check omega overlaps beforehand???
								tmp_is_good_overlap=True
								low_bound=max(1, tmp_ndx2-4) #queryPose1.total_residue()-4)
								high_bound=min(test_mini_pose.total_residue(), tmp_ndx2-1+window_len+4) #queryPose1.total_residue()-1+window_len+4)
								for ires in range(low_bound, high_bound+1):
									tmp_ome=abs(test_mini_pose.omega(ires))
									if((tmp_ome > (0+allowed_omega_deviation)) and (tmp_ome < (180-allowed_omega_deviation))):
										tmp_is_good_overlap=False
										if b_debug:
											print "OMega angle check fail: res/omega", ires, tmp_ome, tmp_counter
										break
								#End experimental code
								
								
								#ToDo: Add carefull checks for:
								#A) Clashes 
								#B) Angles/vs/SS
								#tmp_observerd=0.0
								#tmp_max_expected=0.0
								connection_population_score_min=0.0
								connection_population_score_max=0.0
								connection_population_score_avg=0.0
								if(tmp_is_good_overlap):
									tmp_is_good_overlap = self.check_span_for_ideal_fragments(test_mini_pose, 
																	low_bound, 
																	high_bound,
																	rmsd_pert=1.10, 
																	max_rmsd_limit=max_allowed_fragment_overlap_rmsd_deviation,
																	clustersDB=clustersDB)
								else:
									continue
									
								if(tmp_is_good_overlap):
									#Score this connection
									[connection_population_score_min, 
									 connection_population_score_max,
									 connection_population_score_avg] = self.score_region_conections(test_mini_pose, 
												low_bound,
												high_bound,
												des_assignment_cut_off=population_score_max_rmsd,
												clustersDB=clustersDB,
												b_is_permissive=True)
									
									if (connection_population_score_min<=0):
										tmp_is_good_overlap=False
									
									if b_debug:
										print ("Overlap TESTED: ", tmp_is_good_overlap,
												connection_population_score_min, 
												connection_population_score_max,
												connection_population_score_avg) #tmp_max_expected, 
								else:
									continue
									
								if (tmp_is_good_overlap):
									#Score this solution
									if b_debug:
										tmp_loop_pose.dump_pdb("%s/testmini_%02d_%04d_ll%02d_a%08d_ori.pdb" % (
															 out_frag_tmp_dir, ideal_ss_pose_i, (len(lr_good_poses)), window_len-2, tmp_counter) )
									
									tmp_loop_pose.copy_segment(tmp_loop_pose.total_residue(),
															   test_mini_pose,
															   1,(tmp_ndx2-1)) #queryPose1.total_residue()-1))
									
									tmp_pose_for_score = all_SS_pose.clone()
									if ( (target_num_aa-2) > 0 ):
										tmp_loop_pose_minus_ends.copy_segment((tmp_loop_pose.total_residue()-4),tmp_loop_pose,1,3)
										tmp_pose_for_score.append_pose_by_jump(tmp_loop_pose_minus_ends, 1)
										
									#Score loop
									this_loop_score=self.scorefxn_remodel_cen(tmp_pose_for_score)
									#Loop scoring function ()
									connection_population_score=( (connection_population_score_avg*connection_population_score_avg)
																 /(connection_population_score_avg-connection_population_score_min) )
									print "DEB: ", connection_population_score, connection_population_score_avg, connection_population_score_min
									lr_good_poses.append([tmp_rmsd, 
														  tmp_loop_pose.clone(), 
														  test_mini_pose.clone(),   #tmp_loop_pose_minus_ends.clone(), 
														  this_loop_score, 
														  connection_population_score,
														  delta_resA,
														  delta_resB]) 
									
									this_window_len_loop_num_solutions+=1
									#Debug
									#if True:
									print "TEST NUM SOLUTIONS THIS LEN LOOP: ", this_window_len_loop_num_solutions
									print ("Yes!!! Match LR rmsd/lr-path found: ", delta_resA, delta_resB, len(lr_good_poses), lr_fragment_path, tmp_counter, tmp_rmsd, 
											 tmp_max_pair_dist)
									#tmp_loop_pose.dump_pdb("%s/tmp_fragments/testloop_%02d_%04d.pdb" % (general_out_dir, ideal_ss_pose_i, tmp_counter) )
									test_mini_pose.dump_pdb("%s/testmini_%02d_%04d_ll%02d_a%08d.pdb" % (
															 out_frag_tmp_dir, ideal_ss_pose_i, (len(lr_good_poses)-1), window_len-2, tmp_counter) )
									#tmp_context_pose_for_minimization.dump_pdb("%s/tmp_fragments/testtest_%02d_context_%04d.pdb" % (
									#                                        general_out_dir, ideal_ss_pose_i, tmp_counter) )
									#End Debug
									
									#Turn the flag of sol found!!
									##b_found_solution=True
									
									##hack! Remove
									if( len(lr_good_poses) >= target_num_loops):
										print "Hack break 1a!!! :()"
										break
									elif(this_window_len_loop_num_solutions >= sufficient_num_loops_for_each_len):
										print "Hack break 1b!!! :()"
										break
								if( len(lr_good_poses) >= target_num_loops):
									print "Hack break 2a!!! :()"
									break
								elif(this_window_len_loop_num_solutions >= sufficient_num_loops_for_each_len):
									print "Hack break 2b!!! :()"
									break
							if( len(lr_good_poses) >= target_num_loops):
								print "Hack break 3a!!! :()"
								break
							elif(this_window_len_loop_num_solutions >= sufficient_num_loops_for_each_len):
								print "Hack break 3b!!! :()"
								break
						if( len(lr_good_poses) >= target_num_loops):
							print "Hack break 4a!!! :()"
							break
						elif(this_window_len_loop_num_solutions >= sufficient_num_loops_for_each_len):
							print "Hack break 4b!!! :()"
							break
					if( len(lr_good_poses) >= target_num_loops):
						print "Hack break 5a!!! :()"
						break
					elif(this_window_len_loop_num_solutions >= sufficient_num_loops_for_each_len):
						print "Hack break 5b!!! :()"
						break
							#Debug
							#else:
								#print "This connection cannot happen"
							
							#assert(0==1)
							#End Debug
							
				#Check if we have solutions and append to the result's array
				print "Until now I have found ", len(lr_good_poses), "solutions."
			
				if( len(lr_good_poses) >= sufficient_num_loops):
					print "... and Finally I am saying enough due to the sufficient loops parameter! I have found ", len(lr_good_poses), "solutions."
					break
					
				#if (len(lr_good_poses) >= target_num_loops):
				#    print "and Finally, I have found ", len(lr_good_poses), "solutions."
				#    break
					
					
			#assert(0==1)
			#Stop if it can't find enough solution for some of the loop connections
			if ( len(lr_good_poses) > 0 ): ##b_found_solution:
				good_loop_poses.append(lr_good_poses)
			else:
				print ("ERROR, I couldn't find any suitable solutions (%d)to connect the SS elements: %d<-->%d" % 
					   (len(lr_good_poses), (ideal_ss_pose_i+1), (ideal_ss_pose_i+2) ) )
				print "I'll stop now, sorry"
				break
		print "Done Finding Loops, go ahead!"
		
		return good_loop_poses
		
		
	#Determine "perfection" of fragments
	def check_fragment_distances(self, result_idealized_pose, treshold_factor=1.10, verbose=False):
		
		##Debugline
		result_idealized_pose=result_idealized_pose.clone()
		##End debugline
		
		
		pymol_line="select resid "
		
		print "Determining SS probability using fragments "
		
		aaProbabilies=[]
		for testres in range(0, result_idealized_pose.total_residue()):
			aaProbabilies.append([])
		
		has_bad_fragments=False
		for testres in range(1, result_idealized_pose.total_residue()-2):
				#Compute the span for a 4mer
				testresndx=[testres,testres+3]
				#Calculate the RMSD
				rmsd_result=self.rmsd_calculator.get_broadcast_coordinate_rmsd(
							np.expand_dims( pose_res_to_array(result_idealized_pose, testresndx, atoms=["CA","C","O","N"]), 0), 
							cluster_lr_centers_coordinates.copy() 
							).ravel()
				counter=0
				minrmsd=rmsd_result.min()
				minndx=rmsd_result.argmin()
				if(verbose):
					print "Span [%d-%d] RMSD = %f, Limit= %f, Pass = %d, ndx= %d" % (testresndx[0], testresndx[1], minrmsd, 
											self.cluster_lr_centers_store_data["threshold_distance"][minndx],
											(minrmsd <= (self.cluster_lr_centers_store_data["threshold_distance"][minndx])*treshold_factor),
											minndx)
				#Allow max=1% of variation in the max distance 
				if(minrmsd > (self.cluster_lr_centers_store_data["threshold_distance"][minndx]*treshold_factor)):
					has_bad_fragments=True
					for ires in range(testres, testres+4):
						pymol_line+="%d+"%(ires)
				for ndx_shift in range(0, 4):
					aaProbabilies[(testres-1)+ndx_shift].append(self.cluster_lr_centers_store_data["aa"][minndx][ndx_shift])
					#print (testres-1)+ndx_shift
		if(has_bad_fragments):
			print pymol_line
			
		print "Done"
		return aaProbabilies
		
		
	#Get most probable clusters assignments
	def assign_pose_to_clusters_within_cutoff(self, 
												result_idealized_pose, 
												init_res, 
												end_res, 
												cut_off=0.50, 
												clustersDB=None,
												verbose=False, 
												permissive=True):
		
		if (verbose):
			print "Assigning pose with a cut_off of: ", cut_off, ". Be patient... please"
		##Debugline
		result_idealized_pose=result_idealized_pose.clone()
		##End debugline
		
		assignments=[]
		
		aaProbabilies=[]
		for testres in range(0, result_idealized_pose.total_residue()):
			aaProbabilies.append([])
		
		has_bad_fragments=False
		for testres in range( max(1,init_res), min(result_idealized_pose.total_residue()-2, end_res) ):
				#Compute the span for a 4mer
				testresndx=[testres,testres+3]
				#Calculate the RMSD
				rmsd_result=clustersDB.rmsd_calculator.get_broadcast_coordinate_rmsd(
							np.expand_dims( pose_res_to_array(result_idealized_pose, testresndx, atoms=["CA","C","O","N"]), 0), 
							clustersDB.cluster_lr_centers_coordinates.copy() 
							).ravel()
				counter=0
				
				this_tmp_assignments=np.where(rmsd_result<=cut_off)[0]
				if (len(this_tmp_assignments) > 0):
					assignments.append(this_tmp_assignments)
				else:
					if (permissive):
						assignments.append([])
					else:
						print "ERROR: not assignments within the cutoff for span: ", testresndx, "Die"
						assert(0==1)
						
		if (verbose):
			print "Done"
		return assignments
		
		
	def pose_fragment_population_within(self, 
										test_pose_in, 
										init_res, end_res, 
										cut_off=0.51, 
										clustersDB=None):
		tmp_assignments = self.assign_pose_to_clusters_within_cutoff(test_pose_in, 
																	init_res, 
																	end_res, 
																	cut_off, 
																	clustersDB=clustersDB)
		tmp_data_density = np.zeros((len(tmp_assignments)), int)
		for i in xrange(len(tmp_assignments)):
			tmp_data_density[i]=clustersDB.clusters_lr_population[tmp_assignments[i]].sum()
		return tmp_data_density
		
		
	#Function used to query the  Polar Profile of sequences (see previous cell)
	def get_aa_seq_polar_profile(self, seq):
		polar_profile=""
		for i in seq:
			if (i in POLARaa):
				polar_profile+="P"
			elif(i in HYDROPHIBICaa):
				polar_profile+="H"
			elif(i in NEUTRALaa):
				polar_profile+="N"
			else:
				print "ERROR, I dont know the amminoacid: ", i
				assert(0==1)
		return polar_profile
		
		
	def score_region_conections(self, 
								in_pose, 
								res_design_to_start, 
								res_design_to_end,
								des_assignment_cut_off=1.0,
								clustersDB=None,
								b_is_permissive=False,
								b_debug=False):
		if b_debug:
			print "Start,End,Assign-cutoff: ", res_design_to_start, res_design_to_end, des_assignment_cut_off
			#Sanity Check
			assert(res_design_to_start<res_design_to_end)
			
		#Fix: Hacks to check size, make auto consistent!
		res_design_to_start=max(res_design_to_start,1)
		res_design_to_end=min(res_design_to_end, in_pose.split_by_chain()[1].total_residue()-3)
		
		tmp_assignments_design = self.assign_pose_to_clusters_within_cutoff(in_pose, 
																			res_design_to_start, 
																			res_design_to_end,
																			des_assignment_cut_off,
																			clustersDB=clustersDB,
																			permissive=b_is_permissive)
		tmp_sums=[]
		for assigns in tmp_assignments_design:
			##for assign in tmp_assignments_design:
			if (len(assigns) > 0):
				tmp_sums.append(clustersDB.clusters_lr_population[assigns].sum())
				#tmp_mins.append(np.mean(clustersDB.clusters_lr_population[assigns]))
			else:
				return  0.0
		tmp_sums=np.asarray(tmp_sums)
		tmp_avg=np.average(tmp_sums)
		return tmp_sums.min(), tmp_sums.max(), tmp_avg
		
		
	def check_span_for_ideal_fragments(self, 
									test_pose, 
									start, 
									end, 
									rmsd_pert=1.00, 
									max_rmsd_limit=0.50,
									clustersDB=None):
		for testres in range(start , end-3):
			#Compute the span for a 4mer
			testresndx=[testres,testres+3]
			#Calculate the RMSD
			rmsd_result=clustersDB.rmsd_calculator.get_broadcast_coordinate_rmsd(
						np.expand_dims( pose_res_to_array(test_pose, testresndx, atoms=["CA","C","O","N"]), 0), 
						clustersDB.cluster_lr_centers_coordinates.copy() 
						).ravel()
			
			minndx=rmsd_result.argmin()
			if( (rmsd_result[minndx] > max_rmsd_limit) or 
				(rmsd_result[minndx] > (clustersDB.cluster_lr_centers_store_data["threshold_distance"][minndx]*rmsd_pert))):
				return False
		return True
		
		
	def rmsd_2_np_arrays_wRosettaRotationMatrix(self, crds1, crds2):
		rmsd, rMtx, tVec = rmsd_2_np_arrays_wRotationMatrix(crds1, crds2)
		###Are you kidding me??? Is this the correct way to build arrays inside 
		#of Rosetta core? No indexing?? Who was the creator of this???
		rMtx_xyzM=self.rosetta.numeric.xyzMatrix_double()
		rMtx_xyzM.xx=rMtx[0,0]
		rMtx_xyzM.xy=rMtx[0,1]
		rMtx_xyzM.xz=rMtx[0,2]
		rMtx_xyzM.yx=rMtx[1,0]
		rMtx_xyzM.yy=rMtx[1,1]
		rMtx_xyzM.yz=rMtx[1,2]
		rMtx_xyzM.zx=rMtx[2,0]
		rMtx_xyzM.zy=rMtx[2,1]
		rMtx_xyzM.zz=rMtx[2,2]
		tVec_xyzV=self.rosetta.numeric.xyzVector_double()
		tVec_xyzV.x=tVec[0]
		tVec_xyzV.y=tVec[1]
		tVec_xyzV.z=tVec[2]
		return rmsd, rMtx_xyzM, tVec_xyzV
		
		
	def align_bb_by_ndxs(self, pose1, init_res1, end_res1, pose2, init_res2, end_res2):
		numRes=(end_res1-init_res1+1)
		coorA=np.zeros(((4*numRes),3), float)
		coorB=np.zeros(((4*numRes),3), float)
		
		counter=0
		for res in range (init_res1, (end_res1+1)):
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("CA")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("C")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("O")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("N")[dim])
			counter+=1
		counter=0
		for res in range (init_res2, (end_res2+1)):
			for dim in range(0,3): 
				coorB[counter,dim]=(pose2.residue(res).xyz("CA")[dim])
			counter+=1
			for dim in range(0,3): 
				coorB[counter,dim]=(pose2.residue(res).xyz("C")[dim])
			counter+=1
			for dim in range(0,3): 
				coorB[counter,dim]=(pose2.residue(res).xyz("O")[dim])
			counter+=1
			for dim in range(0,3): 
				coorB[counter,dim]=(pose2.residue(res).xyz("N")[dim])
			counter+=1
			
		#Calculate the RMSD
		rmsdVal, rMtx, tVec = self.rmsd_2_np_arrays_wRosettaRotationMatrix(coorB, coorA)
		pose1.apply_transform_Rx_plus_v(rMtx, tVec)
		
		return rmsdVal
		
		
	def align_bb_1pose_to_2poses_by_ndxs_tips(self, pose1, ndxs1a, ndxs1b, pose2, ndxs2, pose3, ndxs3):
		#CHECK Number of atoms in other RMSD calculations!!! NOW (-1 is for the removed atoms)
		numAtoms=len(ndxs1a)+len(ndxs1b)
		numAtoms=(3*numAtoms)
		coorA=np.zeros((numAtoms,3), float)
		coorB=np.zeros((numAtoms,3), float)
		counter=0
		for res in ndxs1a:
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("CA")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("C")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("N")[dim])
			counter+=1
		for res in ndxs1b:
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("CA")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("C")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("O")[dim])
			counter+=1
		counter=0
		for res in ndxs2:
			for dim in range(0,3): 
				coorB[counter,dim]=(pose2.residue(res).xyz("CA")[dim])
			counter+=1
			for dim in range(0,3): 
				coorB[counter,dim]=(pose2.residue(res).xyz("C")[dim])
			counter+=1
			for dim in range(0,3): 
				coorB[counter,dim]=(pose2.residue(res).xyz("N")[dim])
			counter+=1
		for res in ndxs3:
			for dim in range(0,3): 
				coorB[counter,dim]=(pose3.residue(res).xyz("CA")[dim])
			counter+=1
			for dim in range(0,3): 
				coorB[counter,dim]=(pose3.residue(res).xyz("C")[dim])
			counter+=1
			for dim in range(0,3): 
				coorB[counter,dim]=(pose3.residue(res).xyz("O")[dim])
			counter+=1
		#Apply rotation and translation
		rmsdVal, rMtx, tVec = self.rmsd_2_np_arrays_wRosettaRotationMatrix(coorB, coorA)
		pose1.apply_transform_Rx_plus_v(rMtx, tVec)
		#Also calculate the atom pair's max dist
		#Update matrixes
		counter=0
		for res in ndxs1a:
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("CA")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("C")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("N")[dim])
			counter+=1
		for res in ndxs1b:
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("CA")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("C")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("O")[dim])
			counter+=1
		counter=0
		for res in ndxs2:
			for dim in range(0,3): 
				coorB[counter,dim]=(pose2.residue(res).xyz("CA")[dim])
			counter+=1
			for dim in range(0,3): 
				coorB[counter,dim]=(pose2.residue(res).xyz("C")[dim])
			counter+=1
			for dim in range(0,3): 
				coorB[counter,dim]=(pose2.residue(res).xyz("N")[dim])
			counter+=1
		for res in ndxs3:
			for dim in range(0,3): 
				coorB[counter,dim]=(pose3.residue(res).xyz("CA")[dim])
			counter+=1
			for dim in range(0,3): 
				coorB[counter,dim]=(pose3.residue(res).xyz("C")[dim])
			counter+=1
			for dim in range(0,3): 
				coorB[counter,dim]=(pose3.residue(res).xyz("O")[dim])
			counter+=1
		max_pair_dist=0.000
		for i in range(numAtoms):
			tmp_dist = np.linalg.norm(coorA[i]-coorB[i])
			if(tmp_dist > max_pair_dist):
				max_pair_dist = tmp_dist
				
		return rmsdVal, max_pair_dist
		
		
	def align_by_bb(self, pose1, pose2):
		numRes=pose1.total_residue()
		#D assert (numRes == pose2.total_residue())
		coorA=np.zeros(((4*numRes),3), float)
		coorB=np.zeros(((4*numRes),3), float)
		counter=0
		for res in range (1, (numRes+1)):
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("CA")[dim])
				coorB[counter,dim]=(pose2.residue(res).xyz("CA")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("C")[dim])
				coorB[counter,dim]=(pose2.residue(res).xyz("C")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("O")[dim])
				coorB[counter,dim]=(pose2.residue(res).xyz("O")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("N")[dim])
				coorB[counter,dim]=(pose2.residue(res).xyz("N")[dim])
			counter+=1
		rmsdVal, rMtx, tVec = self.rmsd_2_np_arrays_wRosettaRotationMatrix(coorB, coorA)
		pose1.apply_transform_Rx_plus_v(rMtx, tVec)
		return rmsdVal
		
		
	def align_bb_by_fragment(self, pose1, pose2, init_res, end_res):
		numRes=(end_res-init_res+1)
		coorA=np.zeros(((4*numRes),3), float)
		coorB=np.zeros(((4*numRes),3), float)
		counter=0
		for res in range (init_res, (end_res+1)):
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("CA")[dim])
				coorB[counter,dim]=(pose2.residue(res).xyz("CA")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("C")[dim])
				coorB[counter,dim]=(pose2.residue(res).xyz("C")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("O")[dim])
				coorB[counter,dim]=(pose2.residue(res).xyz("O")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("N")[dim])
				coorB[counter,dim]=(pose2.residue(res).xyz("N")[dim])
			counter+=1
		rmsdVal, rMtx, tVec = self.rmsd_2_np_arrays_wRosettaRotationMatrix(coorB, coorA)
		pose1.apply_transform_Rx_plus_v(rMtx, tVec)
		return rmsdVal
		
		
	def align_bb_1pose_to_2poses_by_ndxs(self, pose1, ndxs1a, ndxs1b, pose2, ndxs2, pose3, ndxs3):
		numAtoms=len(ndxs1a)+len(ndxs1b)
		coorA=np.zeros(((4*numAtoms),3), float)
		coorB=np.zeros(((4*numAtoms),3), float)
		
		counter=0
		for res in ndxs1a:
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("CA")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("C")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("O")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("N")[dim])
			counter+=1
		for res in ndxs1b:
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("CA")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("C")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("O")[dim])
			counter+=1
			for dim in range(0,3): 
				coorA[counter,dim]=(pose1.residue(res).xyz("N")[dim])
			counter+=1
			
		counter=0
		for res in ndxs2:
			for dim in range(0,3): 
				coorB[counter,dim]=(pose2.residue(res).xyz("CA")[dim])
			counter+=1
			for dim in range(0,3): 
				coorB[counter,dim]=(pose2.residue(res).xyz("C")[dim])
			counter+=1
			for dim in range(0,3): 
				coorB[counter,dim]=(pose2.residue(res).xyz("O")[dim])
			counter+=1
			for dim in range(0,3): 
				coorB[counter,dim]=(pose2.residue(res).xyz("N")[dim])
			counter+=1
		for res in ndxs3:
			for dim in range(0,3): 
				coorB[counter,dim]=(pose3.residue(res).xyz("CA")[dim])
			counter+=1
			for dim in range(0,3): 
				coorB[counter,dim]=(pose3.residue(res).xyz("C")[dim])
			counter+=1
			for dim in range(0,3): 
				coorB[counter,dim]=(pose3.residue(res).xyz("O")[dim])
			counter+=1
			for dim in range(0,3): 
				coorB[counter,dim]=(pose3.residue(res).xyz("N")[dim])
			counter+=1
			
		rmsdVal, rMtx, tVec = self.rmsd_2_np_arrays_wRosettaRotationMatrix(coorB, coorA)
		pose1.apply_transform_Rx_plus_v(rMtx, tVec)
		return rmsdVal
		
		
	def check_bb_clashes_hash_vs_pose(self, context_pose_hash, pose):
		target_balls = self.rosetta.utility.vector1_Ball()
		self.rosetta.core.pose.xyzStripeHashPose.extract_pose_balls(pose, target_balls, self.rosetta.core.pose.PoseCoordPickMode.N_CA_C_CB)
		target_clash_map = self.rosetta.utility.map_Size_Size()
		return context_pose_hash.clash_check_residue_pairs(target_balls, target_clash_map)
		##for clash_residue_number in target_clash_map.keys():
		##    print clash_residue_number
		
		
	##def: connects 3 poses as pose1--pose2--pose3. It removes the last residue from pose1, 
	## first and last from pose2, and the first from pose3 
	def connect_3_poses_byNDX_and_minize_w_constraints(self,
														pose1, 
														pose2, 
														pose3, 
														ndx_p1b, 
														ndx_p3a,  
														context_pose,
														sfx=None,
														harmonic_constraint_streght=0.5, 
														cart_constraint_weight=0.5):
		
		#ToDo, check this weight
		sfx.set_weight(self.rosetta.core.scoring.cart_bonded , 0.7)
		sfx.set_weight(self.rosetta.core.scoring.coordinate_constraint , cart_constraint_weight)
		
		##Create a tmp pose that will hold the result
		#target_pose_size = p_ss0.total_residue() +  p_ss1.total_residue() + p_loop0_1.total_residue()-4
		##No need of this beacuse we use from 2nd residue
		#if pose2.residue(1).has_variant_type("LOWER_TERMINUS"):
		#    rosetta.core.pose.remove_variant_type_from_pose_residue(pose2, "LOWER_TERMINUS", 1)
		#if pose2.residue(p_loop0_1.total_residue()).has_variant_type("UPPER_TERMINUS"):
		#    rosetta.core.pose.remove_variant_type_from_pose_residue(pose2, "UPPER_TERMINUS", pose2.total_residue())
		tmp_result_pose = self.rosetta.pose_from_sequence(self.fake_aminoacid_letter_for_design)
		
		#This is the first chunk
		tmp_result_pose.copy_segment(1, pose1, 1, 1)
		for i in range( 1, ndx_p1b-1 ): #range( pose1.total_residue()-1):
			tmp_resi=pose1.residue(i+1)
			tmp_result_pose.append_residue_by_bond( tmp_resi, False);
		#This is the loop/connection in the middle, copy all but not EP
		for i in range( 1, pose2.total_residue()-1 ):
			tmp_resi=pose2.residue(i+1)
			tmp_result_pose.append_residue_by_bond( tmp_resi, False);
		#This is the last chunk
		for i in range( ndx_p3a, pose3.total_residue() ): #range( 1, pose3.total_residue() 
			tmp_resi=pose3.residue(i+1)
			tmp_result_pose.append_residue_by_bond( tmp_resi, False);
		
		#make a copy to store the final result at the ned
		pose_for_return_result=tmp_result_pose.clone()
		ori_pose_size=tmp_result_pose.total_residue()
		
		if (context_pose.total_residue() > 0):
			tmp_result_pose.append_pose_by_jump(context_pose, 1)
		
		#Loop indexes
		loop_start=ndx_p1b #pose1.total_residue()
		loop_end=ndx_p1b+pose2.total_residue()-2   #pose1.total_residue() + pose2.total_residue()-2
		
		#Constraints
		mm_local = self.rosetta.core.kinematics.MoveMap()
		for i in range(loop_start-1, loop_end+1):
			mm_local.set_bb(i, True)
			#Harmonic constraints
			for j in range( 4 ):
				tmp_result_pose.add_constraint( self.CoorCstr(self.rosetta.core.id.AtomID(j+1, i), 
						self.rosetta.core.id.AtomID(1,1), 
						tmp_result_pose.residue(i).atom(j+1).xyz(),
						self.HarmFunc(0.0, harmonic_constraint_streght)) ) 
		
		options_minilbfgs = self.rosetta.core.optimization.MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, True, False, False)
		
		minimizer = self.CartMin()
		score = minimizer.run( tmp_result_pose, mm_local, sfx, options_minilbfgs )
		
		#Remove constraints from the sfx
		sfx.set_weight(self.rosetta.core.scoring.cart_bonded , 0.0)
		sfx.set_weight(self.rosetta.core.scoring.coordinate_constraint , 0.0)
		#Remove constraints form the pose:
		##tmp_result_pose.remove_constraints()
		
		#finaly copy back just the fragment without context
		pose_for_return_result.copy_segment(ori_pose_size, tmp_result_pose, 1, 1)
		return pose_for_return_result, score, score #scorefxn_remodel_cen(pose_for_mini)
		
		
	def minimize_sidechains(self, 
							pose_for_mini, 
							sfx=None):
		#Fist clear the constraint set
		pose_for_mini.remove_constraints()
		
		#Constraints
		mm_local = self.rosetta.core.kinematics.MoveMap()
		mm_local.set_bb(False)
		mm_local.set_chi(True)
		
		sfx.set_weight(self.rosetta.core.scoring.coordinate_constraint , 1.0)
		
		options_minilbfgs = self.rosetta.core.optimization.MinimizerOptions( "lbfgs_armijo_nonmonotone", 1e-5, True, False, False)
		
		minimizer = self.mzr()
		score = minimizer.run( pose_for_mini, mm_local, sfx, options_minilbfgs )
		sfx.set_weight(self.rosetta.core.scoring.coordinate_constraint , 0.0)
		pose_for_mini.remove_constraints()
		
		return score
		
		
	def minimize_all_with_bb_constraints(self,
										pose_for_mini, 
										sfx=None,
										harmonic_constraint_streght=0.5, 
										cart_constraint_weight=1.0):
		#Fist clear the constraint set
		pose_for_mini.remove_constraints()
		
		#Now add new constraints  
		mm_local = self.rosetta.core.kinematics.MoveMap()
		mm_local.set_chi(True)
		
		#Rosetta... Coordinate contraints are a mess, not what they should be... let's think
		###pose_for_mini_reference=pose_for_mini.clone()
		for i in range(1, pose_for_mini.total_residue()+1):
			mm_local.set_bb(i, True)
			#Harmonic constraints
			for j in range( 4 ):
				pose_for_mini.add_constraint( self.CoorCstr(self.rosetta.core.id.AtomID(j+1, i), 
						self.rosetta.core.id.AtomID(1, 1), 
						pose_for_mini.residue(i).atom(j+1).xyz(),
						self.HarmFunc(0.0, harmonic_constraint_streght)) ) 
						
		sfx.set_weight(self.rosetta.core.scoring.coordinate_constraint , cart_constraint_weight)
		#Minimizer options
		options_minilbfgs = self.rosetta.core.optimization.MinimizerOptions( "lbfgs_armijo_nonmonotone", 1e-5, True, False, False)
		#Instantiate the minimizer
		minimizer = self.mzr()
		#Run!
		score = minimizer.run( pose_for_mini, mm_local, sfx, options_minilbfgs )
		sfx.set_weight(self.rosetta.core.scoring.coordinate_constraint , 0.0)
		pose_for_mini.remove_constraints()
		
		return score
		
		
	def relax_with_constraints_to_init(self,
										pose_for_relax, 
										sfx=None):
		#Fist clear the constraint set
		pose_for_relax.remove_constraints()
		sfx.set_weight(self.rosetta.core.scoring.coordinate_constraint , 1.0)
		fastRel=self.FastRel(sfx)
		fastRel.constrain_relax_to_start_coords(True)
		fastRel.apply(pose_for_relax)
		sfx.set_weight(self.rosetta.core.scoring.coordinate_constraint , 0.0)
		pose_for_relax.remove_constraints()
		
		
	def relax_without_constraints_to_init(self,
											pose_for_relax, 
											sfx=None):
		#Fist clear the constraint set
		pose_for_relax.remove_constraints()
		sfx.set_weight(self.rosetta.core.scoring.coordinate_constraint , 0.0)
		
		fastRel=self.FastRel(sfx)
		fastRel.constrain_relax_to_start_coords(False)
		fastRel.apply(pose_for_relax)
		
		sfx.set_weight(self.rosetta.core.scoring.coordinate_constraint , 0.0)
		pose_for_relax.remove_constraints()
		
		
	##def: connects 3 poses as pose1--pose2--pose3. It removes the last residue from pose1, 
	## first and last from pose2, and the first from pose3 
	def connect_3_poses_byNDX_and_add_context(self,
												pose1, 
												pose2, 
												pose3, 
												ndx_p1b, 
												ndx_p3a, 
												context_pose, 
												in_scorefxn,
												harmonic_constraint_streght=0.5, 
												cart_constraint_weight=1.0):
		
		##Create a tmp pose that will hold the result
		#target_pose_size = p_ss0.total_residue() +  p_ss1.total_residue() + p_loop0_1.total_residue()-4
		##No need of this beacuse we use from 2nd residue
		#if pose2.residue(1).has_variant_type("LOWER_TERMINUS"):
		#    rosetta.core.pose.remove_variant_type_from_pose_residue(pose2, "LOWER_TERMINUS", 1)
		#if pose2.residue(p_loop0_1.total_residue()).has_variant_type("UPPER_TERMINUS"):
		#    rosetta.core.pose.remove_variant_type_from_pose_residue(pose2, "UPPER_TERMINUS", pose2.total_residue())
		tmp_result_pose = self.rosetta.pose_from_sequence(self.fake_aminoacid_letter_for_design)
		tmp_result_pose.copy_segment(1, pose1, 1, 1)
		for i in range( 1, ndx_p1b-2): #pose1.total_residue()-2):
			tmp_resi=pose1.residue(i+1)
			tmp_result_pose.append_residue_by_bond( tmp_resi, False);
		for i in range( 0, pose2.total_residue() ):
			tmp_resi=pose2.residue(i+1)
			tmp_result_pose.append_residue_by_bond( tmp_resi, False);
		for i in range( ndx_p3a+1, pose3.total_residue() ):  #2, pose3.total_residue() ):
			tmp_resi=pose3.residue(i+1)
			tmp_result_pose.append_residue_by_bond( tmp_resi, False);
		
		
		#make a copy
		tmp_result_pose_w_context=tmp_result_pose.clone()
		
		if (context_pose.total_residue() > 0):
			tmp_result_pose_w_context.append_pose_by_jump(context_pose, 1)
		
		##Loop indexes
		##loop_start=pose1.total_residue()
		##loop_end=pose1.total_residue()+ pose2.total_residue()-2
		
		###Constraints
		##mm_local = rosetta.core.kinematics.MoveMap()
		##for i in range(loop_start-1, loop_end+1):
		##    mm_local.set_bb(i, True)
		##    #Harmonic constraints
		##    for j in range( 4 ):
		##        pose_for_mini.add_constraint( self.CoorCstr(rosetta.core.id.AtomID(j+1, i), 
		##                rosetta.core.id.AtomID(2, pose_for_mini.total_residue()), tmp_result_pose.residue(i).atom(j+1).xyz(),
		##                self.HarmFunc(0.0, 1.0)) )   
		
		##options_minilbfgs = rosetta.core.optimization.MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, True, False, False)
		
		###minimizer = self.CartMin()
		"""
		score = minimize_all_with_bb_constraints(tmp_result_pose_w_context.clone(), 
												 sfx=in_scorefxn, 
												 harmonic_constraint_streght=harmonic_constraint_streght, 
												 cart_constraint_weight=cart_constraint_weight)     
		"""
		
		score=in_scorefxn(tmp_result_pose_w_context)
		
		#finaly copy back
		##tmp_result_pose.copy_segment(tmp_result_pose.total_residue(), pose_for_mini, 1, 1)
		return tmp_result_pose, tmp_result_pose_w_context, tmp_result_pose_w_context, score
		
		
	##def: connects 3 poses as pose1--pose2--pose3. It removes the last residue from pose1, 
	## first and last from pose2, and the first from pose3 
	def connect_3_poses_byEP_and_add_context(self,
												pose1, 
												pose2, 
												pose3, 
												context_pose, 
												in_scorefxn,
												harmonic_constraint_streght=0.5, 
												cart_constraint_weight=1.0):
		
		##Create a tmp pose that will hold the result
		#target_pose_size = p_ss0.total_residue() +  p_ss1.total_residue() + p_loop0_1.total_residue()-4
		##No need of this beacuse we use from 2nd residue
		#if pose2.residue(1).has_variant_type("LOWER_TERMINUS"):
		#    rosetta.core.pose.remove_variant_type_from_pose_residue(pose2, "LOWER_TERMINUS", 1)
		#if pose2.residue(p_loop0_1.total_residue()).has_variant_type("UPPER_TERMINUS"):
		#    rosetta.core.pose.remove_variant_type_from_pose_residue(pose2, "UPPER_TERMINUS", pose2.total_residue())
		tmp_result_pose = self.rosetta.pose_from_sequence(self.fake_aminoacid_letter_for_design)
		tmp_result_pose.copy_segment(1, pose1, 1, 1)
		#for i in range( 1, pose1.total_residue()-1):
		for i in range( 1, pose1.total_residue()-2):
			tmp_resi=pose1.residue(i+1)
			tmp_result_pose.append_residue_by_bond( tmp_resi, False);
		#for i in range( 1, pose2.total_residue()-1 ):
		for i in range( 0, pose2.total_residue() ):
			tmp_resi=pose2.residue(i+1)
			tmp_result_pose.append_residue_by_bond( tmp_resi, False);
		#for i in range( 1, pose3.total_residue() ):
		for i in range( 2, pose3.total_residue() ):
			tmp_resi=pose3.residue(i+1)
			tmp_result_pose.append_residue_by_bond( tmp_resi, False);
		
		#make a copy
		pose_for_mini=tmp_result_pose.clone()
		
		if (context_pose.total_residue() > 0):
			pose_for_mini.append_pose_by_jump(context_pose, 1)
		
		pose_before_mini=pose_for_mini.clone()
		score = minimize_all_with_bb_constraints(pose_for_mini, 
												 sfx=in_scorefxn, 
												 harmonic_constraint_streght=harmonic_constraint_streght, 
												 cart_constraint_weight=cart_constraint_weight) 
		
		score = minimize_all_with_bb_constraints(pose_for_mini, sfx=in_scorefxn) 
		
		#finaly copy back
		##tmp_result_pose.copy_segment(tmp_result_pose.total_residue(), pose_for_mini, 1, 1)
		return tmp_result_pose, tmp_result_pose, pose_for_mini, score
		
		
	##def: connects 3 poses as pose1--pose2--pose3. It removes the last residue from pose1, 
	## first and last from pose2, and the first from pose3 
	def connect_3_poses_byEP_and_minize_w_constraints(pose1, 
													  pose2, 
													  pose3, 
													  context_pose, 
													  sfx=None,
													  harmonic_constraint_streght=0.5, 
													  cart_constraint_weight=0.5):
		
		#ToDo, check this weight
		sfx.set_weight(self.rosetta.core.scoring.cart_bonded , 0.7)
		sfx.set_weight(self.rosetta.core.scoring.coordinate_constraint , cart_constraint_weight)
		
		
		##Create a tmp pose that will hold the result
		#target_pose_size = p_ss0.total_residue() +  p_ss1.total_residue() + p_loop0_1.total_residue()-4
		##No need of this beacuse we use from 2nd residue
		#if pose2.residue(1).has_variant_type("LOWER_TERMINUS"):
		#    rosetta.core.pose.remove_variant_type_from_pose_residue(pose2, "LOWER_TERMINUS", 1)
		#if pose2.residue(p_loop0_1.total_residue()).has_variant_type("UPPER_TERMINUS"):
		#    rosetta.core.pose.remove_variant_type_from_pose_residue(pose2, "UPPER_TERMINUS", pose2.total_residue())
		tmp_result_pose = self.rosetta.pose_from_sequence(fake_aminoacid_letter_for_design)
		tmp_result_pose.copy_segment(1, pose1, 1, 1)
		for i in range( 1, pose1.total_residue()-1):
			tmp_resi=pose1.residue(i+1)
			tmp_result_pose.append_residue_by_bond( tmp_resi, False);
		for i in range( 1, pose2.total_residue()-1 ):
			tmp_resi=pose2.residue(i+1)
			tmp_result_pose.append_residue_by_bond( tmp_resi, False);
		for i in range( 1, pose3.total_residue() ):
			tmp_resi=pose3.residue(i+1)
			tmp_result_pose.append_residue_by_bond( tmp_resi, False);
		
		#make a copy
		pose_for_return_result=tmp_result_pose.clone()
		ori_pose_size=tmp_result_pose.total_residue()
		
		if (context_pose.total_residue() > 0):
			tmp_result_pose.append_pose_by_jump(context_pose, 1)
		
		#Loop indexes
		loop_start=pose1.total_residue()
		loop_end=pose1.total_residue()+ pose2.total_residue()-2
		
		#Constraints
		mm_local = self.rosetta.core.kinematics.MoveMap()
		for i in range(loop_start-1, loop_end+1):
			mm_local.set_bb(i, True)
			#Harmonic constraints
			for j in range( 4 ):
				tmp_result_pose.add_constraint( self.CoorCstr(self.rosetta.core.id.AtomID(j+1, i), 
						self.rosetta.core.id.AtomID(1,1), 
						tmp_result_pose.residue(i).atom(j+1).xyz(),
						self.HarmFunc(0.0, harmonic_constraint_streght)) )  
		
		options_minilbfgs = self.rosetta.core.optimization.MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, True, False, False)
		
		minimizer = self.CartMin()
		score = minimizer.run( tmp_result_pose, mm_local, sfx, options_minilbfgs )
		
		#Remove constraints from the sfx
		sfx.set_weight(self.rosetta.core.scoring.cart_bonded , 0.0)
		sfx.set_weight(self.rosetta.core.scoring.coordinate_constraint , 0.0)
		#Remove constraints form the pose:
		tmp_result_pose.remove_constraints()
		
		#finaly copy back
		#finaly copy back
		pose_for_return_result.copy_segment(ori_pose_size, tmp_result_pose, 1, 1)
		##tmp_result_pose.copy_segment(tmp_result_pose.total_residue(), pose_for_mini, 1, 1)
		return pose_for_return_result, score, score #scorefxn_remodel_cen(tmp_result_pose)
		
	
	def combine_SS_and_loops(self,
							idealized_pose_array=None,
							loop_poses_array=None,
							entropic_acceptance_ratio_score=0.70,
							general_out_dir="./",
							out_prefix_name="stage1_result"):
		#Connect SSs and its best Loops
		###Parameters
		##min_num_seq_to_match=5
		##entropic_acceptance_ratio_designability=0.95
		##small_fragment_bias_factor=1.0
		####
		pose_buffer1 = idealized_pose_array[0].clone()
		pose_buffer2 = idealized_pose_array[0].clone()
		best_cen_score=np.finfo('d').max
		best_score_hbonds=0.0
		print "Total SS to connect:", len(loop_poses_array)
		
		#Relative indexes of in and out pose
		tmp_counter=0
		original_ss_relative_indexes=[]
		for indx in range(1,idealized_pose_array[0].total_residue()+1):
			original_ss_relative_indexes.append([0,tmp_counter,(indx-1)])
			tmp_counter+=1
		#Loops ss-loops to score and choose
		for i in range(len(loop_poses_array)):
			print "Connecting SS:", i, "<-->", i+1
			#Generate a context pose for optimization of non-local interactions:
			b_is_first=True
			tmp_context_pose = self.rosetta.pose_from_sequence("")
			for j in range(i+2, len(idealized_pose_array)):
				if(b_is_first):
					tmp_context_pose = idealized_pose_array[j].clone()
					b_is_first=False
				else:
					tmp_context_pose.append_pose_by_jump(idealized_pose_array[j],1)
			
			#Sort the loops based on it's local interaction energy. 
			#The structure of the container is: [RMS, pose, pose, score, popscore, delta1, delta2] sort by most populated first = sort(-popscore) 
			#Extract the loops to test
			##loops_to_test = loop_poses_array[i] ###sorted(loop_poses_array[i], key=lambda popscore: -popscore[4])
			loops_to_test = np.copy(loop_poses_array[i]) # sorted(loop_poses_array[i], key=lambda popscore: -popscore[4])
			
			
			"""
			###JUST A TEST REMOVE
			#HACK: score loops again this will be moved to the loop finding routine
			for j in range(len(loops_to_test)):
				[tmp_score_min, tmp_score_max, tmp_score_avg]=score_region_conections(loops_to_test[j][1], 
											1, 
											loops_to_test[j][1].total_residue(),
											des_assignment_cut_off=1.0,
											b_is_permissive=False)
				#Loop scoring function! doh!
				loops_to_test[j][4]= (tmp_score_avg*tmp_score_avg)/(tmp_score_avg-tmp_score_min)
				#print j, loops_to_test[j][4], tmp_score_min, tmp_score_avg
			
			loops_to_test=sorted(loops_to_test, key=lambda popscore: -popscore[4])
			for j in range(len(loops_to_test)):
				###print j,loops_to_test[j][4]
				loops_to_test[j][1].dump_pdb("%s/test/test_%03d_%03d_%03d.pdb"%(general_out_dir, i,j,loops_to_test[j][4])) 
			##break
			###END REMOVE
			"""
			best_population_score=0
			best_loop_len=99999
			has_solution=False
			best_score_hbonds=0.0
			best_score_cen=np.finfo('d').max
			pose_buffer2=pose_buffer1.clone()
			best_pose_with_context=pose_buffer1.clone()
			#Check loops from larger to smaller :)
			for j in reversed(range(len(loops_to_test))): # min(len(loops_to_test), 100)):
				if ( (j%50) == 0 ):
					print "Be patient, until  now I have analyzed: ", j, "of ", len(loops_to_test), "possible loops" 
					
				#####
				#Check clashes    
				#Create a hash of what is there
				#First remove the end (and consider the allowed A delta)
				tmp_pose_size=pose_buffer2.total_residue()-loops_to_test[j][5]-3
				assert(tmp_pose_size > 0)
				tmpline=""
				for kres in range(0, tmp_pose_size):
					tmpline+=self.fake_aminoacid_letter_for_design
				pose_buffer2_for_clash_check=self.rosetta.pose_from_sequence(tmpline)
				pose_buffer2_for_clash_check.copy_segment(tmp_pose_size,pose_buffer2,1,1)
				built_pose_hash = self.rosetta.core.pose.xyzStripeHashPose(pose_buffer2_for_clash_check,
							self.rosetta.core.pose.PoseCoordPickMode.N_CA_C_CB)
				##pose_buffer2_for_clash_check.dump_pdb("%s/test_buffer_%d_%d.pdb" % (general_out_dir, i,j))
				###loops_to_test[j][1].dump_pdb("%s/test_loop_%d_%d.pdb"%(general_out_dir,i,j))
				hasBBclashes = self.check_bb_clashes_hash_vs_pose(built_pose_hash, loops_to_test[j][1])
			   
				if(hasBBclashes):
					#print "Clashes found in loop:", i,"-",j, ". Skipping this."
					continue
				#End Clash check
				#####
				
				#Join poses and get scores
				tmp_ndx2=pose_buffer2.total_residue()-loops_to_test[j][5] #Current minus Delta
				tmp_ndx3=1+loops_to_test[j][6] #Upcoming plus rosetta's stu9id no. 1
				test_pose, test_pose_w_context, test_pose_w_context_mini, score_cen = self.connect_3_poses_byNDX_and_add_context(
																	pose_buffer2,  
																	loops_to_test[j][1],
																	idealized_pose_array[i+1], 
																	tmp_ndx2,
																	tmp_ndx3,
																	tmp_context_pose,
																	self.scorefxn_remodel_cen)
				
				#Normalize score by number of residues:
				#Note, should we minimize first?
				score_cen = score_cen/test_pose.total_residue()
				###print score_cen
				curr_population_score=loops_to_test[j][4]
				curr_loop_len=loops_to_test[j][1].total_residue()-4
				#Debug
				###print "j ,score_cen, curr_population_score, best_population_score: ", j, score_cen, curr_population_score, best_population_score
				##print score_cen, best_score_cen
				##test_pose.dump_pdb("%s/tmppose_%d_%d.pdb"%(general_out_dir,i,j))
				#End Debug
				
				#All pose Hbonding number
				hbond_set_complete_pose = self.rosetta.core.scoring.hbonds.HBondSet()
				test_pose_w_context_mini.update_residue_neighbors();
				self.rosetta.core.scoring.hbonds.fill_hbond_set( test_pose_w_context, False, hbond_set_complete_pose )
				num_hbonds_complete_pose = hbond_set_complete_pose.nhbonds()
				
				#internal loop habond number
				hbond_set_this_loop = self.rosetta.core.scoring.hbonds.HBondSet()
				loops_to_test[j][1].update_residue_neighbors();
				self.rosetta.core.scoring.hbonds.fill_hbond_set( loops_to_test[j][1], False, hbond_set_this_loop )
				num_hbonds_this_loop = hbond_set_this_loop.nhbonds()
				
				num_hbonds_score=-1
				if ( curr_loop_len < 4 ):
					num_hbonds_score = num_hbonds_complete_pose 
				else:    
					num_hbonds_score = num_hbonds_complete_pose + (num_hbonds_this_loop - curr_loop_len)
					###print (num_hbonds_this_loop - curr_loop_len)
				
				###print num_hbonds_score, curr_loop_len
				
				##print j, num_hbonds_score, best_score_hbonds, best_score_cen
				#Algorithm to decide combinations
				##if (score_cen < (best_score_cen) ): #) ):
				##print num_hbonds_score, best_score_hbonds
				if (num_hbonds_score >= best_score_hbonds ):
					#Accept within certain delta
					##if (score_cen < (best_score_cen*entropic_acceptance_ratio_score) ):
					if( curr_population_score >= (best_population_score*entropic_acceptance_ratio_score) ):
							##print curr_population_score
							if(curr_loop_len < best_loop_len):
								best_score_cen = score_cen
								best_score_hbonds = num_hbonds_score
								best_population_score=curr_population_score
								best_loop_len=curr_loop_len
								pose_buffer1=test_pose.clone()
								best_pose_with_context=test_pose_w_context_mini.clone()
								has_solution=True
								#print ("I have found a better (more designable and smaller) loop, #", j, ", stats (len, PospScore, scoreCen): ", 
								#       best_loop_len, best_population_score, best_score_cen, num_hbonds_score)
								
							elif(curr_loop_len == best_loop_len):
								best_score_cen = score_cen
								best_score_hbonds = num_hbonds_score
								best_population_score=curr_population_score
								best_loop_len=curr_loop_len
								pose_buffer1=test_pose.clone()
								best_pose_with_context=test_pose_w_context_mini.clone()
								has_solution=True
								#print ("I have found a better (more designable) loop, #", j, ", stats (len, PospScore, scoreCen): ", 
								#       best_loop_len, best_population_score, best_score_cen, num_hbonds_score)
							
							elif(curr_loop_len > best_loop_len):
								##if((curr_population_score*entropic_acceptance_ratio_score) > (best_population_score)):
									best_score_cen = score_cen
									best_score_hbonds = num_hbonds_score
									best_population_score=curr_population_score
									best_loop_len=curr_loop_len
									pose_buffer1=test_pose.clone()
									best_pose_with_context=test_pose_w_context_mini.clone()
									has_solution=True
									#print ("I have found a better (more designable but longer) loop, #", j, ", stats (len, PospScore, scoreCen): ", 
									#       best_loop_len, best_population_score, best_score_cen, num_hbonds_score)
							else:
								print "I Shouldn't be here! Stopping"
								assert(0==1)
				"""    
				elif( curr_population_score > (best_population_score)):
					 if( score_cen <= (best_score_cen*entropic_acceptance_ratio_score) ):
							best_score_cen = score_cen
							best_score_hbonds = num_hbonds_score
							best_population_score=curr_population_score
							best_loop_len=curr_loop_len
							pose_buffer1=test_pose.clone()
							best_pose_with_context=test_pose_w_context_mini.clone()
							has_solution=True
							print ("I have found a better  loop, #", j, ", stats (len, PospScore, scoreCen, num_hbonds_score): ", 
								   best_loop_len, best_population_score, best_score_cen)
							#ToDo: Make it output various solutions! Not only the first
							#break
				"""
			#If has solution store and continue
			if(has_solution):
				print "Best solution is (numRes, scorePop, scoreE): ", best_loop_len, best_population_score, best_score_cen, best_score_hbonds ## pose_buffer1.total_residue(), best_score_cen
				connected_pose=pose_buffer1.clone()
				#ToDo, complete original_ss_relative_indexes with information about the original pose (missing righ now)
				tmp_counter=0
				for indx in range( (connected_pose.total_residue()-idealized_pose_array[i+1].total_residue()+1), connected_pose.total_residue()+1):
					original_ss_relative_indexes.append([i+1,tmp_counter,(indx-1)])
					tmp_counter+=1
				##original_ss_relative_indexes.append(range(pose_buffer2.total_residue()+loops_to_test[j][1].total_residue()-3, connected_pose.total_residue()+1))
				###print "Relative original indexes: ", original_ss_relative_indexes
				##connected_pose.dump_pdb("%s/merged_test_phase1_C%d.pdb"%(general_out_dir, i))
				##best_pose_with_context.dump_pdb("%s/merged_test_phase1_C%d_wContext.pdb"%(general_out_dir, i))
			else:
				print "Couldn't find a solution, breaking"
				break
		#Conditional return statment
		if has_solution:
			#Dump PDB result
			tmp_outname="%s/%s_result.pdb"%(general_out_dir, out_prefix_name)
			print "Writing non-minimized result pose to: ", tmp_outname
			connected_pose.dump_pdb(tmp_outname)
			print "DONE connecting loops!"
			return True, connected_pose, original_ss_relative_indexes
		else:
			print "Loop connection failed, there is no non-clashing combinations!"
			return False, connected_pose, original_ss_relative_indexes
		
		
	def minimize_calculate_difficulty_and_plot(self,
											connected_pose=None, 
											b_make_plots=None,
											clustersDB=None,
											general_out_dir="./",
											out_prefix_name="stage1_result"):
											
		#if (b_make_plots):
		#Assign to fragments within cut-off
		population_score_max_rmsd=0.8
		tmp_assignments = self.assign_pose_to_clusters_within_cutoff(connected_pose, 
																1, 
																connected_pose.total_residue(), 
																cut_off=population_score_max_rmsd,
																clustersDB=clustersDB)
		
		tmp_data_density_1=np.zeros((2,len(tmp_assignments)), int)
		
		#Debug
		###print connected_pose.total_residue(), len(tmp_assignments)
		#End Debug
		
		for i in xrange(len(tmp_assignments)):
			tmp_data_density_1[0][i]=i
			tmp_data_density_1[1][i]=clustersDB.clusters_lr_population[tmp_assignments[i]].sum()
		
		#Minimize and idealize bond distances!!!
		##connected_pose_mini, score = minimize_with_constraints_cen(connected_pose)
		#ToDo, change name "harmonic_constraint_streght"-> "harmonic_constraint_sd"
		print "Minimizing Pose, Be patient"
		connected_pose_mini = connected_pose.clone()
		score = self.minimize_all_with_bb_constraints(connected_pose_mini, 
												 sfx=self.scorefxn_remodel_cen, 
												 harmonic_constraint_streght=2.0, 
												 cart_constraint_weight=1.0)
		score = self.minimize_all_with_bb_constraints(connected_pose_mini, 
												 sfx=self.scorefxn_remodel_cen, 
												 harmonic_constraint_streght=2.0, 
												 cart_constraint_weight=0.8)
		score = self.minimize_all_with_bb_constraints(connected_pose_mini, 
												 sfx=self.scorefxn_remodel_cen, 
												 harmonic_constraint_streght=2.0, 
												 cart_constraint_weight=0.6)
		score = self.minimize_all_with_bb_constraints(connected_pose_mini, 
												 sfx=self.scorefxn_remodel_cen, 
												 harmonic_constraint_streght=2.0, 
												 cart_constraint_weight=0.4)
		score = self.minimize_all_with_bb_constraints(connected_pose_mini, 
												 sfx=self.scorefxn_remodel_cen, 
												 harmonic_constraint_streght=2.0, 
												 cart_constraint_weight=0.2)
		tmp_rmsd_after_min=self.align_by_bb(connected_pose_mini, connected_pose)
		print "Minimization Done!"
		print "RMSD after minimization = ", tmp_rmsd_after_min
		
		#Dump PDB result
		tmp_outname="%s/%s_result_minimized.pdb"%(general_out_dir,out_prefix_name)
		print "Writing result pose to: ", tmp_outname
		connected_pose_mini.dump_pdb(tmp_outname)
		
		tmp_assignments = self.assign_pose_to_clusters_within_cutoff(connected_pose_mini, 
																1, 
																connected_pose_mini.total_residue(), 
																cut_off=population_score_max_rmsd,
																clustersDB=clustersDB)
		#D assert len(tmp_assignments)
		tmp_data_density_2=np.zeros((2,len(tmp_assignments)), int)
		for i in xrange(len(tmp_assignments)):
			tmp_data_density_2[0][i]=i
			tmp_data_density_2[1][i]+=clustersDB.clusters_lr_population[tmp_assignments[i]].sum()
		#tmp_data_density_2[1] = tmp_data_density_2[1]/2
		min_dif_score=(1.0/tmp_data_density_2[1]).min()
		max_dif_score=(1.0/tmp_data_density_2[1]).max()
		if (b_make_plots):
			#Plot
			fig = plt.figure(figsize=(15,10))
			ax = fig.add_subplot(2,1,1)
			ax.set_yscale('log')
			#plot data before minimization
			ax.plot(tmp_data_density_1[0], 1.0/tmp_data_density_1[1], linewidth=2.0,  color='black')
			#plot data after minimization
			ax.plot(tmp_data_density_2[0], 1.0/tmp_data_density_2[1], linewidth=2.0,  color='red')
			#Set y limit
			ax.set_ylim([0, .001])
			
			#Save to PNG
			tmp_outname="%s/%s_plot_fdifficulty.png"%(general_out_dir,out_prefix_name)
			print "Writing fancy dificulty graph to: ", tmp_outname
			fig.savefig( tmp_outname )
		return connected_pose_mini, min_dif_score, max_dif_score 
		
		
	def read_input_pdbs_stage2(self, 
								stage2_input_pdb_name=None, 
								general_out_dir="./"):
		tmp_pdb_name="%s/%s"%(general_out_dir,stage2_input_pdb_name, )
		print "Reading input PDB from: ", tmp_pdb_name
		connected_pose_tmp = self.rosetta.pose_from_pdb(tmp_pdb_name)
		print "Done reading"
		
		print "In sequence: ", connected_pose_tmp.sequence()
		print "Sequence Len:", connected_pose_tmp.total_residue()
		
		#Idealize pose bond's distances
		#ToDO add RMSD comparision before and after
		print "Idealizing bond lenghts"
		connected_pose_tmp_backup=connected_pose_tmp.clone()
		idealization_sanity_check_b = self.idealize_bond_distances(connected_pose_tmp)
		if not idealization_sanity_check_b:
			print "Something went wrong during the bond idealization routine, please check your input PDB."
			assert (False)
		tmp_rmsd_after_bond_idealization=self.align_by_bb(connected_pose_tmp, connected_pose_tmp_backup)
		print "Idealizing bond lenghts Done!"
		print "RMSD after bond lenght idealization = ", tmp_rmsd_after_bond_idealization
		
		print "DONE"
		return connected_pose_tmp
		
		
	def setup_pseudo_symmetric_packing(self,
									connected_pose_mini_desig=None,
									b_use_symmetric_packing=False, 
									pseudo_symm_starts=-1, 
									pseudo_symm_ends=-1):
		if b_use_symmetric_packing:
			print "Setting symmetry packing booleans"
			#Defines the shift between the residues
			pseudo_symm_phase_shift=pseudo_symm_ends-pseudo_symm_starts+1
			symm_positions_dic={}
			symm_positions=[]
			print "Enabling symmetric packing. Shifts are:"
			for i in range(pseudo_symm_starts, pseudo_symm_ends+1):
				symm_positions_dic[i]=[]
				for j in range(i+pseudo_symm_phase_shift, connected_pose_mini_desig.total_residue()+1,pseudo_symm_phase_shift):
					if i in symm_positions_dic:
						symm_positions_dic[i].append(j)
					#else:
					#    symm_positions_dic[i]=[j]
						
			tmp_symm_positions_dic=symm_positions_dic.copy()
			for i in tmp_symm_positions_dic:
				for j in tmp_symm_positions_dic[i]:
						symm_positions_dic[j]=[i]
						for k in tmp_symm_positions_dic[i]:
							if (j!=k):
								symm_positions_dic[j].append(k)
						
			for i in range(pseudo_symm_starts, pseudo_symm_ends+1-3):
				symm_positions.append(i)
			###print symm_positions
			print "Symmetric positions list: "
			for i in symm_positions_dic:
				print i,"-->", symm_positions_dic[i]
			return symm_positions_dic
		else:
			print "Using non-symmetric packing"
			
			
	def calculate_pose_layers(self,
							connected_pose_mini_desig=None,
							layersDic=None,
							testSetAA=['V'],  #array, a.a. to test
							c_max_Cb_val=[198.00], #array, Constant describing the maximmum Cb SASA per a.a.
							core_cut_off=0.03,
							boundary_cut_off=0.1,
							general_out_dir="./",
							stage2_output_prefix="stage2_"):
		print "Calculating pose layers"
		print "Fisrt calculating pose layer booleans bvased on SASA"
		cb_vals_array=[]
		for testAA in testSetAA:
			print "Testing with poly-a.a.(s):", testAA, ", and corresponding SASA limit(s):", c_max_Cb_val 
			#Make polyAA pose with no repack
			#print "In sequence: ", connected_pose_mini.sequence(), " Len:", connected_pose_mini_desig.total_residue()
			connected_pose_mini_polyAA=self.make_polyAA_pose(connected_pose_mini_desig, targetAA=testAA)
			#print "PolyAA sequence: ", connected_pose_mini_polyAA.sequence(), " Len:", connected_pose_mini_polyAA.total_residue()
			#connected_pose_mini_polyAA.dump_pdb("%s/merged_test_phase1_final_min_c_polVal.pdb"%general_out_dir)
			#print "Done"
			#Calculate SASA
			total_sasa, per_res_sasa_list, per_resAtom_sasa_list_dic = self.calc_pose_sasa_total_per_residue_and_atom(
																		connected_pose_mini_polyAA, probe_raddi=1.5)
			#Calculate Cb SASA value for all but not N-C terminus
			cb_val=np.zeros(connected_pose_mini_polyAA.total_residue(), float)
			#cb_val[0]=0.5
			#cb_val[len(cb_averages_pose)-1]=0.5
			for i in  range(len (per_res_sasa_list) ):  #connected_pose_mini_polyAA.total_residue() ):
				#print per_resAtom_sasa_list_dic[i]
				#if "CA" in per_resAtom_sasa_list_dic[i]:
					cb_val[i]=per_res_sasa_list[i] #per_resAtom_sasa_list_dic[i]["CA"]       
				#else:
				#    assert(0==1)
			cb_vals_array.append(np.copy(cb_val))
			
		###print cb_vals_array
		print "Max for testset: ", testSetAA, "is:", cb_vals_array[0].max()
		##print cb_vals_array[1].max()
		##print cb_vals_array[2].max()
		
		cb_sum=np.zeros(connected_pose_mini_polyAA.total_residue(), float)
		for i in range(connected_pose_mini_polyAA.total_residue() ):
			for j in range(len(testSetAA)):
				cb_sum[i] += cb_vals_array[j][i]/c_max_Cb_val[j]
		cb_exposition_ratio=cb_sum/3
		lineout=""

		#Make assignments -1=ForceALL 0==non-assigned 1==core, 2==inter, 3==surf
		res_layer_assignments=np.zeros(connected_pose_mini_polyAA.total_residue(), int)
		for i in np.where(cb_exposition_ratio<core_cut_off)[0]:
			if(res_layer_assignments[i]==0):
				res_layer_assignments[i]=1 
				
		for i in np.where(cb_exposition_ratio<=boundary_cut_off)[0]:
			if(res_layer_assignments[i]==0):
				res_layer_assignments[i]=2
				
		for i in np.where(cb_exposition_ratio>boundary_cut_off)[0]:
			if(res_layer_assignments[i]==0):
				res_layer_assignments[i]=3
				
		#Here happens the assignment
		core_res=np.where(res_layer_assignments==1)[0]
		boundary_res=np.where(res_layer_assignments==2)[0]
		surface_res=np.where(res_layer_assignments==3)[0]
		
		this_pose_layer_dic_array=[]
		#Check correct size!!
		#D assert((len(core_res)+len(boundary_res)+len(surface_res)) == connected_pose_mini_desig.total_residue())

		for i in range(connected_pose_mini_desig.total_residue()):
			this_pose_layer_dic_array.append({})

		for i in core_res:
			this_pose_layer_dic_array[i]=layersDic.CoreAAdic

		for i in boundary_res:
			this_pose_layer_dic_array[i]=layersDic.BoundaryAAdic
			
		for i in surface_res:
			this_pose_layer_dic_array[i]=layersDic.SurfAAdic
		print "Done calculating pose layer booleans"
		
		print "Generating corresponding rosetta's packer boolenas"
		#Create a layer design task operation bools
		layer_design_task_bools = [] 

		for i in range(1, connected_pose_mini_desig.total_residue()+1):
			##print this_pose_layer_dic_array[i]
			aa_bool_tmp=[]
			for j in range( 1 , 21 ):
				aa_bool_tmp.append(False)
			for aa in  this_pose_layer_dic_array[i-1]:
				aa = self.rosetta.core.chemical.aa_from_oneletter_code( aa )
				aa_bool_tmp[int(aa)-1]=True
			##print i, this_pose_layer_dic_array[i-1]
			aa_bool = self.rosetta.utility.vector1_bool()
			for j in range( 1 , 21 ):
				aa_bool.append(aa_bool_tmp[j-1]==True)
			##layer_design_task.nonconst_residue_task(i).restrict_absent_canonical_aas( aa_bool )
			layer_design_task_bools.append(aa_bool)
		print "Generated %d, rosetta's booleans", len(layer_design_task_bools)
		print "Done generating rosetta's booleans"
		
		#Print layers to a pymol session and the respective pdbfile
		tmp_outname_pdb ="%s/%s_layer_design_input.pdb"%(general_out_dir, stage2_output_prefix)
		###print "Writing pymol graphical representation of layers:::"
		lineout="from pymol import cmd\n"
		lineout+="cmd.load('%s')\n"%tmp_outname_pdb
		lineout+="cmd.hide('all')\n"
		lineout+="cmd.show('sticks')\n"
		if len(core_res > 0):
			lineout+="cmd.select('core', 'resi "
			for i in core_res:
				lineout+="%d+"%(i+1)
			lineout+="')\n"
			lineout+="cmd.color('red', 'core')\n"
		if len(boundary_res > 0):
			lineout+="cmd.select('boundary', 'resi "
			for i in boundary_res:
				lineout+="%d+"%(i+1)
			lineout+="')\n"
			lineout+="cmd.color('orange', 'boundary')\n"
		if len(surface_res > 0):
			lineout+="cmd.select('surface', 'resi "
			for i in surface_res:
				lineout+="%d+"%(i+1)
			lineout+="')\n"
			lineout+="cmd.color('green', 'surface')\n"
		lineout+="cmd.select('empty', 'none')\n"
		lineout+="cmd.show('spheres')\n"
		lineout+="cmd.set('sphere_transparency', '0.7')\n"
		lineout+="cmd.show('cartoon')\n"
		##print lineout
		#Write Pymol Files
		
		tmp_outname_pyscript="%s/%s_layers_description.py"%(general_out_dir, stage2_output_prefix)
		print "Generating pymol layer files for visualization (so you know what is happening!) to files:\n  %s\n  %s" %(tmp_outname_pdb, tmp_outname_pyscript)
		connected_pose_mini_desig.dump_pdb(tmp_outname_pdb)
		#Write pymol python script for fraphical layers representation
		f = open(tmp_outname_pyscript,'w')
		f.write(lineout)
		f.close()
		
		print "If you want to see it...execute ~the-next-command in a terminal: \n   #> pymol %s" % tmp_outname_pyscript
		print "Done generating layers"
		return this_pose_layer_dic_array, layer_design_task_bools
		
		
	# replaces the residue at  <resid>  in  <pose>  with  <new_res>  with repacking
	def mutate_residue( self,
						pose , 
						mutant_positions , 
						mutant_aas , 
						in_task_bools,
						pack_radius = 0.0 ,
						 b_allow_design=False, 
						 pack_scorefxn = '' ):
		
		"""
		Replaces the residue at  <mutant_position>  in  <pose>  with  <mutant_aa>
			and repack any residues within  <pack_radius>  Angstroms of the mutating
			residue's center (nbr_atom) using  <pack_scorefxn>
		note: <mutant_aa>  is the single letter name for the desired ResidueType

		example:
			mutate_residue(pose,[30,31],[A,F])
		See also:
			Pose
			PackRotamersMover
			MutateResidue
			pose_from_sequence
		"""
		#### a MutateResidue Mover exists similar to this except it does not pack
		####    the area around the mutant residue (no pack_radius feature)
		#mutator = MutateResidue( mutant_position , mutant_aa )
		#mutator.apply( test_pose )

		#D assert (len(mutant_positions) == len(mutant_aas))
		#D assert (pose.total_residue() == len(in_task_bools))
		
		if pose.is_fullatom() == False:
			IOError( 'mutate_residue only works with fullatom poses' )

		test_pose = pose.clone()

		# create a standard scorefxn by default
		if not pack_scorefxn:
			pack_scorefxn = self.rosetta.core.scoring.ScoreFunctionFactory.create_score_function("talaris2013") 
		
		task = self.rosetta.standard_packer_task( test_pose )
		
		for i in range(1, len(in_task_bools)+1):
			task.nonconst_residue_task( i ).restrict_absent_canonical_aas( in_task_bools[i-1] )     
		
		b_mutable_positions=np.ones(pose.total_residue(), int)
		for i in mutant_positions:
			b_mutable_positions[i-1]=0
		
		for mutant_aa_i in range(len(mutant_positions)):
			mutant_aa=mutant_aas[mutant_aa_i]
			mutant_position=mutant_positions[mutant_aa_i]
			# the Vector1 of booleans (a specific object) is needed for specifying the
			#    mutation, this demonstrates another more direct method of setting
			#    PackerTask options for design
			aa_bool = self.rosetta.utility.vector1_bool()
			# PyRosetta uses several ways of tracking amino acids (ResidueTypes)
			# the numbers 1-20 correspond individually to the 20 proteogenic amino acids
			# aa_from_oneletter returns the integer representation of an amino acid
			#    from its one letter code
			# convert mutant_aa to its integer representation
			mutant_aa = self.rosetta.core.chemical.aa_from_oneletter_code( mutant_aa )
		
			# mutation is performed by using a PackerTask with only the mutant
			#    amino acid available during design
			# to do this, construct a Vector1 of booleans indicating which amino acid
			#    (by its numerical designation, see above) to allow
			for i in range( 1 , 21 ):
				# in Python, logical expression are evaluated with priority, thus the
				#    line below appends to aa_bool the truth (True or False) of the
				#    statement i == mutant_aa
				aa_bool.append( i == mutant_aa )
		
			# modify the mutating residue's assignment in the PackerTask using the
			#    Vector1 of booleans across the proteogenic amino acids
			task.nonconst_residue_task( mutant_position
				).restrict_absent_canonical_aas( aa_bool )

			# prevent residues from packing by setting the per-residue "options" of
			#    the PackerTask
			if (pack_radius > 0.01):
				center = pose.residue( mutant_position ).nbr_atom_xyz()
				for i in range(1,test_pose.total_residue()+1):
						if (b_mutable_positions[i-1]==1):
							if (center.distance_squared(test_pose.residue( i ).nbr_atom_xyz() ) < (pack_radius**2)):
								b_mutable_positions[i-1]=2
							
		
		for i in np.where(b_mutable_positions==1)[0]:
			task.nonconst_residue_task( i+1 ).prevent_repacking()
		if not b_allow_design:
			for i in np.where(b_mutable_positions==2)[0]:
				task.nonconst_residue_task( i+1 ).restrict_to_repacking()
		###print b_mutable_positions
		##print task
		##print test_pose.total_residue()
		
		# apply the mutation and pack nearby residues
		##print task
		packer = self.rosetta.PackRotamersMover( pack_scorefxn , task )
		packer.apply( test_pose )

		return test_pose
		
		
	# replaces the residue at  <resid>  in  <pose>  with  <new_res>  with repacking
	def mutate_residue_symmetric( self,
								pose, 
								symmetry_dictionary, 
								mutant_positions, 
								mutant_aas , 
								in_task_bools,
								pack_radius = 0.0 , 
								b_allow_design=False, 
								pack_scorefxn = '' ):
		
		"""
		Replaces the residue at  <mutant_position>  in  <pose>  with  <mutant_aa>
			and repack any residues within  <pack_radius>  Angstroms of the mutating
			residue's center (nbr_atom) using  <pack_scorefxn>
		note: <mutant_aa>  is the single letter name for the desired ResidueType

		example:
			mutate_residue(pose,[30,31],[A,F])
		See also:
			Pose
			PackRotamersMover
			MutateResidue
			pose_from_sequence
		"""
		#### a MutateResidue Mover exists similar to this except it does not pack
		####    the area around the mutant residue (no pack_radius feature)
		#mutator = MutateResidue( mutant_position , mutant_aa )
		#mutator.apply( test_pose )
		##print "A", mutant_positions
		##print "B", mutant_aas
		#D assert (len(mutant_positions) == len(mutant_aas))
		#D assert (pose.total_residue() == len(in_task_bools))
		
		tmp_mutant_positions=mutant_positions[:]
		tmp_mutant_aas=mutant_aas[:]
		for i in range(len(mutant_positions)):
			for j in symmetry_dictionary[mutant_positions[i]]:
				tmp_mutant_positions.append(j)
				tmp_mutant_aas=tmp_mutant_aas+mutant_aas[i]
			
		mutant_positions=tmp_mutant_positions[:]
		mutant_aas=tmp_mutant_aas[:]
		
		##print "C", mutant_positions
		##print "D", mutant_aas

		
		
		if pose.is_fullatom() == False:
			IOError( 'mutate_residue only works with fullatom poses' )

		test_pose = pose.clone()

		# create a standard scorefxn by default
		if not pack_scorefxn:
			pack_scorefxn = self.rosetta.core.scoring.ScoreFunctionFactory.create_score_function("talaris2013") 
		
		task = self.rosetta.standard_packer_task( test_pose )
		##for i in range(1, len(in_task_bools)+1):
		##    task.nonconst_residue_task( i ).restrict_absent_canonical_aas( in_task_bools[i-1] )     
		
		b_mutable_positions=np.ones(pose.total_residue(), int)
		for i in mutant_positions:
			b_mutable_positions[i-1]=0
		
		for mutant_aa_i in range(len(mutant_positions)):
			mutant_aa=mutant_aas[mutant_aa_i]
			mutant_position=mutant_positions[mutant_aa_i]
			
			# the Vector1 of booleans (a specific object) is needed for specifying the
			#    mutation, this demonstrates another more direct method of setting
			#    PackerTask options for design
			aa_bool = self.rosetta.utility.vector1_bool()
			# PyRosetta uses several ways of tracking amino acids (ResidueTypes)
			# the numbers 1-20 correspond individually to the 20 proteogenic amino acids
			# aa_from_oneletter returns the integer representation of an amino acid
			#    from its one letter code
			# convert mutant_aa to its integer representation
			mutant_aa = self.rosetta.core.chemical.aa_from_oneletter_code( mutant_aa )
		
			# mutation is performed by using a PackerTask with only the mutant
			#    amino acid available during design
			# to do this, construct a Vector1 of booleans indicating which amino acid
			#    (by its numerical designation, see above) to allow
			for i in range( 1 , 21 ):
				# in Python, logical expression are evaluated with priority, thus the
				#    line below appends to aa_bool the truth (True or False) of the
				#    statement i == mutant_aa
				aa_bool.append( i == mutant_aa )
		
			# modify the mutating residue's assignment in the PackerTask using the
			#    Vector1 of booleans across the proteogenic amino acids
			task.nonconst_residue_task( mutant_position
				).restrict_absent_canonical_aas( aa_bool )

			# prevent residues from packing by setting the per-residue "options" of
			#    the PackerTask
			if (pack_radius > 0.01):
				center = pose.residue( mutant_position ).nbr_atom_xyz()
				for i in range(1,test_pose.total_residue()+1):
						if (b_mutable_positions[i-1]==1):
							if (center.distance_squared(test_pose.residue( i ).nbr_atom_xyz() ) < (pack_radius**2)):
								b_mutable_positions[i-1]=2
							
		
		for i in np.where(b_mutable_positions==1)[0]:
			task.nonconst_residue_task( i+1 ).prevent_repacking()
		if not b_allow_design:
			for i in np.where(b_mutable_positions==2)[0]:
				task.nonconst_residue_task( i+1 ).restrict_to_repacking()
		##print task
		##assert(0==1)
		
		# apply the mutations 
		packer = self.rosetta.PackRotamersMover( pack_scorefxn , task )
		packer.apply( test_pose )

		return test_pose
		
		
	def generate_pose_from_fragment(self,
									in_angles, 
									b_add_nc_fake_caps=True):
		ndx_shift=-1
		tmpline=""
		if b_add_nc_fake_caps:
			ndx_shift=2
			for i in range(0, (len(in_angles)+2)):
				tmpline+=fake_aminoacid_letter_for_design
		else:
			ndx_shift=1
			for i in range(0, (len(in_angles))):
				tmpline+=fake_aminoacid_letter_for_design
		fragment_pose = self.rosetta.pose_from_sequence(tmpline)
		
		for i in range(0, len(in_angles)):
			fragment_pose.set_phi(i+ndx_shift, in_angles[i][0])
			fragment_pose.set_psi(i+ndx_shift, in_angles[i][1])
			fragment_pose.set_omega(i+ndx_shift, in_angles[i][2])
		return fragment_pose
		
		
	def calc_pose_sasa_total_per_residue_and_atom(self,
												in_pose, 
												probe_raddi=1.5):
		#Calculate SASA
		#D assert(in_pose.total_residue() > 0)
		atom_sasa = self.rosetta.core.id.AtomID_Map_Real()
		rsd_sasa = self.rosetta.utility.vector1_Real()
		#define atom_map for main-chain and CB
		atom_map = self.rosetta.core.id.AtomID_Map_bool()
		self.rosetta.core.pose.initialize_atomid_map(atom_map, in_pose, False)
		for ir in range( 1, in_pose.total_residue()+1 ):
				#Atoms 1 to 5     
				rsd = in_pose.residue( ir );
				for j in range (1, rsd.natoms()+1):
						atom = self.rosetta.core.id.AtomID(j, ir)
						atom_map.set( atom, True )
		
		total_SASA=self.rosetta.core.scoring.calc_per_atom_sasa(in_pose, 
													atom_sasa,
													rsd_sasa,
													probe_raddi,
													False,
													atom_map)
		
		per_res_sasa=[]
		for iVal in rsd_sasa:
			per_res_sasa.append(iVal)
		
		per_resAtom_sasa=[]
		for ir in range( 1, in_pose.total_residue()+1 ):
			per_resAtom_sasa.append({})
			rsd = in_pose.residue( ir );
			#Atoms 1 to 5                                                                   
			for j in range (1, 6):
				atom = self.rosetta.core.id.AtomID(j, ir)
				per_resAtom_sasa[ir-1][rsd.atom_name( j ).strip()]=atom_sasa[atom]
			
		return total_SASA, per_res_sasa, per_resAtom_sasa
		
		
	#Make mono-polyAA pose
	def make_polyAA_pose(self,
						in_pose, 
						targetAA='V'):
		pos_for_mut=range(1, in_pose.total_residue()+1)
		res_to_mut=[]
		task_bools_all_true=[]
		for i in range(in_pose.total_residue()):
			res_to_mut.append(targetAA)
			aa_bool_tmp=[]
			for j in range( 1 , 21 ):
				aa_bool_tmp.append(True)
			aa_bool = self.rosetta.utility.vector1_bool()
			for j in range( 1 , 21 ):
				aa_bool.append(aa_bool_tmp[j-1]==True)
			##layer_design_task.nonconst_residue_task(i).restrict_absent_canonical_aas( aa_bool )
			task_bools_all_true.append(aa_bool)
		connected_pose_mini_desig_polyAA=self.mutate_residue( in_pose , 
												 pos_for_mut , res_to_mut , task_bools_all_true )
		return connected_pose_mini_desig_polyAA
		
		
	#Make mono-polyAA pose or random AA pose
	def make_random_sequence_pose(self,
								in_pose, 
								possibleAA=None):
		if ( possibleAA == None ):
			possibleAA=['R','H','K','D','E','S','T','N','Q','A','V','I','L','M','F','Y','W'] #,'P','C','G'
		possibleAA=np.asarray(possibleAA)
		pos_for_mut=range(1, in_pose.total_residue()+1)
		res_to_mut=[]
		task_bools_all_true=[]
		for i in range(in_pose.total_residue()):
			targetAA=np.random.choice(possibleAA)
			res_to_mut.append(targetAA)
			aa_bool_tmp=[]
			for j in range( 1 , 21 ):
				aa_bool_tmp.append(True)
			aa_bool = self.rosetta.utility.vector1_bool()
			for j in range( 1 , 21 ):
				aa_bool.append(aa_bool_tmp[j-1]==True)
			##layer_design_task.nonconst_residue_task(i).restrict_absent_canonical_aas( aa_bool )
			task_bools_all_true.append(aa_bool)
		connected_pose_mini_desig_polyAA=mutate_residue( in_pose , 
												 pos_for_mut , res_to_mut , task_bools_all_true )
		return connected_pose_mini_desig_polyAA
		
		
	#Montecarlo with those sequences!
	def mc_optimize_pose_position(self,
									connected_pose_mini_copy_redes, 
									position, 
									best_possible_sequences_per_position, 
									in_task_bools, 
									sfx_insertion=None, 
									sfx_minimize=None,
									mc_iterations=100, 
									mc_reinsertions=10,
									mc_temp=1.0, 
									min_bb_each=10, 
									pack_radius=5.0,
									b_do_min=True,  
									b_do_final_min=True, 
									b_do_final_relax=True,
									do_allow_design_repacked=False,
									b_Debug=False):
		#Parameter: mutable positions
		mut_pos=range(position+1, position+1+4)
		#Initial Scoring
		connected_pose_mini_copy_redes.remove_constraints()
		startingE_score= sfx_insertion(connected_pose_mini_copy_redes)
		
		if b_Debug:
			print "Start Before, seq, targPos ", startingE_score, connected_pose_mini_copy_redes.sequence(), position

		best_score_until_now=startingE_score
		iteration=0
		b_mc_succed=False
		for iteration in range(1, mc_iterations):
			#Make a copy of the pose
			connected_pose_mini_copy_tmp=connected_pose_mini_copy_redes.clone()
			curr_seq=connected_pose_mini_copy_tmp.sequence()
		
			##Try insertion ~10 times
			for reinsertion_i in range( min( mc_reinsertions, len(best_possible_sequences_per_position) ) ):
				seq=best_possible_sequences_per_position[np.random.randint(0,len(best_possible_sequences_per_position))]
				connected_pose_mini_copy_tmp = self.mutate_residue( connected_pose_mini_copy_tmp ,mut_pos, seq, 
										in_task_bools, pack_radius = pack_radius, 
										b_allow_design=do_allow_design_repacked, pack_scorefxn = sfx_insertion ) 
				#Minimize sidechains after each fragment replacement
				if b_do_min:
					self.minimize_sidechains(connected_pose_mini_copy_tmp, sfx=sfx_minimize)
				tmp_score=sfx_insertion(connected_pose_mini_copy_tmp)
				#If energy is lower accept
				if(tmp_score < (startingE_score)):
					startingE_score=tmp_score
					connected_pose_mini_copy_redes=connected_pose_mini_copy_tmp.clone()
					b_mc_succed=True
					if (b_Debug):
						print position+1, curr_seq[position:position+4], "->", seq, startingE_score, connected_pose_mini_copy_redes.sequence()
					break
				else:
					#Add exp scaling here   e- (epsilon0 - epsilon1 )/KT
					##scale_probability_val=2 #1 ##This is like ~T :)
					if (mc_temp > 0.0): 
						e_val=(startingE_score-tmp_score)/mc_temp #(scipycons.Boltzmann* 200.0)
						boltz_val=np.exp(e_val)
						random_val=np.random.random()
						#print "Metro-vals: ", boltz_val, random_val, startingE_score-tmp_score
						if( (boltz_val > random_val) ):    ##(random_val>=0.95) and 
							if (b_Debug):
								print "Metropolis! ", startingE_score, "vs", tmp_score, "Metro-vals: ", boltz_val, random_val
							startingE_score=tmp_score
							connected_pose_mini_copy_redes=connected_pose_mini_copy_tmp.clone()
							b_mc_succed=True
							break             
			
			#print ( iteration, position+1, curr_seq[position:position+4], "->" , connected_pose_mini_copy_redes.sequence()[position:position+4], 
			#       startingE_score, connected_pose_mini_copy_redes.sequence() )
			
			"""
			if (((iteration % min_bb_each) == 0 )):
				if (b_Debug):
					print "Min all!"
				self.minimize_all_with_bb_constraints(connected_pose_mini_copy_redes, 
												 sfx=sfx_minimize,
												 harmonic_constraint_streght=1.0, 
												 cart_constraint_weight=0.2)
			startingE_score=sfx_insertion(connected_pose_mini_copy_redes)
			"""
			##print "Round M: ", iteration, startingE_score, connected_pose_mini_copy_redes.sequence()  
			
		#Minimization
		"""
		if (b_do_final_min):
			if (b_Debug):
				print "Min all!"
			self.minimize_all_with_bb_constraints(connected_pose_mini_copy_redes, 
												 sfx=sfx_minimize,
												 harmonic_constraint_streght=1.0, 
												 cart_constraint_weight=0.2)
		"""
		#Remove any constraints
		connected_pose_mini_copy_redes.remove_constraints()
		#Fast Real
		if (b_do_final_relax):
			relax_with_constraints_to_init(connected_pose_mini_copy_redes, sfx=sfx_minimize)
		
		startingE_score=sfx_minimize(connected_pose_mini_copy_redes)
		###connected_pose_mini_copy_redes.dump_pdb("%s/tmp_optimizations/test_99_final_fastr_struct.pdb"%general_out_dir)
		
		if (b_Debug):
			print "Final after: ", iteration, startingE_score, connected_pose_mini_copy_redes.sequence() 
			print "\n\nDONE!!!\n\n"
			
		return connected_pose_mini_copy_redes, startingE_score, b_mc_succed
		
		
	#Montecarlo with those sequences!
	def mc_optimize_pose_position_symmetric(self,
									connected_pose_mini_copy_redes, 
									symmetry_dictionary,
									position, 
									best_possible_sequences_per_position, 
									in_task_bools, 
									sfx_insertion=None, 
									sfx_minimize=None, 
									mc_iterations=100, 
									mc_reinsertions=10,
									mc_temp=1.0, 
									min_bb_each=10, 
									pack_radius=5.0,
									b_do_min=True,  
									b_do_final_min=True,
									b_do_final_relax=True,
									b_Debug=False):
		#Parameter: mutable positions
		mut_pos=range(position+1, position+1+4)
		#Initial Scoring
		connected_pose_mini_copy_redes.remove_constraints()
		startingE_score= sfx_insertion(connected_pose_mini_copy_redes)
		
		if b_Debug:
			print "Start Before, seq, targPos ", startingE_score, connected_pose_mini_copy_redes.sequence(), position
		
		best_score_until_now=startingE_score
		iteration=0
		b_mc_succed=False
		for iteration in range(1, mc_iterations):
			#Make a copy of the pose
			connected_pose_mini_copy_tmp=connected_pose_mini_copy_redes.clone()
			curr_seq=connected_pose_mini_copy_tmp.sequence()
			#Select randomly the replace position
			##i = position ##random.randint(0, len(best_possible_sequences_per_position))
		
			##Try insertion ~10 times
			for reinsertion_i in range( min( mc_reinsertions, len(best_possible_sequences_per_position) ) ):
				seq=best_possible_sequences_per_position[np.random.randint(0,len(best_possible_sequences_per_position))]
				connected_pose_mini_copy_tmp = mutate_residue_symmetric( connected_pose_mini_copy_tmp, symmetry_dictionary, 
														mut_pos, seq, in_task_bools, pack_radius = pack_radius, 
														b_allow_design=False, pack_scorefxn = sfx_insertion ) 
				#Minimize sidechains after each fragment replacement
				if b_do_min:
					self.minimize_sidechains(connected_pose_mini_copy_tmp, sfx=sfx_minimize)
				tmp_score=sfx_insertion(connected_pose_mini_copy_tmp)
				#If energy is lower accept
				if(tmp_score < (startingE_score)):
					startingE_score=tmp_score
					connected_pose_mini_copy_redes=connected_pose_mini_copy_tmp.clone()
					b_mc_succed=True
					if (b_Debug):
						print position+1, curr_seq[position:position+4], "->", seq, startingE_score, connected_pose_mini_copy_redes.sequence()
					break
				else:
					#Add exp scaling here  e- (epsilon0 - epsilon1 )/KT
					##scale_probability_val=2 #1 ##This is like ~T :)
					if (mc_temp > 0.0):
						e_val=(startingE_score-tmp_score)/mc_temp #(scipycons.Boltzmann* 200.0)
						boltz_val=np.exp(e_val)
						random_val=np.random.random()
						#print "Metro-vals: ", boltz_val, random_val, startingE_score-tmp_score
						if( (boltz_val > random_val) ):    ##(random_val>=0.95) and 
							if (b_Debug):
								print "Metropolis! ", startingE_score, "vs", tmp_score, "Metro-vals: ", boltz_val, random_val
							startingE_score=tmp_score
							connected_pose_mini_copy_redes=connected_pose_mini_copy_tmp.clone()
							b_mc_succed=True
							break             
			
			#print ( iteration, position+1, curr_seq[position:position+4], "->" , connected_pose_mini_copy_redes.sequence()[position:position+4], 
			#       startingE_score, connected_pose_mini_copy_redes.sequence() )
			
			
			if (((iteration % min_bb_each) == 0 )):
				if (b_Debug):
					print "Min all!"
				self.minimize_all_with_bb_constraints(connected_pose_mini_copy_redes, 
												 sfx=sfx_minimize,
												 harmonic_constraint_streght=1.0, 
												 cart_constraint_weight=0.2)
			startingE_score=sfx_insertion(connected_pose_mini_copy_redes)
			##print "Round M: ", iteration, startingE_score, connected_pose_mini_copy_redes.sequence()  
			
		#Minimization
		if (b_do_final_min):
			if (b_Debug):
				print "Min all!"
			self.minimize_all_with_bb_constraints(connected_pose_mini_copy_redes, 
												 sfx=sfx_minimize,
												 harmonic_constraint_streght=1.0, 
												 cart_constraint_weight=0.2)
		#Remove any constraints
		connected_pose_mini_copy_redes.remove_constraints()
		#Fast Real
		if (b_do_final_relax):
			relax_with_constraints_to_init(connected_pose_mini_copy_redes, sfx=sfx_minimize)
		
		startingE_score=sfx_minimize(connected_pose_mini_copy_redes)
		###connected_pose_mini_copy_redes.dump_pdb("%s/tmp_optimizations/test_99_final_fastr_struct.pdb"%general_out_dir)
		
		if (b_Debug):
			print "Final after: ", iteration, startingE_score, connected_pose_mini_copy_redes.sequence() 
			print "\n\nDONE!!!\n\n"
			
		return connected_pose_mini_copy_redes, startingE_score, b_mc_succed
		
		
	def calculate_possible_aa_by_res_using_fragments(self,
										connected_pose_mini_desig=None,
										this_pose_layer_dic_array=None,
										clustersDB=None,
										min_tar_num_seq_match=10,
										to_center_cut_off_dist=0.8,
										general_out_dir="./",
										stage2_output_prefix="stage2_"):
		#Calculate posible aa matching the sequence profile
		#Settings
		print "Calculating optimum PIKAAs, this might take a while, be patient."
		##End Settings
		optimized_assignments_design = self.assign_pose_to_clusters_within_cutoff(connected_pose_mini_desig, 
																			1, 
																			connected_pose_mini_desig.total_residue(),
																			cut_off=to_center_cut_off_dist,
																			clustersDB=clustersDB,
																			permissive=False)

		optimized_pikaa=[]
		for iPos in range(connected_pose_mini_desig.total_residue()):
			optimized_pikaa.append({})
			
		for iPos in range(connected_pose_mini_desig.total_residue()-3):
			#print "Calculating opt PIKAA for pos:", iPos, "of:", connected_pose_mini_desig.total_residue()-3
			#Find the ndx for the corresponding assignments
			corresponding_assignment_ndxs=[]
			for cluster in optimized_assignments_design[iPos]:
				for ndx in clustersDB.aa_assignments_ndx_dic[cluster]:
					corresponding_assignment_ndxs.append(ndx)
					
			#Find the corresponding sequences and make a dictionary
			possible_sequences = clustersDB.assignments_lr_store_data["aaSequence"][corresponding_assignment_ndxs]
			#print "Analyzing #seq: ", len(possible_sequences)
			
			
			possible_sequences_dic={}
			for i in possible_sequences:
				if not i in possible_sequences_dic:
					possible_sequences_dic[i]=1
				else:
					possible_sequences_dic[i]+=1
					
			tmp_seq_array=[]
			for seq in possible_sequences_dic:
				if (clustersDB.aa_seq_dic[seq][0] >= min_tar_num_seq_match):
					tmp_seq_array.append(seq)
			tmp_seq_array=np.asarray(tmp_seq_array)

			#Prune by layer
			tmp_frag_layer_dic_array = this_pose_layer_dic_array[iPos:iPos+4] 
			best_possible_sequences=[]
			for iSeq in range(len(tmp_seq_array)):
				if(compare_seq_2polar_profile(tmp_seq_array[iSeq], tmp_frag_layer_dic_array)):
					best_possible_sequences.append(tmp_seq_array[iSeq])

			for seq in best_possible_sequences:
				for jPos in range(4):            
					if seq[jPos] not in optimized_pikaa[(iPos+jPos)]:
						optimized_pikaa[(iPos+jPos)][seq[jPos]]=1
					else:
						optimized_pikaa[(iPos+jPos)][seq[jPos]]+=1
						
		#for i in range( len(optimized_pikaa) ):
		#    print i, optimized_pikaa[i]
		print "Done calculating optimum PIKAAs,"
		
		
		#Now out a simple layer resfile for rosetta
		outfilename="%s/%s_layers_resfile.resfile"%(general_out_dir, stage2_output_prefix)
		print "Generating a general layer profile resfile to:", outfilename
		aa_line_out=[]
		for i in range(connected_pose_mini_desig.total_residue()):
			tmplineout="%d A PIKAA " %(i+1)
			for aa in optimized_pikaa[i]:
				tmplineout += aa
			aa_line_out.append(tmplineout)
		lines_out = "ALLAA # This default command applies to all residues that are not given specific commands\nstart\n"
		for line in aa_line_out:
			lines_out+=line+"\n"
			#print line
		f = open(outfilename,'w')
		f.write(lines_out)
		f.close()
		print "DONE"
		
		#Print a pose a resfile enforcing our optimization on loopsand restricting everything else to the previously calculated layers
		outfilename="%s/%s_layers_lockedLoops_resfile.resfile"%(general_out_dir, stage2_output_prefix)
		print "Generating a loop-locking layer profile resfile:", outfilename
		tmp_dssp = self.rosetta.core.scoring.dssp.Dssp(connected_pose_mini_desig)
		tmp_dssp_ss =tmp_dssp.get_dssp_secstruct()
		tmp_dssp_ss=np.asarray(list(tmp_dssp_ss))
		print "DSSP predictio to be used to Lock loop (L):", tmp_dssp_ss
		aa_line_out=[]
		for i in range(connected_pose_mini_desig.total_residue()):
			tmplineout="%d A PIKAA " %(i+1)
			for aa in optimized_pikaa[i]:
				tmplineout += aa
			aa_line_out.append(tmplineout)
			
		final_aaseq=connected_pose_mini_desig.sequence()
		tmp_ndx=np.where(tmp_dssp_ss=='L')[0]
		print "Lock loop (L) indexes in resfile ", tmp_ndx
		print "First and last 2 residues have been avoided to be locked in purpose"
		for i in tmp_ndx:
			if ( (i > 1) and (i < connected_pose_mini_desig.total_residue()-2) ):
				aa_line_out[i] = "%d A NATAA" % (i+1) #%s" %(i+1, final_aaseq[i])
		lines_out = "ALLAA # This default command applies to all residues that are not given specific commands\nstart\n"
		for line in aa_line_out:
			lines_out+=line+"\n"
			##print line
		f = open(outfilename,'w')
		f.write(lines_out)
		f.close()
		print "Done generating lresfiles, you can use them later with rosetta if you like."
		
		
	def optimize_sequence_using_clustered_fragments_stats(self,
											in_pose=None,
											clustersDB=None,
											this_pose_layer_dic_array=None,
											layer_design_task_bools=None,
											design_quality_entropy_cut_off=100000.0, #Adaptive initial target quality ratio
											min_target_design_quality_entropy_cut_off=2.0, #Adaptive final target quality ratio
											design_quality_entropy_look_within_percentage=0.1, #in the scale of 1.0 to 0.0
											out_cycles=20, #AutoOptOption
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
											num_cycles_without_improvement_tol=10, #(in_cycles*2), #FineTuneOption
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
											general_out_dir="./",
											stage2_tmp_dir_name="stage2_seqOptimization_tmp",
											locked_aa_indexes_identities_array=[], #Opt. ex: locked_aa_indexes_identities_array.append([2, 'E'])
											b_lock_coreSS_aa=False, #Do lock a.a. that belong to SScore during the optimization?
											b_use_symmetric_packing=False,
											b_Debug=False  ):
		#ToDo: Missing to implement options:
		#mc_temperature_minimization_start=100.0
		
		#DoNoT modify the original input, use a copy"
		connected_pose_mini_desig=in_pose.clone()
		
		#Score average and output
		tmp_score = self.scorefxn_tal(connected_pose_mini_desig)
		connected_pose_mini_desig.dump_pdb("%s/%s/merged_test_phase1_designed_01_optStart.pdb"%(general_out_dir,stage2_tmp_dir_name))
		target_score_design_reference=tmp_score
		print "Initial score: ", target_score_design_reference, tmp_score
		
		#Take care of user locked a.a. positions
		#-1 because our internal convention starts numbering from 0
		for i in range(len(locked_aa_indexes_identities_array)):
			locked_aa_indexes_identities_array[i][0]=locked_aa_indexes_identities_array[i][0]-1

		global_avoid_refine_positions=[]
		for aadefinition in (locked_aa_indexes_identities_array):
			#Remember: Rosetta seq pos convention is +1
			#print aadefinition[0]+1, aadefinition[1]
			connected_pose_mini_desig = self.mutate_residue( connected_pose_mini_desig , [(aadefinition[0]+1)] , [aadefinition[1]] ,
										[], pack_radius=0.0, 
										b_allow_design=False, pack_scorefxn=self.scorefxn_tal )
			global_avoid_refine_positions.append(aadefinition[0])
		print "Will avoid to design locked global positions: ", global_avoid_refine_positions
		print "New sequence after user enforced aa(s): ", connected_pose_mini_desig.sequence()
		
		#Now use DSSP to add constraints to all the SS if desired
		
		global_avoid_refine_ss_at=[]
		if b_lock_coreSS_aa:
			tmp_avoid_refine_at=[]
			tmp_dssp = self.rosetta.core.scoring.dssp.Dssp(connected_pose_mini_desig)
			tmp_dssp_ss =tmp_dssp.get_dssp_secstruct()
			for i in range(connected_pose_mini_desig.total_residue()-3):
				if(tmp_dssp_ss[i] != 'L'):
					tmp_avoid_refine_at.append(i)
			  
			#Fill gaps
			tmp_avoid_refine_at_2=[]
			for i in range(len(tmp_avoid_refine_at)-1):
				if(tmp_avoid_refine_at[i] >= (tmp_avoid_refine_at[i+1]-3)):
					for j in range(tmp_avoid_refine_at[i], tmp_avoid_refine_at[i+1]):
						tmp_avoid_refine_at_2.append(j)
			tmp_avoid_refine_at=tmp_avoid_refine_at_2[:]
				
			#shrink by +- 1
			isFirst=True
			for i in range(len(tmp_avoid_refine_at)-1):
				if isFirst and (tmp_avoid_refine_at[i] >= (tmp_avoid_refine_at[i+1]-1)):
					isFirst=False
				elif ((not isFirst) and (tmp_avoid_refine_at[i] >= (tmp_avoid_refine_at[i+1]-1))): 
					global_avoid_refine_ss_at.append(tmp_avoid_refine_at[i])
				else:
					isFirst=True
					
			print "Just an informative pymol selection for you (Do not design at):", global_avoid_refine_ss_at
			pymolline="select resi "
			for i in global_avoid_refine_ss_at:
				pymolline=pymolline+str(i)+"+"
			print pymolline
			
		#Automatic sequence optimizator based on clustered fragments information
		#ToDo: ADD SS based restriction of Prolines (i.e., for example, cannot have them in the middle of sheets)
		#ToDo: Improve packing, not good enough
		#Improve relaxation, now it can drift from input too much
		print "Starting complicated regions sequence optimization. Be patient"
		if b_use_symmetric_packing:
			print "NOTE: You have activated Pseudo-symmetric packing"
		else:
			print "Using not symmetring packing"
			
		list_of_optimized_positions=[]

		#mc_temperature_minimization_start=100.0

		#Make a container to store the best result.
		connected_pose_mini_desig_best_optimization=connected_pose_mini_desig.clone()
		best_worst_ratio_value=np.finfo('d').max
		##best_worst_ratio_cycle_num=-1
		best_optimization_list_of_optimized_positions=[]

		#Step1
		#Assign structure to clusters within certain cutoff
		tmp_assignments_design = self.assign_pose_to_clusters_within_cutoff(connected_pose_mini_desig, 
																			1, 
																			connected_pose_mini_desig.total_residue(), 
																			cut_off=target_max_cluster_distance, 
																			clustersDB=clustersDB,
																			permissive=False)

		#Booleans for trimed sequences
		alrready_tried_seq=[]
		for i in range(connected_pose_mini_desig.total_residue()):
			alrready_tried_seq.append({})
			
		#Main optimization loop
		num_cycles_without_improvement=0;
		woirst_position_ndx=-1
		#Save a reference (will need it later)
		connected_pose_mini_desig_initial_reference=in_pose.clone()
		best_outer_design_cycle_pose=in_pose.clone()
		best_total_score=self.scorefxn_tal(connected_pose_mini_desig_initial_reference)
		
		system_total_residue=connected_pose_mini_desig.total_residue()
		per_res_total_energies=np.zeros(system_total_residue, float)
		b_allow_change_restricted_positions=False
		opt_cycle_number=-1
		for outer_optimization_cycle in range(out_cycles):
			avoid_refine_positions=[]
			#Step1
			#re-Assign structure to clusters within certain cutoff
			if b_reassign_each_out_cycle:
				tmp_assignments_design = self.assign_pose_to_clusters_within_cutoff(connected_pose_mini_desig, 
																					1, 
																					connected_pose_mini_desig.total_residue(), 
																					cut_off=target_max_cluster_distance, 
																					clustersDB=clustersDB,
																					permissive=False)
			for inner_optimization_cycle in range(in_cycles):
				opt_cycle_number+=1
				print outer_optimization_cycle, "-", inner_optimization_cycle

				#ToDo: maybe only do this on energy optimization cycles
				#Get per residue Energies for minimization:
				curr_score = self.scorefxn_tal(connected_pose_mini_desig)
				print "E curr score: ", curr_score
				for i in range(system_total_residue):
					per_res_total_energies[i] = connected_pose_mini_desig.energies().residue_total_energy(i+1)
					
				#Step2
				#ToDO: modify to use dictionaries
				#First figure out how many fragments correspond to the current sequence of each fragment
				#Calculate the number of fragments
				tmp_aaseq=connected_pose_mini_desig.sequence()
				print "In a.a. sequence: ", tmp_aaseq
				num_fragment_sequences_match_per_position = np.zeros((len(tmp_aaseq)-3,2), int)
				for i in range(len(tmp_aaseq)-3):
					num_fragment_sequences_match_per_position[i,0]=i
					try:
						num_fragment_sequences_match_per_position[i,1]= clustersDB.aa_seq_dic[tmp_aaseq[i:i+4]][0] 
					except KeyError, e:
						num_fragment_sequences_match_per_position[i,1]=0
						print "Error: The sequence ", tmp_aaseq[i:i+4], "doesn't exist in the database.. skipping"
						
				#Now plot
				##fig = plt.figure(figsize=(15,10))
				##ax = fig.add_subplot(2,1,1)
				##ax.plot(num_fragment_sequences_match_per_position[:,0], -num_fragment_sequences_match_per_position[:,1], c='blue')
				
				#Step3
				#ToDO: modify to use dictionaries
				#Now figure out if our clusters contain such sequences per fragment and how many of them
				tmp_aaseq=connected_pose_mini_desig.sequence()
				num_fragment_sequences_matching_structure = np.zeros((len(tmp_aaseq)-3,2), int)
				for ipos in range(len(tmp_assignments_design)):
					#Remember to -1
					try:
						tmp_clusters_seq=(clustersDB.assignments_lr_store_data["assignment"][clustersDB.aa_seq_2ndx_dic[tmp_aaseq[ipos:ipos+4]]]-1)  
					except KeyError, e:
						print "Error: The sequence ", tmp_aaseq[ipos:ipos+4], "doesn't exist in the database.. skipping"
						tmp_clusters_seq=[]
					 
					num_fragment_sequences_matching_structure[ipos,0]=ipos
					num_fragment_sequences_matching_structure[ipos,1]=np.in1d(tmp_clusters_seq, (tmp_assignments_design[ipos])).sum()
				#Now plot total seq vs contained sequences
				if (b_generate_quality_plot):
					fig = plt.figure(figsize=(15,10))
					ax = fig.add_subplot(2,1,1)
					ax.axis([num_fragment_sequences_match_per_position[:,0].min(), 
							 num_fragment_sequences_match_per_position[:,0].max(), 
							 -500, 0.5])
					ax.plot(num_fragment_sequences_match_per_position[:,0], -num_fragment_sequences_match_per_position[:,1], c='blue')   
					ax.plot(num_fragment_sequences_matching_structure[:,0], -num_fragment_sequences_matching_structure[:,1], c='red')
					fig.savefig("%s/%s/plot_diff_c%04d_%04d_p%04d.png"% (general_out_dir,
																		stage2_tmp_dir_name,
																		outer_optimization_cycle,
																		inner_optimization_cycle,
																		(woirst_position_ndx+1)) ) 
				
				#Step4
				#Now plot the difference (This is entropy at 1kT)
				small_val_pert=0.00001
				seq_entropy_ratio = (num_fragment_sequences_match_per_position[:,1])/(
													num_fragment_sequences_matching_structure[:,1]+small_val_pert)
				##print seq_entropy_ratio
				##assert 0==1
				#Plot ratio_entropy
				if (b_generate_quality_plot):
					fig = plt.figure(figsize=(15,10))
					ax = fig.add_subplot(2,1,1)
					ax.axis([num_fragment_sequences_match_per_position[:,0].min(), 
							 num_fragment_sequences_match_per_position[:,0].max(), 
							 0.99, 100])
					ax.set_yscale('log')
					ax.plot(num_fragment_sequences_match_per_position[:,0], seq_entropy_ratio, c='blue') 
					fig.savefig("%s/%s/plot_qual_c%04d_%04d_p%04d.png"% (general_out_dir,
																		stage2_tmp_dir_name,
																		outer_optimization_cycle,
																		inner_optimization_cycle,
																		(woirst_position_ndx+1)) )
				#Plot per residue energy
				if (b_generate_quality_plot):
					fig = plt.figure(figsize=(15,10))
					ax = fig.add_subplot(2,1,1)
					tmp_resids=np.asarray(range(len(per_res_total_energies)))+1
					ax.axis([tmp_resids.min(), 
							 tmp_resids.max(), 
							 -5, #per_res_total_energies.min(), 
							 +5]) #per_res_total_energies.max()])
					ax.set_yscale('linear')
					ax.plot(tmp_resids, per_res_total_energies, c='blue') 
					fig.savefig("%s/%s/plot_rener_c%04d_%04d_p%04d.png"% (general_out_dir,
																		stage2_tmp_dir_name,
																		outer_optimization_cycle,
																		inner_optimization_cycle,
																		(woirst_position_ndx+1)) )
				  
				#This will avoid ploting in the next round if there is no improvement
				##b_generate_quality_plot=False
				
				#Step5
				tmp_seq_entropy_ratio = np.copy(seq_entropy_ratio)
				
				#Debug
				##for i in range(len(tmp_seq_entropy_ratio)):
				##    print (i, tmp_seq_entropy_ratio[i], num_fragment_sequences_match_per_position[:,1][i], 
				##        num_fragment_sequences_matching_structure[:,1][i]+small_val_pert)
				#END DEBUG
					
				#Avoid designing user constrained positions
				###print "Avoid global positions", global_avoid_refine_positions
				for pos in global_avoid_refine_positions:
					#for j in range(4):
					tmp_seq_entropy_ratio[pos]=-1.0
					per_res_total_energies[pos]=-9999.0
						
				#Avoid designing SS???
				###print "Avoid global SSs", global_avoid_refine_ss_at
				for pos in global_avoid_refine_ss_at:
					tmp_seq_entropy_ratio[pos]=-1.0
					per_res_total_energies[pos]=-9999.0
				
				
				#Asses the target quality
				if (tmp_seq_entropy_ratio.max() < design_quality_entropy_cut_off):
					tmp_delta=1.0
					#tmp_new_tar_quality=min((design_quality_entropy_cut_off-tmp_delta), tmp_seq_entropy_ratio.max())
					tmp_new_tar_quality=tmp_seq_entropy_ratio.max()-tmp_delta
					if (tmp_new_tar_quality < min_target_design_quality_entropy_cut_off):
						tmp_new_tar_quality=min_target_design_quality_entropy_cut_off
					print ("Target quality reached:", design_quality_entropy_cut_off, 
						   "Max now:",  tmp_seq_entropy_ratio.max(),
						   " increasing stringency to ", tmp_new_tar_quality)
					design_quality_entropy_cut_off=tmp_new_tar_quality
				else:
					tmp_new_tar_quality=design_quality_entropy_cut_off
					tmp_new_tar_quality-=0.05
					if (tmp_new_tar_quality < min_target_design_quality_entropy_cut_off):
						tmp_new_tar_quality=min_target_design_quality_entropy_cut_off
					print "Auto increasing qual stringency from, to:", design_quality_entropy_cut_off, tmp_new_tar_quality
					design_quality_entropy_cut_off=tmp_new_tar_quality
					
					
						 
				#Step6
				#Avoid designing previously designed positions
				###print "Avoid previously optimized positions", avoid_refine_positions
				##if False:
				for i in avoid_refine_positions:
					tmp_seq_entropy_ratio[i]=-1.0
					
						
				#Step 6b
				#Avoid designing regions that are alrready good
				#ToDo: Remove the "span stuff, it is not usefull.!
				avoid_alrready_optimal_positions=[]
				optimizable_positions=[]
				tmp_spanned_design_quality_entropy_cut_off=(design_quality_entropy_cut_off-
										(design_quality_entropy_cut_off*design_quality_entropy_look_within_percentage))
				for i in range(len(tmp_seq_entropy_ratio)):
					tmp_low_span=max(0,i-1)
					tmp_hi_span=min(i+1, len(tmp_seq_entropy_ratio))
					tmp_hi_span+=1
					#decide if work on higest energy or higest seq disagreement
					if b_allow_change_restricted_positions:
						#Allow all
						#Remove -1 and -9999 hack
						#Get positions worst than average
						if(per_res_total_energies[i] > -9999.0):
								optimizable_positions.append(i)
						else:
							avoid_alrready_optimal_positions.append(i)
							tmp_seq_entropy_ratio[i]=-1.0
					else:
						if (tmp_seq_entropy_ratio[tmp_low_span:tmp_hi_span].max() < tmp_spanned_design_quality_entropy_cut_off):
							avoid_alrready_optimal_positions.append(i)
							tmp_seq_entropy_ratio[i]=-1.0
						else:
							optimizable_positions.append(i)
						
				optimizable_positions=np.asarray(optimizable_positions)
				##assert(0==1)
						
				###if (len(avoid_alrready_optimal_positions) > 0):
				###    print "Avoiding alrready optimal positions: ", avoid_alrready_optimal_positions
				
				print "At quality :", b_allow_change_restricted_positions, design_quality_entropy_cut_off, "-", tmp_spanned_design_quality_entropy_cut_off, ", the optimizable positions: ", optimizable_positions+1
				if( len(optimizable_positions) == 0 ):
					print "NOTE: All designable positions are alrready optimal, EnForcing END of loop NOW!!!"
					break
				##print tmp_seq_entropy_ratio
				##assert(0==1)
				#Identify where is the worst problem, but ?avoid repetition based on previous worked positions?
				#ToDo, also identify those positions with too few sequences
				if b_allow_change_restricted_positions:
					print "Rand-Choosing worst position by energy"
					#This is energy minimization
					tmp_energy_average=np.average(per_res_total_energies[optimizable_positions])
					print "Looking for worst residue position by energy > :", tmp_energy_average, "the optimizable positions: ", optimizable_positions+1
					woirst_position_ndx=np.random.choice(optimizable_positions)
					#Limit up to the last 4mer fragment (n-3)
					woirst_position_ndx=min(len(per_res_total_energies)-3,woirst_position_ndx)
					print "TEST, position with worst energy: ", woirst_position_ndx+1, per_res_total_energies[woirst_position_ndx]
				else:
					#This is fragment quality improvement
					print "rand-Choosing worst position by quality"
					woirst_position_ndx=np.random.choice(optimizable_positions)
				#Append to the list to avoid refining during the inner cycle
				avoid_refine_positions.append(woirst_position_ndx)
				
				##woirst_position_ndx=tmp_seq_entropy_ratio.argmax()
				woirst_position_ratio=(1.0/tmp_seq_entropy_ratio[woirst_position_ndx])
				##if b_Debug:
				print ("Now working in pos/seq/ratio: ", (woirst_position_ndx+1), 
					   tmp_aaseq[woirst_position_ndx:woirst_position_ndx+4], 
					   "Ratio: ", woirst_position_ratio)
				#Add such position to the list so we do not re try it in the next internal iteration:
				#Don't move this further down in the code!!!
				avoid_refine_positions.append(woirst_position_ndx)
				
				
				#Find the ndx for the corresponding assignments
				corresponding_assignment_ndxs=[]
				for cluster in tmp_assignments_design[woirst_position_ndx]:
					for ndx in clustersDB.aa_assignments_ndx_dic[cluster]:
						corresponding_assignment_ndxs.append(ndx)
						
				#Find the corresponding sequences and make a dictionary
				possible_sequences = (clustersDB.assignments_lr_store_data["aaSequence"][corresponding_assignment_ndxs])
				if b_Debug:
					print "Analyzing #seq: ", len(possible_sequences)
				
				possible_sequences_dic={}
				for i in possible_sequences:
					if not i in possible_sequences_dic:
						possible_sequences_dic[i]=1
					else:
						possible_sequences_dic[i]+=1
			
				#possible_sequences_unique = np.unique(clustersDB.assignments_lr_store_data["aaSequence"][corresponding_assignment_ndxs])
				##print "Analyzing #seqU: ", len(possible_sequences_dic)
				
				#NOTE: Best possible sequence are those with many sequences here but in which the sequence 
				##is not too represented outside of these cluster
				##min_target_quality=woirst_position_ratio
				#while ( (len(best_possible_sequences) < 10) and (target_quality > 0.10 ) ):

				#ToDo: Make thise numbers adaptive
				min_target_quality=max((woirst_position_ratio),tar_seq_match_min_ratio) #i.e. at least increase by double??
				##if b_Debug:
				print "Looking for optimizing sequences with min seq/qual/temp: ", tar_num_seq_match, min_target_quality, mc_temperature
				b_generate_quality_plot=False
				min_number_seq_2_pass_phase1=1
				min_number_seq_2_pass_phase_2=1
				top_best_possible_sequences=[]
				cut_ratio=0.9
				
				##print "****BLASHA: ", alrready_tried_seq
				#Fix this while
				best_possible_sequences=[]
				best_possible_sequences_ratios=[]
				tmp_seq_array_vals=[]
				tmp_seq_array=[]
				for seq in possible_sequences_dic:
					#Avoid repeating the same sequence that is already there and previously tried sequences
					if(seq == tmp_aaseq[woirst_position_ndx:woirst_position_ndx+4]):
						continue
					elif ( (not allow_frag_seq_retry) and 
						  (seq in alrready_tried_seq[woirst_position_ndx]) ):
						##print "REPEATED!!!", seq
						continue
					if (clustersDB.aa_seq_dic[seq][0] >= tar_num_seq_match):
						#Calculate the ratio
						tmp_ratio_val=(float(possible_sequences_dic[seq])/clustersDB.aa_seq_dic[seq][0])
						#Debug
						tmp_seq_array_vals.append(tmp_ratio_val)
						tmp_seq_array.append(seq)
						#End Debug
						#if( (tmp_ratio_val >= target_quality) and (tmp_ratio_val >= woirst_position_ratio )):
						#    tmp_frag_layer_dic_array = this_pose_layer_dic_array[woirst_position_ndx:woirst_position_ndx+4] 
						#    if(compare_seq_2polar_profile(seq, tmp_frag_layer_dic_array)):
						#    ##if(True):
						#        best_possible_sequences.append(seq)
						#        best_possible_sequences_ratios.append(tmp_ratio_val)
				
				tmp_seq_array_vals=np.asarray(tmp_seq_array_vals)
				tmp_seq_array=np.asarray(tmp_seq_array)
				tmp_ndx_best=np.where(tmp_seq_array_vals >= min_target_quality)
				tmp_seq_array_best = tmp_seq_array[tmp_ndx_best]
				tmp_seq_array_vals_best =tmp_seq_array_vals[tmp_ndx_best]

				best_possible_sequences=[]
				best_possible_sequences_ratios=[]
				tmp_frag_layer_dic_array = this_pose_layer_dic_array[woirst_position_ndx:woirst_position_ndx+4] 
				###print tmp_frag_layer_dic_array
				for i_seq in range(len(tmp_seq_array_best)):
					if(compare_seq_2polar_profile(tmp_seq_array_best[i_seq], tmp_frag_layer_dic_array)):
						best_possible_sequences.append(tmp_seq_array_best[i_seq])
						best_possible_sequences_ratios.append(tmp_seq_array_vals_best[i_seq])
				###print "SAIUHGAISUH", best_possible_sequences, best_possible_sequences_ratios
				best_possible_sequences=np.asarray(best_possible_sequences)
				best_possible_sequences_ratios=np.asarray(best_possible_sequences_ratios)


				##print "Warning A: not improving sequences found.. skipping to next", target_quality, tar_num_seq_match, len(best_possible_sequences)
				if ( (len(best_possible_sequences) < min_number_seq_2_pass_phase1) ):
					b_generate_quality_plot=False
					if b_Debug:
						print "Warning A: not improving sequences found.. skipping to next position"
					avoid_refine_positions.append(woirst_position_ndx)
					if (mc_temperature > 0):
						mc_temperature+=5
						if (mc_temperature > max_mc_temperature):
							mc_temperature=max_mc_temperature
					#Experimental adaptive expected num seq change
					tar_seq_match_min_ratio-=0.05
					if (tar_seq_match_min_ratio < min_allowed_tar_seq_match_min_ratio):
						tar_seq_match_min_ratio=min_allowed_tar_seq_match_min_ratio
					tar_num_seq_match-=5
					if (tar_num_seq_match < min_allowed_tar_num_seq_match):
						tar_num_seq_match=min_allowed_tar_num_seq_match
					#End Experimental 
					continue  ##??
				
				#Step8 choose sequences
				#Debug HACK, remove
				cut_ratio=0.2
				min_cut_ratio_2_pass=0.05
				#Also detect here if it can't improve anymore in that position
				top_best_possible_sequences=[]
				#while ( (len(top_best_possible_sequences) < min_number_seq_2_pass_phase_2)
				#           and (cut_ratio >= min_cut_ratio_2_pass)):
				if True:
					##print "TEST:", min_number_seq_2_pass_phase_2, len(top_best_possible_sequences),cut_ratio, min_cut_ratio_2_pass
					top_best_possible_sequences=[]
					for seq in best_possible_sequences:
						
						#print "A", seq
						tmp_init_pos=max(0, woirst_position_ndx-3)
						tmp_end_pos=min(woirst_position_ndx+7, len(tmp_aaseq))
						tmp_aaseq_query_old=tmp_aaseq[tmp_init_pos:tmp_end_pos]
						tmp_aaseq_query_new=tmp_aaseq[tmp_init_pos:woirst_position_ndx]
						tmp_aaseq_query_new+=seq
						tmp_aaseq_query_new+=tmp_aaseq[woirst_position_ndx+4:tmp_end_pos]
						
						#for j in range(woirst_position_ndx-3, woirst_position_ndx+6):
						b_is_good=True
						for j in range(len(tmp_aaseq_query_old)-3):
							try:
								#Remember to -1
								tmp_clusters_seq=(clustersDB.assignments_lr_store_data["assignment"][clustersDB.aa_seq_2ndx_dic[tmp_aaseq_query_new[j:j+4]]]-1)
							except KeyError, e:
									b_is_good=False
									print "Error: The sequence ", tmp_aaseq_query_new[j:j+4], "doesn't exist in the database.. skipping"
									break
							#This try is because the input could be random crap
							aaseq_query_old_count=float(0.00001)
							try:
								aaseq_query_old_count=clustersDB.aa_seq_dic[tmp_aaseq_query_old[j:j+4]][0]
							except KeyError, e:
								if b_Debug:
									print "Error: The sequence ", tmp_aaseq_query_old[j:j+4], "doesn't exist in the database.. using pseudo val:",aaseq_query_old_count 
							"""
							#Experimental, don't let the old ratio to be too good:
							old_ratio = (float( num_fragment_sequences_matching_structure[tmp_init_pos+j][1])
										 /aaseq_query_old_count)
							old_ratio = min(max_allowed_tar_seq_match_min_ratio, old_ratio)
							"""

							#####print old_ratio, new_ratio                    
							new_ratio = (float( np.in1d(tmp_clusters_seq, (tmp_assignments_design[tmp_init_pos+j])).sum())
										 /clustersDB.aa_seq_dic[tmp_aaseq_query_new[j:j+4]][0])
							##if (new_ratio < old_ratio):
							if (new_ratio < (1.0/design_quality_entropy_cut_off)):
								###print tmp_init_pos+j, old_ratio,  new_ratio
								b_is_good=False
								break
									
						if(b_is_good):
							top_best_possible_sequences.append(seq)
					cut_ratio-=0.05
				if(len(top_best_possible_sequences) > 0):
					if b_Debug:
						print ("Top Best Compatible SeqA ", cut_ratio+0.05, len(top_best_possible_sequences), top_best_possible_sequences, "Of: ",
							len(best_possible_sequences), tar_num_seq_match ) ##, (target_quality+0.1) )
				
				
				
				b_generate_quality_plot=False
				if b_Debug:
					print ("Found: Top Best Compatible SeqB ", cut_ratio+0.05, len(top_best_possible_sequences), top_best_possible_sequences)
				if (len(top_best_possible_sequences) > 0):
					if b_Debug:
						print "Making mutations and scoring!"
				else:
					avoid_refine_positions.append(woirst_position_ndx)
					if b_Debug:
						print "Warning!!! B not improving sequences found.. skipping this position: ", woirst_position_ndx, "of: ", len(best_possible_sequences)
					#Experimental adaptive expected temperature/ratio change
					if (mc_temperature > 0):
						mc_temperature+=1
						if (mc_temperature > max_mc_temperature):
							mc_temperature=max_mc_temperature
					tar_seq_match_min_ratio-=0.01
					if (tar_seq_match_min_ratio < min_allowed_tar_seq_match_min_ratio):
						tar_seq_match_min_ratio=min_allowed_tar_seq_match_min_ratio
					tar_num_seq_match-=1
					if (tar_num_seq_match < min_allowed_tar_num_seq_match):
						tar_num_seq_match=min_allowed_tar_num_seq_match
					#End Experimental 
					continue
				
				
				connected_pose_mini_design2=connected_pose_mini_desig.clone()
				####min_score=9999    ###self.scorefxn_tal(connected_pose_mini_design2)*0.8
				#print "Initi Score = ", min_score
				
				###Loop optimization:Optimize such position using the NEW MC fragments packer:
				curr_score=self.scorefxn_tal(connected_pose_mini_design2)
			
				print "E Current score vs best:", curr_score, best_total_score
				
				##Check minimization flag
				b_do_min_rotamer_trials=False
				if(design_quality_entropy_cut_off < quality_to_start_rotamer_minimization):
					print "Allowing repack minimization now"
					b_do_min_rotamer_trials=True
					
				b_do_design_packed=False   
				if(design_quality_entropy_cut_off < quality_to_start_structure_relaxation):
					b_do_design_packed=True
					
				#Non symmetric and symmetric packing routines
				best_pose_design=connected_pose_mini_design2.clone() 
				old_score= self.scorefxn_tal(best_pose_design)
				new_score=99999.0
				b_mc_succed=False
				if not b_use_symmetric_packing:
					best_pose_design, new_score, b_mc_succed = self.mc_optimize_pose_position(connected_pose_mini_design2, 
																woirst_position_ndx, 
																top_best_possible_sequences, 
																layer_design_task_bools, 
																sfx_insertion=self.scorefxn_tal_low_farep, # self.scorefxn_tal, self.scorefxn_soft_rep, self.scorefxn_tal_low_farep
																sfx_minimize=self.scorefxn_tal,
																mc_iterations=min(10,(2*len(top_best_possible_sequences))), 
																mc_reinsertions=1,
																min_bb_each=1000,
																mc_temp=mc_temperature,
																pack_radius=8.0,
																b_do_min=b_do_min_rotamer_trials, 
																b_do_final_min=True, 
																b_do_final_relax=False,
																do_allow_design_repacked=False,  ##True, #b_do_design_packed,
																b_Debug=False)
				elif b_use_symmetric_packing:
					best_pose_design, new_score, b_mc_succed = self.mc_optimize_pose_position_symmetric(connected_pose_mini_design2,
															symm_positions_dic,
															woirst_position_ndx, 
															top_best_possible_sequences, 
															layer_design_task_bools, 
															sfx_insertion=self.scorefxn_tal_low_farep, 
															sfx_minimize=self.scorefxn_tal, 
															mc_iterations=min(50,(2*len(top_best_possible_sequences))), 
															mc_reinsertions=1,
															min_bb_each=1000,
															mc_temp=mc_temperature,
															pack_radius=5.0,
															b_do_min=b_do_min_rotamer_trials, 
															b_do_final_min=False,
															b_do_final_relax=False,
															b_Debug=False)
				else:
					print "I shouldn't be here. Error. Terminating"
					assert(0==1)

				new_seq=best_pose_design.sequence()
				print ("Got position/seq: ", (woirst_position_ndx+1), 
					   new_seq[woirst_position_ndx:woirst_position_ndx+4])
			
				#Append current sequence to the dictionary
				alrready_tried_seq[woirst_position_ndx][best_pose_design.sequence()[woirst_position_ndx:woirst_position_ndx+4]]=woirst_position_ndx
				
				#This will avoid ploting in the next round if there is no improvement
				##new_score= self.scorefxn_tal(best_pose_design)
				
				#0.2 means 80% of.
				if ( b_mc_succed ): ##(new_score<old_score) or (abs(new_score-old_score) <= abs(old_score*0.05)) ):
					#Update pose and output
					connected_pose_mini_desig=best_pose_design.clone()
					#Experimental structure recovery
					"""
					current_pose_seq=connected_pose_mini_desig.sequence()
					print "Recovering original structure "
					connected_pose_mini_desig = self.mutate_residue( connected_pose_mini_desig_initial_reference ,(np.arange(len(current_pose_seq))+1) , 
															   current_pose_seq , [], pack_radius = 7.0 , 
															   b_allow_design=False, pack_scorefxn = self.scorefxn_tal )
					"""
					print "Minimizing with constraints!"
					self.minimize_all_with_bb_constraints(connected_pose_mini_desig, 
													sfx=self.scorefxn_tal_low_farep_sstrand, 
													harmonic_constraint_streght=0.1, 
													cart_constraint_weight=1.0)
					#END Experimental recovery  
					
					print "E Updating sequence and all. Energies (new vs old vs best):", new_score, "vs", old_score, "vs", best_total_score, "\n"
					#Step10
					tmp_rmsd_after_packing=self.align_by_bb(connected_pose_mini_desig_best_optimization, connected_pose_mini_desig_initial_reference)
					print "RMSD after packing:", tmp_rmsd_after_packing
					
					#This is the structure that will be ploted in the next iteration
					connected_pose_mini_desig.dump_scored_pdb("%s/%s/merged_test_phase1_designed_opt_c%04d_%04d_p%04d.pdb"%(
																							general_out_dir,
																							stage2_tmp_dir_name,
																							outer_optimization_cycle,
																							inner_optimization_cycle,
																							(woirst_position_ndx+1)), self.scorefxn_tal)
					#Allow to plot in the next iteration
					b_generate_quality_plot=True
					#Save a list of the optimized positions
					#Experimental (raise the bar):
					tar_num_seq_match+=2 #5
					if(tar_num_seq_match > max_allowed_tar_num_seq_match):
						tar_num_seq_match = max_allowed_tar_num_seq_match
					mc_temperature-=10 #20
					if (mc_temperature < min_allowed_mc_temperature):
						mc_temperature=min_allowed_mc_temperature
						
					tar_seq_match_min_ratio+=0.05
					if (tar_seq_match_min_ratio > max_allowed_tar_seq_match_min_ratio):
						tar_seq_match_min_ratio=max_allowed_tar_seq_match_min_ratio
						
					#Save the best Max quality until now
					if (tmp_seq_entropy_ratio.max() < design_quality_entropy_cut_off):
						print "Output best-worst quality design until now: ", tmp_seq_entropy_ratio.max(), design_quality_entropy_cut_off
						#TEST
						
						
						current_pose_seq=connected_pose_mini_desig.sequence()
						connected_pose_mini_desig_best_maxquality = connected_pose_mini_desig.clone()
						"""
						if b_reassing_and_minimize_before_output_struct:
							print "Recovering original structure for output purposes"
							connected_pose_mini_desig_best_maxquality = self.mutate_residue( connected_pose_mini_desig_initial_reference ,(np.arange(len(current_pose_seq))+1) , 
																	   current_pose_seq , [], pack_radius = 7.0 , 
																	   b_allow_design=False, pack_scorefxn = self.scorefxn_tal )
							#Minimize with cst
							print "Minimizing with constraints!"
							self.minimize_all_with_bb_constraints(connected_pose_mini_desig_best_maxquality, 
															sfx=self.scorefxn_tal, 
															harmonic_constraint_streght=0.5, 
															cart_constraint_weight=0.7)
							#Get RMSD
							tmp_rmsd_after_min=self.align_by_bb(connected_pose_mini_desig_best_maxquality, connected_pose_mini_desig_initial_reference)
							print "RMSD after minimization:", tmp_rmsd_after_min
							print "Done recovering"
						"""
						#ENDTEST
						connected_pose_mini_desig_best_maxquality.dump_scored_pdb(
												"%s/%s/merged_test_phase1_designed_c%04d_%04d_p%04d_BestTargetMaxQual_%05.1f.pdb"%(
																							general_out_dir,
																							stage2_tmp_dir_name,
																							outer_optimization_cycle,
																							inner_optimization_cycle,
																							(woirst_position_ndx+1),
																							design_quality_entropy_cut_off), self.scorefxn_tal)
						
					#Save the best seqQual until now (without considering restricted regions)
					avg_seq_entropy=np.average(tmp_seq_entropy_ratio[np.where(tmp_seq_entropy_ratio > -1.0)])
					#Hacked 10000 number, remove
					if ( avg_seq_entropy < best_worst_ratio_value ): ###and (tmp_seq_entropy_ratio.max() < 10000)):
						print "Updating best average design until now: ", best_worst_ratio_value, tmp_seq_entropy_ratio.max()
						num_cycles_without_improvement=0
						best_worst_ratio_value=avg_seq_entropy
						best_optimization_list_of_optimized_positions=list(list_of_optimized_positions)
						##best_worst_ratio_cycle_num=opt_cycle_number-1
						
						#TEST
						current_pose_seq=connected_pose_mini_desig.sequence()
						connected_pose_mini_desig_best_optimization = connected_pose_mini_desig.clone()
						"""
						if b_reassing_and_minimize_before_output_struct:
							print "Recovering original structure for output purposes"
							connected_pose_mini_desig_best_optimization = self.mutate_residue( connected_pose_mini_desig_initial_reference ,(np.arange(len(current_pose_seq))+1) , 
																	   current_pose_seq , [], pack_radius = 7.0 , 
																	   b_allow_design=False, pack_scorefxn = self.scorefxn_tal )
							#Minimize with cst
							print "Minimizing with constraints!"
							self.minimize_all_with_bb_constraints(connected_pose_mini_desig_best_optimization, 
															 sfx=self.scorefxn_tal,
															 harmonic_constraint_streght=0.5, 
															 cart_constraint_weight=0.7)
							#Get RMSD
							tmp_rmsd_after_min=self.align_by_bb(connected_pose_mini_desig_best_optimization, connected_pose_mini_desig_initial_reference)
							print "RMSD after minimization:", tmp_rmsd_after_min
							print "Done recovering"
							#ENDTEST
						"""

						#CommentedForTest
						#connected_pose_mini_desig_best_optimization=connected_pose_mini_desig.clone()
						
						##if (best_worst_ratio_cycle_num >= 0):
						connected_pose_mini_desig_best_optimization.dump_scored_pdb(
												"%s/%s/merged_test_phase1_designed_c%04d_%04d_p%04d_BestAvgSeq%05.1f.pdb"%(
																						general_out_dir,
																						stage2_tmp_dir_name,
																						outer_optimization_cycle,
																						inner_optimization_cycle,
																						(woirst_position_ndx+1),
																						avg_seq_entropy), self.scorefxn_tal)
					else:
						num_cycles_without_improvement+=1
						if ( num_cycles_without_improvement > num_cycles_without_improvement_tol):
							##print "WARNING, I have reached the max number of cycles without improvement in worst position. Ending now"
							print "WARNING, I have reached the max number of cycles without improvement in worst position. Reducing the cut!"
							tar_num_seq_match-=5
							if (tar_num_seq_match < min_allowed_tar_num_seq_match):
								tar_num_seq_match=min_allowed_tar_num_seq_match
							print "New cut: ", tar_num_seq_match
							##num_cycles_without_improvement=0
							#Experimental temperature change:
							if (mc_temperature > 0):
								mc_temperature+=5
								if (mc_temperature > max_mc_temperature):
									mc_temperature=max_mc_temperature
							#Experimental adaptive expected ratio change
							tar_seq_match_min_ratio-=0.05
							if (tar_seq_match_min_ratio < min_allowed_tar_seq_match_min_ratio):
								tar_seq_match_min_ratio=min_allowed_tar_seq_match_min_ratio
							#Reset the failure counter (this is equivalent to an MC minimizer counter)
							num_cycles_without_improvement=0
								
							##connected_pose_mini_desig.dump_scored_pdb(
							##                        "%s/%s/merged_test_phase1_designed_optFinal_Forced.pdb"%(general_out_dir, 
							##																					stage2_tmp_dir_name, self.scorefxn_tal))
							###break
						
					##Out Poses Section###
					#Save the best score until now 
					#ToDo: Change name "old_score" to something less missleading
					if (old_score < best_total_score):
						print "Updating best Energy-wise design until now: ", best_total_score, old_score
						best_total_score=old_score
						#TEST
						
						current_pose_seq=connected_pose_mini_desig.sequence()
						connected_pose_mini_desig_best_energy = connected_pose_mini_desig.clone()
						"""
						if b_reassing_and_minimize_before_output_struct:
							print "Recovering original structure for output purposes"
							connected_pose_mini_desig_best_energy = self.mutate_residue( connected_pose_mini_desig_initial_reference ,(np.arange(len(current_pose_seq))+1) , 
																	   current_pose_seq , [], pack_radius = 7.0 , 
																	   b_allow_design=False, pack_scorefxn = self.scorefxn_tal )
							#Minimize with cst
							print "Minimizing with constraints!"
							self.minimize_all_with_bb_constraints(connected_pose_mini_desig_best_energy, 
															 sfx=self.scorefxn_tal, 
															 harmonic_constraint_streght=0.5, 
															 cart_constraint_weight=0.7)
							#Get RMSD
							tmp_rmsd_after_min=self.align_by_bb(connected_pose_mini_desig_best_energy, connected_pose_mini_desig_initial_reference)
							print "RMSD after minimization:", tmp_rmsd_after_min
							print "Done recovering"
							#ENDTEST
						"""
						
						connected_pose_mini_desig_best_energy.dump_scored_pdb("%s/%s/merged_test_phase1_designed_c%04d_%04d_p%04d_BestEner%05.1f.pdb"%(
																							general_out_dir,
																							stage2_tmp_dir_name,
																							outer_optimization_cycle,
																							inner_optimization_cycle,
																							(woirst_position_ndx+1),
																							best_total_score), 
																							self.scorefxn_tal)
						
				else:
					print "E The improvement is not good, try again later :(!\n", new_score, "vs", old_score
					if (mc_temperature > 0):
						mc_temperature+=1
						if (mc_temperature > max_mc_temperature):
							mc_temperature=max_mc_temperature
					#Experimental adaptive expected temperature/ratio change
					tar_seq_match_min_ratio-=0.05
					if (tar_seq_match_min_ratio < min_allowed_tar_seq_match_min_ratio):
						tar_seq_match_min_ratio=min_allowed_tar_seq_match_min_ratio
					tar_num_seq_match-=5
					if (tar_num_seq_match < min_allowed_tar_num_seq_match):
						tar_num_seq_match=min_allowed_tar_num_seq_match
					#End Experimental 
					continue
				
				
				for ii in range(4):
					list_of_optimized_positions.append(woirst_position_ndx+ii+1)
					if b_use_symmetric_packing:
						for kk in  symm_positions_dic[woirst_position_ndx+ii+1]:
							list_of_optimized_positions.append(kk)
					
				print "DONE cycle", opt_cycle_number, ", out-in: ", outer_optimization_cycle, "-", inner_optimization_cycle
			
			print "End of inner cycle..."
			#if (mc_temperature <= mc_temperature_minimization_start):
			##print best_outer_design_cycle_pose.sequence()
			
			###if(design_quality_entropy_cut_off < quality_to_start_rotamer_minimization):
			if True:
				"""
				#HACK Temporarily disabled
				##First recover the original structure!!! (Avoid drifting)
				print "Recovering original structure and Aplying current sequence"
				current_pose_seq=connected_pose_mini_desig.sequence()
				connected_pose_mini_desig = connected_pose_mini_desig_initial_reference.clone()
				connected_pose_mini_desig = self.mutate_residue( connected_pose_mini_desig ,(np.arange(len(current_pose_seq))+1) , 
														   current_pose_seq , [], pack_radius = 7.0 , 
														   b_allow_design=False, pack_scorefxn = self.scorefxn_tal )
				print "Done recovering"
				

				#Align the non-minimized pose based on the input reference
				connected_pose_mini_desig_before_mini=connected_pose_mini_desig.clone()
				self.align_by_bb(connected_pose_mini_desig_before_mini, connected_pose_mini_desig_initial_reference)
				
				#Minimize with cst
				print "Now minimizing with constraints!"
				self.minimize_all_with_bb_constraints(connected_pose_mini_desig, 
												 sfx=self.scorefxn_tal, 
												 harmonic_constraint_streght=0.5, 
												 cart_constraint_weight=0.7)
				"""
				#Get RMSD
				tmp_rmsd_after_min=self.align_by_bb(connected_pose_mini_desig, connected_pose_mini_desig_initial_reference)
				print "RMSD after minimization:", tmp_rmsd_after_min
				
					
				if (tmp_rmsd_after_min < out_cycle_max_minimization_rmsd):
					print "Minimization looks good!:",  tmp_rmsd_after_min
					best_outer_design_cycle_pose=connected_pose_mini_desig.clone()
					
					if(design_quality_entropy_cut_off <= quality_to_start_structure_relaxation):
						print "MC Quality reached the threshold: ", quality_to_start_structure_relaxation
						print "Now relaxing with csts"
						connected_pose_relax_desig=connected_pose_mini_desig.clone()
						relax_with_constraints_to_init(connected_pose_relax_desig, sfx=self.scorefxn_tal) 
						tmp_rmsd_after_rel=self.align_by_bb(connected_pose_relax_desig, connected_pose_mini_desig_initial_reference)
						if (tmp_rmsd_after_rel < out_cycle_max_relax_rmsd):
							#Save changes if they don't break the pose too much            
							print "Minimization/relaxation looks good! Accepting changes (min/rel: %0.4f / %0.4f)" %(tmp_rmsd_after_min, tmp_rmsd_after_rel)
							##Out Disabled temporarily
							##connected_pose_mini_desig_before_mini.dump_scored_pdb("%s/tmp_optimizations/design_opt_endOfInnerCycle_num%05d_nonmin.pdb"%(
							##                                                                general_out_dir, outer_optimization_cycle), self.scorefxn_tal )
							##connected_pose_mini_desig.dump_scored_pdb("%s/tmp_optimizations/design_opt_endOfInnerCycle_num%05d_min.pdb"%(
							##                                                                general_out_dir, outer_optimization_cycle), self.scorefxn_tal )
							##connected_pose_mini_desig.dump_scored_pdb("%s/tmp_optimizations/design_opt_endOfInnerCycle_num%05d_rel.pdb"%(
							##                                                                general_out_dir, outer_optimization_cycle), self.scorefxn_tal )
							#Update stuff
							best_outer_design_cycle_pose=connected_pose_relax_desig.clone()
							connected_pose_mini_desig=connected_pose_relax_desig.clone()
						
						else: 
							print "The changes break the pose at the Relaxation step: ", tmp_rmsd_after_min, tmp_rmsd_after_rel
							print "Rejecting  this cycle changes and recovering last min best"
							connected_pose_mini_desig=best_outer_design_cycle_pose.clone()
							##Out Disabled temporarily
							##connected_pose_mini_desig.dump_scored_pdb("%s/tmp_optimizations/design_opt_endOfInnerCycle_num%05d_min.pdb"%(
							##                                                                general_out_dir, outer_optimization_cycle), self.scorefxn_tal )
					else:
						None
						##Out Disabled temporarily
						##connected_pose_mini_desig_before_mini.dump_scored_pdb("%s/tmp_optimizations/design_opt_endOfInnerCycle_num%05d_nonmin.pdb"%(
						##                                                                    general_out_dir, outer_optimization_cycle), self.scorefxn_tal )
						##connected_pose_mini_desig.dump_scored_pdb("%s/tmp_optimizations/design_opt_endOfInnerCycle_num%05d_min.pdb"%(
						##                                                                    general_out_dir, outer_optimization_cycle), self.scorefxn_tal )
				else:
					print "RMSD to high, breaking the cycle/optimization and all, sorry"
					break
					##assert(0==1)
					#print "The changes break the pose at the Minimization step: ", tmp_rmsd_after_min
					#print "Rejecting  this cycle changes and recovering last best"
					#connected_pose_mini_desig=best_outer_design_cycle_pose.clone()
					
		 
			###Oscillator allowing eaach other outer cycle to optimizae any position (allows faster convergence)
			if b_allow_alternate_good_and_bad_sequence_regions_design_each_in_cycle:
				if b_allow_change_restricted_positions:    
					b_allow_change_restricted_positions=False
					#elif (design_quality_entropy_cut_off): # < quality_to_start_rotamer_minimization):
					#    b_allow_change_restricted_positions=True
				else:
					b_allow_change_restricted_positions=True
				
			
		#Relax final structure and output
		##print "Executing Fast Relax"
		##fastR.apply(connected_pose_mini_desig)
		#num_cycles_scoring=2
		#tmp_scores=np.zeros(num_cycles_scoring, float)
		#for i in range(num_cycles_scoring):
		final_score = self.scorefxn_tal(connected_pose_mini_desig)
		print "Finished Optimization Cycles. Final score = ", final_score, "Ref:", target_score_design_reference
		#connected_pose_mini_desig.dump_scored_pdb(
		#				"%s/%s/merged_test_phase1_designed_099_optFinal.pdb"%(general_out_dir, 
		#																		stage2_tmp_dir_name, 
		#																		self.scorefxn_tal) )
		#Make unique list of optimized positions
		list_of_optimized_positions=np.sort(np.unique(list_of_optimized_positions))
		best_optimization_list_of_optimized_positions=np.sort(np.unique(best_optimization_list_of_optimized_positions))
			
		#Print something
		print "Finished optimization :), the optimized positions were:", list_of_optimized_positions
		print "I'll see you again soon (9r0b = pp.pp%)"
		return connected_pose_mini_desig, connected_pose_mini_desig.sequence(), tmp_score,  best_optimization_list_of_optimized_positions, 
		
		
	def check_sequence_qual_using_clustered_fragments_stats(self,
											in_pose=None,
											clustersDB=None,
											permissive=False,
											target_max_cluster_distance=0.8,
											out_gen_prefix="",
											b_generate_quality_plot=False,
											b_Debug=False  ):
		
		#DoNoT modify the original input, use a copy"
		#Step1
		connected_pose_mini_desig=in_pose.clone()
		tmp_aaseq=connected_pose_mini_desig.sequence()
		
		#Step2
		#Assign structure to clusters within certain cutoff
		tmp_assignments_design = self.assign_pose_to_clusters_within_cutoff(connected_pose_mini_desig, 
																			1, 
																			connected_pose_mini_desig.total_residue(), 
																			cut_off=target_max_cluster_distance, 
																			clustersDB=clustersDB,
																			permissive=False)
																			
																			
		num_fragment_sequences_match_per_position = np.zeros((len(tmp_aaseq)-3,2), int)
		for i in range(len(tmp_aaseq)-3):
			num_fragment_sequences_match_per_position[i,0]=i
			try:
				num_fragment_sequences_match_per_position[i,1]= clustersDB.aa_seq_dic[tmp_aaseq[i:i+4]][0] 
			except KeyError, e:
				num_fragment_sequences_match_per_position[i,1]=0
				print "Error: The sequence ", tmp_aaseq[i:i+4], "doesn't exist in the database.. skipping"
				
		
		#Step3
		#ToDO: modify to use dictionaries
		#Now figure out if our clusters contain such sequences per fragment and how many of them
		
		num_fragment_sequences_matching_structure = np.zeros((len(tmp_aaseq)-3,2), int)
		for ipos in range(len(tmp_assignments_design)):
			#Remember to -1
			try:
				tmp_clusters_seq=(clustersDB.assignments_lr_store_data["assignment"][clustersDB.aa_seq_2ndx_dic[tmp_aaseq[ipos:ipos+4]]]-1)  
			except KeyError, e:
				print "Error: The sequence ", tmp_aaseq[ipos:ipos+4], "doesn't exist in the database.. skipping"
				tmp_clusters_seq=[]
			num_fragment_sequences_matching_structure[ipos,0]=ipos
			num_fragment_sequences_matching_structure[ipos,1]=np.in1d(tmp_clusters_seq, (tmp_assignments_design[ipos])).sum()
			
		#Now plot total seq vs contained sequences
		if (b_generate_quality_plot):
			fig = plt.figure(figsize=(15,10))
			ax = fig.add_subplot(2,1,1)
			ax.axis([num_fragment_sequences_match_per_position[:,0].min(), 
					 num_fragment_sequences_match_per_position[:,0].max(), 
					 -500, 0.5])
			ax.plot(num_fragment_sequences_match_per_position[:,0], -num_fragment_sequences_match_per_position[:,1], c='blue')   
			ax.plot(num_fragment_sequences_matching_structure[:,0], -num_fragment_sequences_matching_structure[:,1], c='red')
			fig.savefig("%s_plot_diff.png"% (out_gen_prefix) ) 
		
		#Step4
		#Now plot the difference (This is entropy at 1kT)
		small_val_pert=0.00001
		seq_entropy_ratio = (num_fragment_sequences_match_per_position[:,1])/(
											num_fragment_sequences_matching_structure[:,1]+small_val_pert)
		##print seq_entropy_ratio
		##assert 0==1
		#Plot ratio_entropy
		if (b_generate_quality_plot):
			fig = plt.figure(figsize=(15,10))
			ax = fig.add_subplot(2,1,1)
			ax.axis([num_fragment_sequences_match_per_position[:,0].min(), 
					 num_fragment_sequences_match_per_position[:,0].max(), 
					 0.99, 100])
			ax.set_yscale('log')
			ax.plot(num_fragment_sequences_match_per_position[:,0], seq_entropy_ratio, c='blue') 
			fig.savefig("%s_plot_qual.png"% (out_gen_prefix))
																
																
		return [num_fragment_sequences_match_per_position, 
				num_fragment_sequences_matching_structure, 
				seq_entropy_ratio]
				



    

