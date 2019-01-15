#The Protein Idealizer (a.k.a. The De Novo Protein Mimic Designer) by D.A.S.
#Authors: Daniel-Adriano Silva, Enrique Marcos, Javier Castellanos and David Baker.
##Contributions:
### D.-A.S.: Original Idea, Algoritmic and Code Implementation.
### E.M.: Advice on priciples for de novo protein design.
### J.C.: Advice on Automation
### D.B.: Advisor, P.I.

#The Baker Lab. UW, Seattle, USA. 2014-2018"
#Date: Dec/16/2018

print "De Novo Protein Mimic Designer, by D.A.S."

print "\n Sure! There are indeed more than a billion ways to implement the ideas here.\nNevertheless, what is really fundamental (and novel) here is the concept and what it proofs (i.e.'same structure equals same function'\nD.A.S."

print "Citation: \n  Daniel-Adriano Silva*, Shawn Yu*, Umut Y. Ulge*, Jamie B. Spangler, et. al., De novo design of potent and selective mimics of IL-2 and IL-15, Nature, 2019. https://doi.org/10.1038/s41586-018-0830-7"

print "HI D. Let's start!"

import sys
if not "/software/pyrosetta2/latest" in sys.path: 
    print "ADDING"
    sys.path=sys.path+["/software/pyrosetta2/latest"]
print sys.path


print "Loading the libraries for the Protein Idealizer (a.k.a. The Protein Mimetics Designer), by D.A.S."


#Import libs
from os import path 
import os 
import numpy as np
import itertools 
import time 
import pandas 
import pickle 
import hashlib
import random
import collections

from operator import itemgetter

from scipy import constants as scipycons
from scipy.stats import mode as scipy_mode

import matplotlib
if "DISPLAY" not in os.environ:
    print "***Enabling Agg for mapplotlib since there is no DIPLAY env"
    matplotlib.use('Agg') #!!!Only enabled for command-line run
import matplotlib.pyplot as plt



#Import pyRosetta
import pyrosetta

#Other general aliases
##FastRel=pyrosetta.rosetta.protocols.relax.FastRelax
CartMin=pyrosetta.rosetta.core.optimization.CartesianMinimizer
HarmonicFunc=pyrosetta.rosetta.core.scoring.func.HarmonicFunc
CoorCstr=pyrosetta.rosetta.core.scoring.constraints.CoordinateConstraint
##mzr=pyrosetta.rosetta.core.optimization.AtomTreeMinimizer

#import interface_fragment_matching as interface_fragment_matching
import pyRMSD.RMSDCalculator as fast_rmsd_calc
#from interface_fragment_matching.fragment_fitting.rmsd_calc import rmsd_calc

print pyrosetta.__file__


#Definitions
def calculate_disembodied_ss_spans(inpose=None, 
                                   min_allowed_ss_size_dic={}, 
                                   discard_positions=[], #1-ndx
                                   force_h_positions=[], #1-ndx
                                   force_e_positions=[], #1-ndx
                                   b_break_by_chain=True,
                                   keep_chains=[],
                                   b_detect_breaks=True,
                                   propagate_neighbors_NC_ss=False,
                                   detect_E_bulges=False,
                                   NC_threshold_dist=1.7,
                                   bdebug=False):
    
    #break_seq_chars={'L','G','C'}
    allowed_ss_types={'E','H','K'} #min_allowed_ss_size_dic.keys
    gap_code='G'
    keep_intact_code='K'
    tmp_dssp_instance = pyrosetta.rosetta.core.scoring.dssp.Dssp(inpose)
    tmp_dssp_sequence = np.ascontiguousarray(list(tmp_dssp_instance.get_dssp_reduced_IG_as_L_secstruct()),str)
        
    if detect_E_bulges:
        ##print tmp_dssp_sequence
        for ires in np.where(tmp_dssp_sequence == 'E' )[0]:
            if ( (inpose.phi(ires+1) > 0.0) or 
                 (inpose.psi(ires+1) < 0.0) ):
                #print ires, inpose.phi(ires+1) , inpose.psi(ires+1)
                tmp_dssp_sequence[ires]='B'

    for i in discard_positions:
        tmp_dssp_sequence[i-1]='C'
        
    for i in force_h_positions:
        tmp_dssp_sequence[i-1]='H'
    
    for i in force_e_positions:
        tmp_dssp_sequence[i-1]='E'
    
     ##print tmp_dssp_sequence
    if propagate_neighbors_NC_ss:
        if ((tmp_dssp_sequence[1] == 'H') or ( tmp_dssp_sequence[1] == 'E')):
            tmp_dssp_sequence[0]=tmp_dssp_sequence[1]
        elif((tmp_dssp_sequence[2] == 'H') or (tmp_dssp_sequence[2] == 'E')):
            tmp_dssp_sequence[0]=tmp_dssp_sequence[2]
            tmp_dssp_sequence[1]=tmp_dssp_sequence[2]
        if ((tmp_dssp_sequence[-2] == 'H') or (tmp_dssp_sequence[-2] == 'E')):
            tmp_dssp_sequence[-1]=tmp_dssp_sequence[-2]
        elif((tmp_dssp_sequence[-3] == 'H') or (tmp_dssp_sequence[-3] == 'E')):
            tmp_dssp_sequence[-1]=tmp_dssp_sequence[-3]
            tmp_dssp_sequence[-2]=tmp_dssp_sequence[-3]
        for ires in xrange(0,len(tmp_dssp_sequence)-1):
            if (tmp_dssp_sequence[ires]=='L' and tmp_dssp_sequence[ires+1]!='L'):
                tmp_dssp_sequence[ires]= tmp_dssp_sequence[ires+1]
        for ires in xrange(1,len(tmp_dssp_sequence)):
            if (tmp_dssp_sequence[ires-1]=='L' and tmp_dssp_sequence[ires]!='L'):
                tmp_dssp_sequence[ires-1]= tmp_dssp_sequence[ires]
                
    tmp_dssp_sequence_oriProp=np.copy(tmp_dssp_sequence)
                
    if bdebug:
        print "ORIPROP:", tmp_dssp_sequence_oriProp 
        #print "GAPPED :", tmp_dssp_sequence
    
    
    system_size=len(tmp_dssp_sequence_oriProp)
    contiguous_SS_spans=[]
    contiguous_SS_spans_type=[]
    #numpy 'or' operator is different (i.e. '|') to python or, remember always!!!
    curr_ss_type_startNdx=0
    curr_ss_type=tmp_dssp_sequence_oriProp[curr_ss_type_startNdx]
    ##(THIS usually means poor resolution of those residues or a break in the chain, so we will ignore them)
    
    chain_end_ndxs={}
    if b_detect_breaks:
        print "Auto detecting chain breaks"
        gap_residues=check_protein_for_chain_breaks(inpose, 
                                                treshold_dist=NC_threshold_dist)
        for ires in gap_residues:
            chain_end_ndxs[ires]=1    

    if b_break_by_chain:
        print "Also considering chain ends! as TERS"
        inpose_by_chain=inpose.split_by_chain()
        #print inpose_by_chain
        start=1
        for ichain in xrange(1,len(inpose_by_chain)+1):
            end=start+inpose_by_chain[ichain].total_residue()-1
            chain_end_ndxs[end]=1 
            if ichain in keep_chains:
                ##print start,end
                for jres in xrange(start,end+1):
                    tmp_dssp_sequence_oriProp[jres-1]=keep_intact_code
                    tmp_dssp_sequence[jres-1]=keep_intact_code
            start+=inpose_by_chain[ichain].total_residue()
    print chain_end_ndxs
    
    for ires in xrange(system_size): #1,system_size):
        tmp_span=[]
        #Terminate span
        if ( (tmp_dssp_sequence_oriProp[ires] != curr_ss_type) or
             (ires in gap_residues) or
             (ires in chain_end_ndxs)
            ):
            if curr_ss_type in allowed_ss_types:
                #if((ires+1) < len(tmp_dssp_sequence_oriProp)):
                    #if (tmp_dssp_sequence_oriProp[ires]==tmp_dssp_sequence_oriProp[ires+1]):
                    #    continue
                #print ires, tmp_dssp_sequence_oriProp[ires], curr_ss_type
                #print tmp_dssp_sequence
                if ( (ires-curr_ss_type_startNdx)>= 
                      min_allowed_ss_size_dic[curr_ss_type] ):
                    contiguous_SS_spans.append([curr_ss_type_startNdx,
                                                ires-1])
                    contiguous_SS_spans_type.append(curr_ss_type)
            #Update!!!
            curr_ss_type=tmp_dssp_sequence_oriProp[ires]
            curr_ss_type_startNdx=ires
            #if (ires in gap_residues):
            #    tmp_dssp_sequence[[(ires-1),ires]]=gap_code
            #    tmp_dssp_sequence_oriProp[[(ires-1),ires]]=gap_code
        elif ((ires+1)==system_size):
            if curr_ss_type in allowed_ss_types:
                if ( (ires-curr_ss_type_startNdx+1)>= 
                     min_allowed_ss_size_dic[curr_ss_type] ):
                    contiguous_SS_spans.append([curr_ss_type_startNdx,
                                                ires])
                    contiguous_SS_spans_type.append(curr_ss_type)
            
    if bdebug:
        print "GAPPED :", tmp_dssp_sequence  
    
    if (len(contiguous_SS_spans)==0):
        return False, tmp_dssp_sequence_oriProp, np.asarray(contiguous_SS_spans)+1, np.asarray(contiguous_SS_spans_type)
        
    return True, tmp_dssp_sequence_oriProp, np.asarray(contiguous_SS_spans)+1, np.asarray(contiguous_SS_spans_type)


def pose_atoms_to_nparray(inpose=None, 
                          target_atoms=["CA"]):
    pose_size=inpose.total_residue()
    coors = np.zeros( ( (pose_size*len(target_atoms)), 3 ),float)
    #Extract the coordinates
    tmp_index=0
    for resi in range (pose_size):
        for atom in (target_atoms):
            for dim in range(0,3): 
                coors[tmp_index][dim] = inpose.residue(resi+1).xyz(atom)[dim]
            tmp_index+=1
    return coors


def check_protein_for_chain_breaks(inpose, 
                                   left_residue_atom="C", 
                                   right_residue_atom="N", 
                                   treshold_dist=1.40):
    
    ###R1-C--X--N-R2
    broken_positions=[]
    #WARNING: +1 Rosetta numbering
    for i in xrange(inpose.total_residue()-1):
        #WARNING: +1 Rosetta numbering
        if ((inpose.residue(i+1).xyz( left_residue_atom ) - inpose.residue(i+2).xyz( right_residue_atom )).norm() >= treshold_dist):
            #This is zero indexed again
            #broken_positions.append((i)) 
            broken_positions.append((i+1)) 
        
    return np.asarray(broken_positions)

def extract_pose_spans_as_chains(inpose=None, 
                                 spans=[]): #recives ndx-0
    inpose_pdbinfo=inpose.pdb_info()
    result_pose=pyrosetta.pose_from_sequence("")
    isFirst=True
    for ispan in spans:
        size_seq=ispan[1]-ispan[0]+1
        #seq=""
        #for res in xrange(size_seq):
        #seq="A"*size_seq
        if isFirst:
            result_pose=pyrosetta.pose_from_sequence("")
            #result_pose.copy_segment(1, inpose, 1, ispan[0]+1)
            for jres in xrange(ispan[0],ispan[0]+size_seq):
                result_pose.append_residue_by_bond(inpose.residue(jres))
            #result_pose.copy_segment( size_seq, 
            #                         inpose, 
            #                         1, 
            #                         ispan[0]+1)
            isFirst=False
          
        
        else:
            #tmp_pose=pyrosetta.pose_from_sequence("A")
            #tmp_pose.copy_segment(1, inpose, 1, ispan[0]+1)
            tmp_pose=pyrosetta.pose_from_sequence("")
            for jres in xrange(ispan[0]+1,ispan[0]+size_seq+1):
                tmp_pose.append_residue_by_bond(inpose.residue(jres))
                
            result_pose.append_pose_by_jump(tmp_pose,
                                            tmp_pose.total_residue())
            """
            result_pose=pyrosetta.pose_from_sequence(seq)
            result_pose.copy_segment(ori_pose_size, inpose, 1, ispan[0]+1)
            for jres in xrange(ispan[0]+2,ispan[0]+size_seq+1):
                result_pose.append_residue_by_bond(inpose.residue(jres))
            result_pose.append_pose_by_jump(tmp_pose, 1)
            """
        
        print "oping PDBinfor LABELS"
        #Copy all PDBinfo labels
        curr_ndx_shift=1
        for ires in xrange( ispan[0], ispan[1]+1 ): 
            for jlabel in inpose_pdbinfo.get_reslabels(ires):
                result_pose.pdb_info().add_reslabel(curr_ndx_shift+ires-ispan[0], 
                                                        jlabel)
        curr_ndx_shift+=size_seq #alt: tmp_result_pose.total_residue()

        
    return result_pose

#Empty
def zero_center_ordered_arange(target, 
                               increment):
    result_array=[]
    for i in np.arange(0.0,target,increment):
        for j in (1.0,-1.0):
            k=i*j
            if not(k in result_array):
                result_array.append(k)
    
    return result_array

#Add some type of pre-hashing here to speed things up (a-lot)!!
def sample_helix_fitting_parameters(target_pose=None,
                                    avg_phi=-63.0,
                                    avg_psi=-41.0,
                                    avg_omega=180.0,
                                    max_perturbation_angle_phi=5.0,
                                    max_perturbation_angle_psi=5.0,
                                    max_perturbation_angle_omega=0.0,
                                    max_phi_curvature_pert=5.0,
                                    max_psi_curvature_pert=5.0,
                                    perturbation_step_phi=1.0,
                                    perturbation_step_psi=1.0,
                                    perturbation_step_omega=1.0,
                                    perturbation_step_phi_curvature=1.0,
                                    perturbation_step_psi_curvature=1.0,
                                    curv_phi_pitch_aa=4,
                                    curv_psi_pitch_aa=4,
                                    convergence_RMSD=0.5,
                                    fake_aa_for_design='V',
                                    debug=False):
    
    #Time printing
    start_time = time.time()
    
    
    target_h_size=target_pose.total_residue()
    tmp_fake_h_seq=""
    
    #+2 due to the hanging bonds-residues
    for ires in xrange(target_h_size+2):
        tmp_fake_h_seq+=fake_aa_for_design
    tmp_fake_h_pose=pyrosetta.pose_from_sequence(tmp_fake_h_seq)
    
    #Sample the angle perturbations
    tmp_rmsd_results=[]
    #ToDo:
    ##1. Make it faster
    ##2. Fix angle perturbation possible circular bug.
    for l_phi_pitch_pert in zero_center_ordered_arange(max_phi_curvature_pert+perturbation_step_phi_curvature,
                                        perturbation_step_phi_curvature):
        
        for m_psi_pitch_pert in zero_center_ordered_arange(max_psi_curvature_pert+perturbation_step_psi_curvature,
                            perturbation_step_psi_curvature):
            
            for i_phi_angle in zero_center_ordered_arange(max_perturbation_angle_phi+perturbation_step_phi, 
                                       perturbation_step_phi):
                #Circular angle function
                tmp_phi=avg_phi+i_phi_angle
                if (tmp_phi>180.0):
                    while (tmp_phi>180.0):
                        tmp_phi-=360.0
                elif (tmp_phi<-180.0):
                    while (tmp_phi<-180.0):
                        tmp_phi+=360.0
                
                
                for j_psi_angle in zero_center_ordered_arange(max_perturbation_angle_psi+perturbation_step_psi,
                                           perturbation_step_psi):
                    #Circular angle function
                    tmp_psi=avg_psi+j_psi_angle
                    if (tmp_psi>180.0):
                        while (tmp_psi>180.0):
                            tmp_psi-=360.0
                    elif (tmp_psi<-180.0):
                        while (tmp_psi<-180.0):
                            tmp_psi+=360.0
                    for k_omega_angle in zero_center_ordered_arange(max_perturbation_angle_omega+perturbation_step_omega,
                                           perturbation_step_omega):
                        #Circular angle function
                        tmp_omega=avg_omega+k_omega_angle
                        if (tmp_omega>180.0):
                            while (tmp_omega>180.0):
                                tmp_omega-=360.0
                        elif (tmp_omega<-180.0):
                            while (tmp_omega<-180.0):
                                tmp_omega+=360.0
                        
                        #Apply the angles to the pose
                        for ires in xrange(target_h_size):
                            internal_num=ires+2
                            rosetta_num=ires+1
                            if (curv_phi_pitch_aa > 0):
                                if ((rosetta_num%curv_phi_pitch_aa)==0):
                                    tmp_fake_h_pose.set_phi(internal_num,
                                                        tmp_phi+l_phi_pitch_pert)
                                else:
                                    tmp_fake_h_pose.set_phi(internal_num,
                                                            tmp_phi)
                            else:
                                tmp_fake_h_pose.set_phi(internal_num,
                                                        tmp_phi)
                            
                            if (curv_psi_pitch_aa > 0):
                                if ((rosetta_num%curv_psi_pitch_aa)==0):
                                    tmp_fake_h_pose.set_psi(internal_num,
                                                        tmp_psi+m_psi_pitch_pert)
                                else:
                                    tmp_fake_h_pose.set_psi(internal_num,
                                                            tmp_psi)
                            else:
                                    tmp_fake_h_pose.set_psi(internal_num,
                                                            tmp_psi)
                            tmp_fake_h_pose.set_omega(internal_num,
                                                    tmp_omega)
                        
                        #Calculate results
                        rmsd_result=rmsd_atoms_by_ndxs(tmp_fake_h_pose, 
                                                     2, 
                                                     target_h_size+1, 
                                                     target_pose, 
                                                     1, 
                                                     target_h_size,
                                                      atoms=["CA","C","O","N"])
                        
                        tmp_rmsd_results.append(np.array([rmsd_result,
                                                          tmp_phi,
                                                          tmp_psi,
                                                          tmp_omega,
                                                          l_phi_pitch_pert,
                                                          m_psi_pitch_pert,
                                                          curv_phi_pitch_aa,
                                                          curv_psi_pitch_aa]))
                        if debug:
                            print "RMSD:", rmsd_result,l_phi_pitch_pert,m_psi_pitch_pert,i_phi_angle,j_psi_angle,k_omega_angle, tmp_rmsd_results[-1]
                        if ( rmsd_result <= convergence_RMSD ):
                            print "Stopping parameters search since I reached the convergence RMSD parameter: %0.2f " % convergence_RMSD
                            return np.asarray(tmp_rmsd_results)
                       
    print " >Full SS-grid sampling completed. Total elapsed Time=%0.1f s ###"%(time.time() - start_time)        
    return np.asarray(tmp_rmsd_results)
        
    
#Empty
def check_pose_bonds_are_correct(in_pose,
                                 NC_lowBound=1.1,
                                 NC_highBound=1.5):
    for ires in range(in_pose.total_residue()-1):
        this_NCdist=(in_pose.residue(ires+2).xyz("N")-in_pose.residue(ires+1).xyz("C")).norm()
        ##print ires+1,"-",ires+2, this_NCdist
        if ((this_NCdist < NC_lowBound) or (this_NCdist > NC_highBound)):
            return False, this_NCdist
    return True, 0.0


#A copy of the pose minus N- C- res
def return_pose_copy_minus_ends(in_pose,
                                fake_aminoacid_letter_for_design="A"):
    tmpline=""
    for ires in range(in_pose.total_residue()-2):
        tmpline+=fake_aminoacid_letter_for_design
    out_pose_minus_ends=pyrosetta.pose_from_sequence(tmpline)
    ##( (Pose)arg1, (int)size, (Pose)src, (int)begin, (int)src_begin)
    out_pose_minus_ends.copy_segment(in_pose.total_residue()-2,
                                                            in_pose,
                                                            1,
                                                            2)
    return out_pose_minus_ends

def rmsd_atoms_by_ndxs(pose1, 
                     init_res1, 
                     end_res1, 
                     pose2, 
                     init_res2, 
                     end_res2,
                     atoms=["CA","C","O","N"]):
        numRes=(end_res1-init_res1+1)
        coorsAB=np.zeros((2,(len(atoms)*numRes),3), float)
        #coorB=np.zeros(((len(atoms)*numRes),3), float)

        counter=0
        for res in range (init_res1, (end_res1+1)):
            for atom in atoms:
                for dim in range(0,3): 
                    coorsAB[0,counter,dim]=(pose1.residue(res).xyz(atom)[dim])
                counter+=1

        counter=0
        for res in range (init_res2, (end_res2+1)):
            for atom in atoms:
                for dim in range(0,3): 
                    coorsAB[1,counter,dim]=(pose2.residue(res).xyz(atom)[dim])
                counter+=1

        #print coorsAB
        rmsd_calculator=fast_rmsd_calc.RMSDCalculator("QCP_OMP_CALCULATOR",coorsAB)
        #Calculate the RMSD
        #rmsdVal = rmsd_2_np_arrays(coorB, coorA)

        return rmsd_calculator.pairwise(0,1)
    

def rmsd_and_rosettaRotations_2_np_arrays(crds1, 
                                        crds2):
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
            rMtx_xyzM=pyrosetta.rosetta.numeric.xyzMatrix_double_t()
            rMtx_xyzM.xx=rMtx[0,0]
            rMtx_xyzM.xy=rMtx[0,1]
            rMtx_xyzM.xz=rMtx[0,2]
            rMtx_xyzM.yx=rMtx[1,0]
            rMtx_xyzM.yy=rMtx[1,1]
            rMtx_xyzM.yz=rMtx[1,2]
            rMtx_xyzM.zx=rMtx[2,0]
            rMtx_xyzM.zy=rMtx[2,1]
            rMtx_xyzM.zz=rMtx[2,2]
            tVec_xyzV=pyrosetta.rosetta.numeric.xyzVector_double_t()
            tVec_xyzV.x=tVec[0]
            tVec_xyzV.y=tVec[1]
            tVec_xyzV.z=tVec[2]
            return np.sqrt(rmsd_sq), rMtx_xyzM, tVec_xyzV
        
        
def rmsd_and_rotations_2_np_arrays(crds1, 
                                        crds2):
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
            #tVec=COM1-(np.dot(COM2, np.linalg.inv(rMtx)))
            ###Are you kidding me??? Is this the correct way to build arrays inside 
            #of Rosetta core? No indexing?? Who was the creator of this???
            """rMtx_xyzM=pyrosetta.rosetta.numeric.xyzMatrix_double_t()
            rMtx_xyzM.xx(rMtx[0,0])
            rMtx_xyzM.xy(rMtx[0,1])
            rMtx_xyzM.xz(rMtx[0,2])
            rMtx_xyzM.yx(rMtx[1,0])
            rMtx_xyzM.yy(rMtx[1,1])
            rMtx_xyzM.yz(rMtx[1,2])
            rMtx_xyzM.zx(rMtx[2,0])
            rMtx_xyzM.zy(rMtx[2,1])
            rMtx_xyzM.zz(rMtx[2,2])
            tVec_xyzV=pyrosetta.rosetta.numeric.xyzVector_double_t()
            tVec_xyzV.x=tVec[0]
            tVec_xyzV.y=tVec[1]
            tVec_xyzV.z=tVec[2]"""
            return np.sqrt(rmsd_sq), rMtx, COM1, COM2 #rMtx_xyzM, tVec_xyzV
        
        
def align_atoms_by_ndxs(pose1, 
                         init_res1, 
                         end_res1, 
                         pose2, 
                         init_res2, 
                         end_res2,
                         atoms=["CA","C","O","N"]):
        numRes=(end_res1-init_res1+1)
        coorA=np.zeros(((len(atoms)*numRes),3), float)
        coorB=np.zeros(((len(atoms)*numRes),3), float)

        counter=0
        for res in range (init_res1, (end_res1+1)):
            for atom in atoms:
                for dim in range(0,3): 
                    coorA[counter,dim]=(pose1.residue(res).xyz(atom)[dim])
                counter+=1

        counter=0
        for res in range (init_res2, (end_res2+1)):
            for atom in atoms:
                for dim in range(0,3): 
                    coorB[counter,dim]=(pose2.residue(res).xyz(atom)[dim])
                counter+=1

        #Calculate the RMSD
        rmsdVal, rMtx, tVec = rmsd_and_rosettaRotations_2_np_arrays(coorB, coorA)
        pose1.apply_transform_Rx_plus_v(rMtx, tVec)

        return rmsdVal

dlooper_alg=None
class dlooper():
    ##########
    def parse_dadriano_clustered_fragment_table(self,
                                                table_name,
                                                atom_types_arr=["CA", "C", "O", "N"]):

        pert_pseudocount=0.000001

        def get_fragment_entry_dtype(num_fragment_atoms):
                frag_size=num_fragment_atoms/len(atom_types_arr)
                return np.dtype([("assignment", int), 
                                        ("size", int),
                                        ("avg_distance", float), 
                                        ("threshold_distance", float), 
                                        ("coordinates", float, (num_fragment_atoms, 3)), 
                                        ("bb", float, (frag_size, 3)), 
                                        ("ss", float, (frag_size, 3)),
                                        ("aa", float, (frag_size, 20)), 
                                        ("aa_w", float, (frag_size, 20)) ])

        fragment_table = pandas.read_csv(table_name, index_col=False,  delim_whitespace=True)
        num_coordinates = max([int(c.split("_")[1]) for c in fragment_table.columns if c.startswith("X")])
        fragment_entries = np.zeros(len(fragment_table), get_fragment_entry_dtype(num_coordinates))
        fragment_entries["size"]=fragment_table.NumCon
        fragment_entries["assignment"] = fragment_table.Clust
        fragment_entries["avg_distance"] = fragment_table.AvgRad
        fragment_entries["threshold_distance"] = fragment_table.MaxRad

        for i in xrange(1, num_coordinates + 1):
                fragment_entries["coordinates"][:,i - 1,0] = fragment_table["X_%s" % i]
                fragment_entries["coordinates"][:,i - 1,1] = fragment_table["Y_%s" % i]
                fragment_entries["coordinates"][:,i - 1,2] = fragment_table["Z_%s" % i]

        fragment_aa_size=(num_coordinates/len(atom_types_arr))
        """aaA_1 aaR_1 aaN_1 aaD_1 aaC_1 aaE_1 aaQ_1 aaG_1 aaH_1 aaI_1 aaL_1 aaK_1 aaM_1 aaF_1 aaP_1 aaS_1 aaT_1 aaW_1 aaY_1 aaV_1"""
        for i in xrange(1, (fragment_aa_size+1) ):
                #BB angles position (center only)
                fragment_entries["bb"][:,i - 1,0 ] = fragment_table["Phi_%s" % i]
                fragment_entries["bb"][:,i - 1,1 ] = fragment_table["Psi_%s" % i]
                fragment_entries["bb"][:,i - 1,2 ] = fragment_table["Ome_%s" % i]
                #SS prob/position
                fragment_entries["ss"][:,i - 1,0 ] = fragment_table["ssH_%s" % i] + fragment_table["ssG_%s" % i] + fragment_table["ssI_%s" % i]
                fragment_entries["ss"][:,i - 1,1 ] = fragment_table["ssB_%s" % i] + fragment_table["ssE_%s" % i] 
                fragment_entries["ss"][:,i - 1,2 ] = fragment_table["ssT_%s" % i] + fragment_table["ssS_%s" % i] + fragment_table["ssL_%s" % i]
                #AA prob/position
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

        return fragment_entries,fragment_aa_size #fragment_spec, fragment_spec.fragment_data_to_location_and_atom(fragment_entries), fragment_aa_size

    def parse_dadriano_clustered_fragment_assignments(self,
                                                        table_name,
                                                        frag_size):
        def get_assignment_entry_dtype():
                return np.dtype([("name", str, 200), 
                                 ("span", int, 2) , 
                                 ("assignment", int), 
                                 ("distance", float), 
                                 ("aaSequence", str, frag_size)])
        assignment_table = pandas.read_csv(table_name, index_col=False,  sep='\t|-| ')
        assignment_entries = np.zeros(len(assignment_table), get_assignment_entry_dtype())
        assignment_entries["name"] = assignment_table.Structure
        assignment_entries["span"][:,0] = assignment_table.SpanI
        assignment_entries["span"][:,1] = assignment_table.SpanE
        assignment_entries["assignment"] = assignment_table.Assignment
        assignment_entries["distance"] = assignment_table.Distance
        assignment_entries["aaSequence"] = assignment_table.aaSequence
        return assignment_entries

    ##########
    
    def convert_ss_to_ints(self,
                                    dssp_seq=[]):
        ss_as_ints=np.zeros(len(dssp_seq), int)
        ss_as_ints.fill(-1)
        ss_as_ints[np.where(dssp_seq=='H')[0]]=0
        ss_as_ints[np.where(dssp_seq=='E')[0]]=1
        ss_as_ints[np.where(dssp_seq=='L')[0]]=2
        if -1 in ss_as_ints:
                print "Something is wrong in your ss_def, probably non-standar codes (standard: H,S,L) are present: "
                print dssp_seq
                assert(0==1)
        return ss_as_ints
    ###
    
    
    
    def rmsd_and_rosettarotations_2_np_arrays(self,
                                        crds1, 
                                        crds2):
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
            rMtx_xyzM=pyrosetta.rosetta.numeric.xyzMatrix_double_t()
            rMtx_xyzM.xx(rMtx[0,0])
            rMtx_xyzM.xy(rMtx[0,1])
            rMtx_xyzM.xz(rMtx[0,2])
            rMtx_xyzM.yx(rMtx[1,0])
            rMtx_xyzM.yy(rMtx[1,1])
            rMtx_xyzM.yz(rMtx[1,2])
            rMtx_xyzM.zx(rMtx[2,0])
            rMtx_xyzM.zy(rMtx[2,1])
            rMtx_xyzM.zz(rMtx[2,2])
            tVec_xyzV=pyrosetta.rosetta.numeric.xyzVector_double_t()
            tVec_xyzV.x=tVec[0]
            tVec_xyzV.y=tVec[1]
            tVec_xyzV.z=tVec[2]
            return np.sqrt(rmsd_sq), rMtx_xyzM, tVec_xyzV
        
    
    
    
    
    ###
    ##def: connects 3 poses as pose1--pose2--pose3. It removes the last residue from pose1, 
    def connect_3_poses_byNDX_and_minize_w_constraints(self,
                                                                        pose1=None, 
                                                                        pose2=None, 
                                                                        pose3=None, 
                                                                        ndx_p1b=[], #In rosetta +1 numbering scheme
                                                                        ndx_p2a=[], #In rosetta +1 numbering scheme
                                                                        ndx_p2b=[], #In rosetta +1 numbering scheme
                                                                        ndx_p3a=[], #In rosetta +1 numbering scheme
                                                                        context_pose=None,
                                                                        sfx=None,
                                                                        harmonic_constraint_streght=0.5, 
                                                                        cart_constraint_weight_ramp=[0.5],
                                                                        aminoacid_letter_for_design="V"):
        
        ##Create the score function
        #scorefxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function("soft_rep")
        #scorefxn.set_weight(pyrosetta.rosetta.core.scoring.cart_bonded , 0.7)
        #scorefxn.set_weight(pyrosetta.rosetta.core.scoring.coordinate_constraint , cart_constraint_weight)
        
        #ToDo, check this weight
        sfx.set_weight(pyrosetta.rosetta.core.scoring.cart_bonded , 0.7)
        
        
        #scorefxn_remodel_cen = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function("remodel_cen")
        
        ##Create a tmp pose that will hold the result
        #target_pose_size = p_ss0.total_residue() +  p_ss1.total_residue() + p_loop0_1.total_residue()-4
        ##No need of this beacuse we use from 2nd residue
        #if pose2.residue(1).has_variant_type(pyrosetta.rosetta.core.chemical.VariantType.LOWER_TERMINUS_VARIANT):
        #    pyrosetta.rosetta.core.pose.remove_variant_type_from_pose_residue(pose2, pyrosetta.rosetta.core.chemical.VariantType.LOWER_TERMINUS_VARIANT, 1)
        #if pose2.residue(p_loop0_1.total_residue()).has_variant_type(pyrosetta.rosetta.core.chemical.VariantType.UPPER_TERMINUS_VARIANT):
        #    pyrosetta.rosetta.core.pose.remove_variant_type_from_pose_residue(pose2, pyrosetta.rosetta.core.chemical.VariantType.UPPER_TERMINUS_VARIANT, pose2.total_residue())
        tmp_result_pose = pyrosetta.pose_from_sequence(self.fake_aminoacid_letter_for_design)
        
        #This is the first chunk
        tmp_result_pose.copy_segment(1, pose1, 1, 1)
        #Copy all PDBinfo labels
        pose1_pdbinfo=pose1.pdb_info()
        for jlabel in pose1_pdbinfo.get_reslabels(1):
                    tmp_result_pose.pdb_info().add_reslabel(tmp_result_pose.total_residue(), jlabel)
        #Now from the 2nd to the end!
        for i in range( 2, ndx_p1b+1 ): #range( pose1.total_residue()-1):
                tmp_resi=pose1.residue(i)
                tmp_result_pose.append_residue_by_bond( tmp_resi, 
                                                       False);
                #Copy all PDBinfo labels
                for jlabel in pose1_pdbinfo.get_reslabels(i):
                    tmp_result_pose.pdb_info().add_reslabel(tmp_result_pose.total_residue(), jlabel)
        #Remove pyrosetta.rosetta.core.chemical.VariantType.UPPER_TERMINUS_VARIANT
        ##if tmp_result_pose.residue(tmp_result_pose.total_residue()).has_variant_type(pyrosetta.rosetta.core.chemical.VariantType.UPPER_TERMINUS_VARIANT):
        ##        pyrosetta.rosetta.core.pose.remove_variant_type_from_pose_residue(tmp_result_pose, pyrosetta.rosetta.core.chemical.VariantType.UPPER_TERMINUS_VARIANT, 
        ##                                                                tmp_result_pose.total_residue())
        #This is the loop/connection in the middle, copy all but not EP
        pose2_pdbinfo=pose2.pdb_info()
        #Do not copy first and last since they are the overlaps (that might erase hotspots from the target-SS
        for i in range( ndx_p2a+1, ndx_p2b ):
            tmp_resi=pose2.residue(i)
            tmp_result_pose.append_residue_by_bond( tmp_resi, 
                                                   False);
            #Copy all PDBinfo labels
            for jlabel in pose2_pdbinfo.get_reslabels(i):
                tmp_result_pose.pdb_info().add_reslabel(tmp_result_pose.total_residue(), jlabel)
        #Remove pyrosetta.rosetta.core.chemical.VariantType.UPPER_TERMINUS_VARIANT again!
        if tmp_result_pose.residue(tmp_result_pose.total_residue()).has_variant_type(pyrosetta.rosetta.core.chemical.VariantType.UPPER_TERMINUS_VARIANT):
                pyrosetta.rosetta.core.pose.remove_variant_type_from_pose_residue(tmp_result_pose, pyrosetta.rosetta.core.chemical.VariantType.UPPER_TERMINUS_VARIANT, 
                                                                        tmp_result_pose.total_residue())
        #This is now the last chunk
        pose3_pdbinfo=pose3.pdb_info()
        for i in range( ndx_p3a, pose3.total_residue()+1 ): #range( 1, pose3.total_residue() 
            tmp_resi=pose3.residue(i)
            tmp_result_pose.append_residue_by_bond( tmp_resi, 
                                                   False);
            #Copy all PDBinfo labels
            for jlabel in pose3_pdbinfo.get_reslabels(i):
                tmp_result_pose.pdb_info().add_reslabel(tmp_result_pose.total_residue(), jlabel)
    
        #make a copy to store the final result at the ned
        pose_for_return_result=tmp_result_pose.clone()
        ori_pose_size=tmp_result_pose.total_residue()
        
        if (context_pose.total_residue() > 0):
                tmp_result_pose.append_pose_by_jump(context_pose, 1)
        
        #Loop indexes
        #Note: If ndx_p2a>ndx_p2b we will be in problems! Can fix with an assertion: assert(ndx_p2a<=ndx_p2b)
        loop_start=ndx_p1b #pose1.total_residue()
        loop_end=ndx_p1b+(ndx_p2b-ndx_p2a)   #pose1.total_residue() + pose2.total_residue()-2
        
        #Constraints
        mm_local = pyrosetta.rosetta.core.kinematics.MoveMap()
        #Minimize -1,+1a.a. of the loop's endpoints
        for i in range(loop_start-1, loop_end+2):
                mm_local.set_bb(i, True)
                #Harmonic constraints
                for j in range( 4 ):
                    tmp_result_pose.add_constraint( CoorCstr(pyrosetta.rosetta.core.id.AtomID(j+1, i), 
                                pyrosetta.rosetta.core.id.AtomID(1,1), 
                                tmp_result_pose.residue(i).atom(j+1).xyz(),
                                HarmonicFunc(0.0, harmonic_constraint_streght)) ) 
        
        options_minilbfgs = pyrosetta.rosetta.core.optimization.MinimizerOptions( "lbfgs_armijo_nonmonotone", 
                                                                                            0.001, 
                                                                                            True, 
                                                                                            False, 
                                                                                            False)
        
        #Minimization loop
        minimizer = CartMin()
        for iweight in cart_constraint_weight_ramp:
                sfx.set_weight(pyrosetta.rosetta.core.scoring.coordinate_constraint , iweight)
                minimizer.run( tmp_result_pose, mm_local, sfx, options_minilbfgs )
        
        #Remove constraints from the sfx
        sfx.set_weight(pyrosetta.rosetta.core.scoring.cart_bonded , 0.0)
        sfx.set_weight(pyrosetta.rosetta.core.scoring.coordinate_constraint , 0.0)
        #Remove constraints form the pose:
        tmp_result_pose.remove_constraints()
        
        #finaly copy back just the fragment without context
        pose_for_return_result.copy_segment(ori_pose_size, tmp_result_pose, 1, 1)
    
        return pose_for_return_result, sfx(tmp_result_pose)
    ###
    
    def check_span_for_ideal_fragments(self,
                                        test_pose, 
                                        start, 
                                        end,
                                        rmsd_pert=1.1, 
                                        max_rmsd_limit=1.1,
                                        spans_min_population_cut=30,
                                        default_allowed_blind_raddi=0.5,
                                        cluster_fragment_len=1,
                                        cluster_centers_coordinates=None,
                                        cluster_centers_store_data=None,
                                        rmsd_calculator=None):
        #D assert(start>=1)
        #D assert(end<=test_pose.total_residue())
        ##minndx=-1
        ##
        #ToDo: check this!
        max_population=0 #_sum=0
        rmsd_sum=0
        ##max_rmsd_limit_sum=0
        num_fragments=0
        
        end_minus_fragment=end-cluster_fragment_len+1
        if (end_minus_fragment < start):
            print "***OHOH!!! This shouldn't happen, I can't asses this SS-LL-SS connection because it is too small."
            print "*** Therefore I'll have to skip it for now as a bad result, but it might be actually good! "
            print "*** Please fix me ASAP by adding variable size fragment alignment (or something)"
            return False, [9999.99, 
                          0.0,
                          -1]
        
        for testres in xrange(start , (end_minus_fragment+1)):
                #Compute the span for a Nmer
                testresndx=[testres,(testres+cluster_fragment_len-1)]
                

                #print "DASMDASM",testresndx, self.pose_res_to_array(test_pose, testresndx, atoms=["CA","C","O","N"]), cluster_centers_coordinates.shape                
                tmp_rmsdcalculator=rmsd_calculator.RMSDCalculator( "QCP_OMP_CALCULATOR",np.insert(cluster_centers_coordinates, 
                                                                                                  0, 
                                                                                                  self.pose_res_to_array(test_pose, testresndx, atoms=["CA","C","O","N"]), 
                                                                                                  axis=0) 
                                                                 )
                rmsd_result=tmp_rmsdcalculator.oneVsTheOthers(0)
                
                #Calculate the RMSD
                ##print test_pose.total_residue(), testresndx
                #rmsd_result=rmsd_calculator.get_broadcast_coordinate_rmsd(
                #                np.expand_dims( self.pose_res_to_array(test_pose, testresndx, atoms=["CA","C","O","N"]), 0), 
                #                cluster_centers_coordinates.copy() 
                #                ).ravel()
                
                minndx=rmsd_result.argmin()
                #By deafult allow this to pass
                #ToDO: double check this function
                rmsd_minndx_upperlim=cluster_centers_store_data["threshold_distance"][minndx]*rmsd_pert
                rmsd_target_val=max_rmsd_limit #max(default_allowed_blind_raddi, rmsd_minndx_upperlim)
                if( ( rmsd_result[minndx] <= rmsd_target_val ) ):
                    best_ndxs=np.where(rmsd_result <= rmsd_target_val)[0]
                    tmp_max_population=np.sum(cluster_centers_store_data["size"][best_ndxs]) #max(cluster_centers_store_data["size"][best_ndxs])
                    if (tmp_max_population < spans_min_population_cut):
                        return False, [rmsd_result[minndx],
                                       rmsd_target_val,
                                       tmp_max_population]
                    else: #This is a good result
                        rmsd_sum+=rmsd_result[minndx]
                        if (tmp_max_population>max_population):
                            max_population=tmp_max_population
                        ##max_rmsd_limit_sum+=max_rmsd_limit #This for WT??F Daniel???
                        num_fragments+=1
                else: #if( (rmsd_result[minndx] > rmsd_target_val) or (rmsd_result[minndx] > max_rmsd_limit) ):
                    return False, [rmsd_result[minndx], 
                                   rmsd_target_val,
                                   0]
                    
        return True, [(rmsd_sum*1.0/num_fragments), 
                      max_rmsd_limit,
                      max_population ]
    
    
    
    ###
    def pose_res_to_array(self,
                                    pose, 
                                    ndxs, 
                                    atoms=["CA","C","O","N"]):
        coors = np.zeros(( (len(atoms)*(ndxs[1]-ndxs[0]+1)),3),float)
        #Extract the coordinates
        tmp_index=0
        for resi in range (ndxs[0], ndxs[1]+1):
                for atom in (atoms):
                    for dim in range(0,3): 
                        coors[tmp_index][dim] = pose.residue(resi).xyz(atom)[dim]
                    tmp_index+=1
        return coors
    
    ####
    #Clash checking by hash
    def check_pm_clashes_hash_vs_pose(self,
                                    context_pose_hash,
                                    inpose,
                                    pick_mode=pyrosetta.rosetta.core.pose.PoseCoordPickMode_N_CA_C_CB):

        target_balls = pyrosetta.rosetta.utility.vector1_numeric_geometry_hashing_Ball()
        pyrosetta.rosetta.core.pose.xyzStripeHashPose.extract_pose_balls(inpose,
                                                            target_balls,
                                                            pick_mode)
        target_clash_map = pyrosetta.rosetta.std.map_unsigned_long_unsigned_long()
        context_pose_hash.clash_check_residue_pairs(target_balls,
                                                    target_clash_map)

        return bool(target_clash_map), (np.array(list(target_clash_map.items()))-1)
    
    
    
    ###
    
    def calculate_ss_assignments_by_dssp(self,
                                         in_pose_by_chain=None,
                                         b_flatten_ss_assignments=False ):

        #Get the SS assignments
        assert (len(in_pose_by_chain) >=1)
        if self.b_debug_print:
            print "Detecting existing SS in the input: "
        in_pose_for_loops=in_pose_by_chain[0].clone()
        for indx in xrange(1, len(in_pose_by_chain)):
            in_pose_for_loops.append_pose_by_jump(in_pose_by_chain[indx].clone(), 1)
        tmp_dssp_instance=pyrosetta.rosetta.core.scoring.dssp.Dssp(in_pose_for_loops)
        tmp_dssp_sequence=np.ascontiguousarray(list(tmp_dssp_instance.get_dssp_reduced_IG_as_L_secstruct()))
        ##print tmp_dssp_sequence
        curr_ndx=0
        count_chains=0
        pose_by_chain_ss=[]
        for ipose in in_pose_by_chain:
                #tmp_dssp_instance = pyrosetta.rosetta.core.scoring.dssp.Dssp(ipose)
                #tmp_dssp_sequence = np.ascontiguousarray(list(tmp_dssp_instance.get_dssp_reduced_IG_as_L_secstruct()))
                #Use clusters to get SS?    
                tmp_sub_dssp_seq=tmp_dssp_sequence[curr_ndx:(curr_ndx+ipose.total_residue())]
                curr_ndx+=ipose.total_residue()
                #Now Propagate extremes
                tmp_sub_dssp_seq=self.convert_ss_to_ints(dssp_seq=tmp_sub_dssp_seq)  
                #HACKS TO ASSIGN THE SS "properly"
                tmp_sub_dssp_seq[0]= tmp_sub_dssp_seq[1]
                tmp_sub_dssp_seq[-1]= tmp_sub_dssp_seq[-2]
                this_ss_mode = int(scipy_mode(tmp_sub_dssp_seq)[0])
                #Assign the SS mode as the whole (hacky).
                #ToDo: writte something better
                if b_flatten_ss_assignments:
                    tmp_sub_dssp_seq[np.where(tmp_sub_dssp_seq != this_ss_mode)]=this_ss_mode
                    
                if (((tmp_sub_dssp_seq[0]) == 2 and (count_chains>=1)) or 
                     (tmp_sub_dssp_seq[-1] == 2) and (count_chains<len(in_pose_by_chain))):
                    if ( (not (0 in tmp_sub_dssp_seq)) and b_flatten_ss_assignments ):
                        print "WE HAVE A PROBLEM HUSTON. your SS", count_chains+1,"is all a loop:", tmp_sub_dssp_seq
                        print "I'll try to fix it automatically by converting this to a strand (Don't be surprised if this mess your model!)"
                        tmp_sub_dssp_seq.fill(1)
                        print "Now SS-assignment is:", tmp_sub_dssp_seq
                    elif( (0 in tmp_sub_dssp_seq) and b_flatten_ss_assignments ):
                        print "WE HAVE A PROBLEM HUSTON. your SS", count_chains+1,"is all a loop and has few helix fragments:", tmp_sub_dssp_seq
                        print "I can't fix this, I'll countinue but this might be DOOMED"
                    elif(not b_flatten_ss_assignments):
                        print "You decided to trust DSSP assignments, this might be unwise. Anyway I'll go with this assignment:", tmp_sub_dssp_seq
                        print "Dear DASM,\n PleasE consider fixing Me by using a more flexible SS definition, like a probabilistic one!!! \n Regards,\n Your D'Looper Hacky Code"
                        #assert(0==1)  
                pose_by_chain_ss.append(tmp_sub_dssp_seq[:])
                count_chains+=1
        if self.b_debug_print:
            print "After some magic now it loos like (0=H, 1=E, 2=L):"
            for issndx in xrange(len(pose_by_chain_ss)):
                    print "Chain %d, SS:"%(issndx+1), pose_by_chain_ss[issndx]
        return pose_by_chain_ss
        #
        #
        #END SS
    
    def find_closure_loops(self,
                                in_pose_for_loops=None,
                                pose_by_chain_ss=[],
                                allowVariateSSsizeBy=0,
                                only_as_input_conectivity_order_mode=True,
                                b_use_global_context=False,
                                global_context_pose_hash=None,
                                num_global_context_clashes_initial=0,
                                b_minimize_loops=True,
                                max_num_loop_solutions=99999):
        #Grow this by one becouse 1 == none ;)
        allowVariateSSsizeBy+=1
        
        #Check missleading options???
        #if((not b_flatten_ss_assignments) and
        #   (not b_generate_loop_loop_closures)):
        #    print "!!!Warning!!! You settled b_flatten_ss_assignments and b_generate_loop_loop_closures to FALSE. This might be a MISTAKE. You have been WARNED"
        
        #Check if using global context
        if (b_use_global_context):
            assert (global_context_pose_hash != None)
        
        #By chain and Number of chains???
        in_pose_by_chain=in_pose_for_loops.split_by_chain()
        in_num_chains=len(in_pose_by_chain)
        
        if self.b_debug_print:
            print "Input pose to reloop has number of chains: ", in_num_chains
        if in_num_chains < 2:
            print "Error. I have to STOP. The input should have at least two chains"
            assert (0==1)
        
        #Hash the ss context minus ends in an array
        context_pose_hash_array=[]
        #keep the nparray of the atoms used
        in_pose_by_chain_bb_xyz=[]
        tmp_count=0
        for ichainpose in in_pose_by_chain:
                tmp_count+=1
                #testname="%s/test_%d.pdb"%(out_path,tmp_count)
                #ichainpose.dump_pdb(testname)
                tmpline=""
                for ires in range(ichainpose.total_residue()-2):
                    tmpline+=self.fake_aminoacid_letter_for_design
                tmp_pose_minus_ends=pyrosetta.pose_from_sequence(tmpline)
                ##( (Pose)arg1, (int)size, (Pose)src, (int)begin, (int)src_begin)
                tmp_pose_minus_ends.copy_segment(ichainpose.total_residue()-2,
                                                                        ichainpose,
                                                                        1,
                                                                        2)
                in_pose_by_chain_bb_xyz.append(pose_atoms_to_nparray(ichainpose, ["CA","C","O","N"]))
                context_pose_hash_array.append(pyrosetta.rosetta.core.pose.xyzStripeHashPose(
                                                        tmp_pose_minus_ends,
                                                        pyrosetta.rosetta.core.pose.PoseCoordPickMode_CA, #Maybe N_CA_C 
                                                        radius=2.0))
        
        
        #Pre-build a distance matrix assesment of possible reconnections
        closure_by_dist_matrix=np.zeros((in_num_chains, 
                                                    in_num_chains, 
                                                    allowVariateSSsizeBy, 
                                                    allowVariateSSsizeBy), float)
        closure_by_dist_matrix_allowed_variations=[]
        closure_by_dist_matrix.fill(9999.9)
        
        ##print closure_by_dist_matrix
        
        i_labels=[]
        for i in range(len(in_pose_by_chain)):
                for iChainDelta in range(allowVariateSSsizeBy):
                    this_indx=(in_pose_by_chain[i+1].total_residue()-iChainDelta)
                    for j in range ( len(in_pose_by_chain)):
                        if (i != j):
                                for jChainDelta in range(allowVariateSSsizeBy):
                                    this_jndx=(1+jChainDelta)
                                    #Check that we don't get SS too small 
                                    ##ToDo1: compare to a dictionary of SSsizes??
                                    ##ToDo2: break at hotspots??
                                    
                                    tmp_val_dist_ij=( in_pose_by_chain[i+1].residue(this_indx).xyz( "CA" ) - 
                                                            in_pose_by_chain[j+1].residue(this_jndx).xyz( "CA" ) ).norm()
                                    #print("Debug", this_indx, this_jndx,tmp_val_dist_ij)
                                    
                                    if (tmp_val_dist_ij < closure_by_dist_matrix[i][j][iChainDelta][jChainDelta]):
                                        closure_by_dist_matrix[i][j][iChainDelta][jChainDelta]=tmp_val_dist_ij
        
        #The main loop finding loop :)
        #Find possible loop closures for all different combinations desired
        #Create an array to hold the results
        closure_results_dic={}
        for i in range(len(in_pose_by_chain)):
                for j in range (i+1, len(in_pose_by_chain)):
                    closure_results_dic[i,j]=[]
                    closure_results_dic[j,i]=[]
                    
        #Create an array to hold the potential closures possible
        binary_closure_matrix=np.zeros((in_num_chains, in_num_chains), bool)
        binary_closure_matrix.fill(False)
        total_num_loop_solutions=0
        
        for ichndx in range(len(in_pose_by_chain)):
            for iChainDelta in range(allowVariateSSsizeBy):
                for jchndx in range(len(in_pose_by_chain)):
                    for jChainDelta in range(allowVariateSSsizeBy):
                        ###print "TEST:", ichndx, jchndx
                        if (ichndx==jchndx):
                            continue
                        if(only_as_input_conectivity_order_mode and (ichndx != (jchndx-1))):
                            continue
                        #Here starts the action
                        possible_closure_ndxs=[]
                        #Find by NC distances the loops that are near to the solution
                        if (closure_by_dist_matrix[ichndx][jchndx][iChainDelta][jChainDelta] < self.nc_points_max_dist_threshold):
                            
                            possible_closure_ndxs=np.where((self.nScS_NCpoints_distances-closure_by_dist_matrix[ichndx][jchndx][iChainDelta][jChainDelta]) < 
                                                                self.loop_NCpoints_distance_max_threshold_to_try)[0]
 

                            #possible_closure_ndxs=np.where(abs(self.nScS_NCpoints_distances-closure_by_dist_matrix[ichndx][jchndx][iChainDelta][jChainDelta]) < 
                            #                                    self.loop_NCpoints_distance_max_threshold_to_try)[0]
                        else: #If not possible to close just quit
                            continue

                        #This array contains: chainA, resA, chainB, resB
                        AB_closureMode=True
                        closure_order=[ichndx+1, 
                                            in_pose_by_chain[ichndx+1].total_residue()-iChainDelta, 
                                            jchndx+1, 
                                            1+jChainDelta]
                        if (ichndx>jchndx):
                            AB_closureMode=False

                            
                        #print("A", self.nScS_NCpoints_distances)
                        #print("B", closure_by_dist_matrix[ichndx][jchndx][iChainDelta][jChainDelta])
                        #print("C", self.loop_NCpoints_distance_max_threshold_to_try)
                        #print("D", self.nScS_NCpoints_distances-closure_by_dist_matrix[ichndx][jchndx][iChainDelta][jChainDelta])
                        if self.b_debug_print:
                            print ("Now I'll try, number of potential closure modes for : ", 
                                    "C-%d D_%d"%(closure_order[0],closure_order[1]), 
                                    "C-%d D_%d"%(closure_order[2],closure_order[3]), 
                                    "=", 
                                    len(possible_closure_ndxs))#, 
                                    #", ndxs: ",
                                    #possible_closure_ndxs)
                            print "Dist closure_by_dist_matrix[ichndx][jchndx][iChainDelta][jChainDelta]", closure_by_dist_matrix[ichndx][jchndx][iChainDelta][jChainDelta] 
                            print "Corder:", closure_order

                        #print pose_by_chain_ss[closure_order[0]-1][closure_order[1]-1], 
                        #print pose_by_chain_ss[closure_order[2]-1][closure_order[3]-1]

                        #print "Total possibilities:", len(possible_closure_ndxs)
                        for indx in possible_closure_ndxs:
                            
                            ##print "loop, NC side:", self.nScS_NCpoints_ndxs[indx]
                            #These are the indexes for the fragments, must be shifted by 2 (+1 rosetta indexing, +1 extra residue)
                            rmsd_cmp_array_a_ndxs=np.arange(0, 
                                                            (self.nScS_NCpoints_ndxs[indx][0]+1))+2
                            rmsd_cmp_array_b_ndxs=np.arange(self.nScS_NCpoints_ndxs[indx][1], 
                                                            (self.cluster_lr_fragment_size))+2
                            
                            #loop_len here is actually the real loop len without the EP
                            loop_len=rmsd_cmp_array_b_ndxs[0]-rmsd_cmp_array_a_ndxs[-1]-1

                            #Note: not sure if we should handle the 0(zero)-case
                            #Note: this might be moved to the loop generation stage
                            if self.b_debug_print>1:
                                print("Loop len to try:", loop_len)
                            if(loop_len < max(0,self.min_allowed_loop_len)):
                                continue
                            elif(loop_len > self.max_allowed_loop_len):
                                continue

                            ##print "loop pose nc res ndxs: ", rmsd_cmp_array_a_ndxs, rmsd_cmp_array_b_ndxs
                            a_side_diff=len(rmsd_cmp_array_a_ndxs)
                            b_side_diff=len(rmsd_cmp_array_b_ndxs)
                            #These are the indexes for the pose N-side
                            rmsd_cmp_array_c_ndxs=[]
                            rmsd_cmp_array_d_ndxs=[]
                            #ToDO: Double check this! Possibly is wrong
                            #if AB_closureMode:
                            rmsd_cmp_array_c_ndxs=np.arange(closure_order[1]-(a_side_diff-1), closure_order[1]+1)
                            rmsd_cmp_array_d_ndxs=np.arange(closure_order[3], closure_order[3]+b_side_diff)

                            #ToDO: Remove negative numbers from c and d
                            #Note: This is to correct the ndxs if we overflow the SS fragment, might be possible to write it more beautifully :)
                            min_align_ndx=1
                            max_align_index_b=len(pose_by_chain_ss[closure_order[2]-1])
                            if ( rmsd_cmp_array_c_ndxs[0] < min_align_ndx ):
                                tmp_good_a_c_ndxs=np.where(rmsd_cmp_array_c_ndxs >= min_align_ndx)[0]
                                ##print "SKIPPING_C:", min_align_ndx, rmsd_cmp_array_a_ndxs, rmsd_cmp_array_c_ndxs, tmp_good_a_c_ndxs
                                rmsd_cmp_array_a_ndxs=rmsd_cmp_array_a_ndxs[tmp_good_a_c_ndxs]
                                rmsd_cmp_array_c_ndxs=rmsd_cmp_array_c_ndxs[tmp_good_a_c_ndxs]
                                ##print "NEW ACndxs:", rmsd_cmp_array_a_ndxs, rmsd_cmp_array_c_ndxs
                                #continue
                            if (rmsd_cmp_array_d_ndxs[-1] > max_align_index_b):
                                tmp_good_b_d_ndxs=np.where(rmsd_cmp_array_d_ndxs <= max_align_index_b)[0]
                                ##print "SKIPPING_D:", max_align_index_b, rmsd_cmp_array_b_ndxs, rmsd_cmp_array_d_ndxs, tmp_good_b_d_ndxs  
                                rmsd_cmp_array_b_ndxs=rmsd_cmp_array_b_ndxs[tmp_good_b_d_ndxs]
                                rmsd_cmp_array_d_ndxs=rmsd_cmp_array_d_ndxs[tmp_good_b_d_ndxs]
                                #print "NEW CDndxs:", rmsd_cmp_array_b_ndxs, rmsd_cmp_array_d_ndxs
                                #continue
                            
                            #Get the poses into the tmp_pointers
                            SS_pose_nside=pose_by_chain_ss[closure_order[0]-1][rmsd_cmp_array_c_ndxs-1]
                            SS_pose_cside=pose_by_chain_ss[closure_order[2]-1][rmsd_cmp_array_d_ndxs-1]                  

                            #This is to compare the SSs type matching!
                            SS_loop_nside=self.nScS_poses_ss[indx][rmsd_cmp_array_a_ndxs-2]
                            SS_loop_cside=self.nScS_poses_ss[indx][rmsd_cmp_array_b_ndxs-2]
                            #Hack, bypass SS comparision for now
                            #if ( not np.array_equal(SS_pose_nside, SS_loop_nside) or 
                            #    (not np.array_equal(SS_pose_cside, SS_loop_cside) ) ):
                            #    continue
                            
                            
                            tmp_ca_c_o_n_ndxs_a=[]
                            for kk in ((rmsd_cmp_array_a_ndxs-1)*4): #4==4atoms==CA,C,O,N
                                for ll in xrange(4):
                                    tmp_ca_c_o_n_ndxs_a.append(kk+ll)
                            tmp_ca_c_o_n_ndxs_b=[]
                            for kk in ((rmsd_cmp_array_b_ndxs-1)*4): #4==4atoms==CA,C,O,N
                                for ll in xrange(4):
                                    tmp_ca_c_o_n_ndxs_b.append(kk+ll)
                            tmp_ca_c_o_n_ndxs_c=[]
                            for kk in ((rmsd_cmp_array_c_ndxs-1)*4): #4==4atoms==CA,C,O,N
                                for ll in xrange(4):
                                    tmp_ca_c_o_n_ndxs_c.append(kk+ll)
                            tmp_ca_c_o_n_ndxs_d=[]
                            for kk in ((rmsd_cmp_array_d_ndxs-1)*4): #4==4atoms==CA,C,O,N
                                for ll in xrange(4):
                                    tmp_ca_c_o_n_ndxs_d.append(kk+ll)
                            
                            
                            tmp_coorA=np.vstack([self.nScS_poses_bbxyz[indx][tmp_ca_c_o_n_ndxs_a],
                                                                            self.nScS_poses_bbxyz[indx][tmp_ca_c_o_n_ndxs_b]])
                            tmp_coorB=np.vstack([in_pose_by_chain_bb_xyz[closure_order[0]-1][tmp_ca_c_o_n_ndxs_c],
                                                                            in_pose_by_chain_bb_xyz[closure_order[2]-1][tmp_ca_c_o_n_ndxs_d]])

                            
                            
                            #Thist aligns A->B
                            nc_points_rmsd,rMatA,tVecA,tVecB=rmsd_and_rotations_2_np_arrays(tmp_coorA, tmp_coorB)
                            if self.b_debug_print>1:
                                print "nc_points_rmsd:",nc_points_rmsd
                            if ( nc_points_rmsd > self.nc_points_rmsd_dist_threshold ) :
                                continue
                            """    
                            #OLD SLOW RMSD CODE (since it extracts the XYZ every time and Rosetta Sucks on that)
                            ##RMSD checking routine
                            test_loop_pose=self.nScS_poses[indx].clone()
                            nc_points_rmsd=self.align_bb_1pose_to_2poses_by_ndxs(test_loop_pose, 
                                                                    rmsd_cmp_array_a_ndxs, 
                                                                    rmsd_cmp_array_b_ndxs, 
                                                                    in_pose_by_chain[closure_order[0]], 
                                                                    rmsd_cmp_array_c_ndxs, 
                                                                    in_pose_by_chain[closure_order[2]], 
                                                                    rmsd_cmp_array_d_ndxs)
                            """
                            
                            #tmpCoorA and tmp_coorB return translated, so rotateA and you get perAtomDist
                            perAtom_dist_ABsides=np.linalg.norm(np.dot(tmp_coorA,rMatA)-tmp_coorB,axis=1)
                            #If per atom distance is too large skip
                            #ToDo; could make the threshold adaptive respect to the RMSD as in: min((nc_points_rmsd*2.0),self.nc_points_individual_msd_dist_threshold)
                            if self.b_debug_print>1:
                                print "perAtom_dist_ABsides", perAtom_dist_ABsides
                            if (max(perAtom_dist_ABsides)>self.nc_points_individual_msd_dist_threshold):
                                continue
                                
                            """
                            #OLD SLOW MSD CODE (since it extracts the XYZ every time and Rosetta Sucks on that)
                            #MAybe need a second RMSD alignment here using only the NC-points
                            if ( nc_points_rmsd > self.nc_points_rmsd_dist_threshold ) :
                                continue
     
                            #MSD checking routine
                            msdAside=self.msd_bb_by_ndxs(test_loop_pose, 
                                                        rmsd_cmp_array_a_ndxs[-1], 
                                                        rmsd_cmp_array_a_ndxs[-1], 
                                                        in_pose_by_chain[closure_order[0]], 
                                                        rmsd_cmp_array_c_ndxs[-1], 
                                                        rmsd_cmp_array_c_ndxs[-1])

                            msdBside=self.msd_bb_by_ndxs(test_loop_pose, 
                                                        rmsd_cmp_array_b_ndxs[0], 
                                                        rmsd_cmp_array_b_ndxs[0], 
                                                        in_pose_by_chain[closure_order[2]], 
                                                        rmsd_cmp_array_d_ndxs[0], 
                                                        rmsd_cmp_array_d_ndxs[0]) 
                            #Check MSDs
                            #print "The MSDs:", msdAside, msdBside
                            if ( (msdAside > 
                                    self.nc_points_individual_msd_dist_threshold) or 
                                (msdBside > 
                                    self.nc_points_individual_msd_dist_threshold) ):
                                continue
                            """
                            
                            ##print "Good:", nc_points_rmsd, max(perAtom_dist_ABsides)
                            
                            #Now copy the pose and translate it to where it should be as rosettaLikes
                            #ToDo: move this to a function
                            test_loop_pose=self.nScS_poses[indx].clone()
                            rMatA=rMatA.T
                            tmp_tVec=tVecB-(np.dot(tVecA, np.linalg.inv(rMatA)))
                            tmp_rMtx_xyzM=pyrosetta.rosetta.numeric.xyzMatrix_double_t()
                            tmp_rMtx_xyzM.xx=rMatA[0,0]
                            tmp_rMtx_xyzM.xy=rMatA[0,1]
                            tmp_rMtx_xyzM.xz=rMatA[0,2]
                            tmp_rMtx_xyzM.yx=rMatA[1,0]
                            tmp_rMtx_xyzM.yy=rMatA[1,1]
                            tmp_rMtx_xyzM.yz=rMatA[1,2]
                            tmp_rMtx_xyzM.zx=rMatA[2,0]
                            tmp_rMtx_xyzM.zy=rMatA[2,1]
                            tmp_rMtx_xyzM.zz=rMatA[2,2]
                            tmp_tVec_xyzV=pyrosetta.rosetta.numeric.xyzVector_double_t()
                            tmp_tVec_xyzV.x=tmp_tVec[0]
                            tmp_tVec_xyzV.y=tmp_tVec[1]
                            tmp_tVec_xyzV.z=tmp_tVec[2]
                            test_loop_pose.apply_transform_Rx_plus_v(tmp_rMtx_xyzM, 
                                                                     tmp_tVec_xyzV)
                            
                            ##in_pose_by_chain[closure_order[0]].dump_pdb("/home/dadriano/tmp/testA.pdb")
                            ##in_pose_by_chain[closure_order[2]].dump_pdb("/home/dadriano/tmp/testB.pdb")
                            ##test_loop_pose.dump_pdb("/home/dadriano/tmp/testC%04d.pdb"%indx)
                            #break
                            
                            #Check for clashes:
                            #We have already removed the 2-np-edges from the context_pose_hash_array
                            if (loop_len > 0): #If loop len is == 0 then this test will return clashes for sure
                                tmpline=""
                                for ires in range(loop_len):
                                    tmpline+=self.fake_aminoacid_letter_for_design
                                test_loop_pose_minus_ends=pyrosetta.pose_from_sequence(tmpline)
                                ##( (Pose)arg1, (int)size, (Pose)src, (int)begin, (int)src_begin)

                                #Double check this (loop_len)
                                test_loop_pose_minus_ends.copy_segment(loop_len,
                                                                        test_loop_pose,
                                                                        1,
                                                                        rmsd_cmp_array_a_ndxs[-1]+1)

                                hasBBclashes=False
                                for ihash in context_pose_hash_array:
                                    hasBBclashes, target_clash_map = self.check_pm_clashes_hash_vs_pose(ihash, 
                                                                                                        test_loop_pose_minus_ends,
                                                                                                        pick_mode=pyrosetta.rosetta.core.pose.PoseCoordPickMode_CA #Maybe N_CA_C
                                                                                                        )
                                    if (hasBBclashes): 
                                        break
                                if (hasBBclashes):
                                    if self.b_debug_print:
                                        print "Clashes from loop (ll=%d) with the target SSs:\n"%loop_len, test_loop_pose, rmsd_cmp_array_a_ndxs[-1]+1 
                                        print " Fail clashCheck: ", "clashloop_c%02d_%02d-c%02d_%02d_ll%04d_ci%05d_loop.pdb"%(closure_order[0],
                                                                                                                closure_order[1],
                                                                                                                closure_order[2],
                                                                                                                closure_order[3],
                                                                                                                loop_len,
                                                                                                                indx) 
                                    if self.b_debug_out:
                                        test_loop_pose_minus_ends.dump_pdb("%s/clashloop_c%02d_%02d-c%02d_%02d_ll%04d_ci%05d_loop.pdb"%(self.debug_out_path,
                                                                                                                closure_order[0],
                                                                                                                closure_order[1],
                                                                                                                closure_order[2],
                                                                                                                closure_order[3],
                                                                                                                loop_len,
                                                                                                                indx) )
                                    continue

                            #ToDO:
                            ##1. Move here check Clashes agains the Global, will speed up everything 
                            #Connect and minimize
                            [test_mini_pose, 
                             score_merge]= self.connect_3_poses_byNDX_and_minize_w_constraints(
                                                                                            pose1=in_pose_by_chain[closure_order[0]],  
                                                                                            pose2=test_loop_pose, 
                                                                                            pose3=in_pose_by_chain[closure_order[2]], 
                                                                                            ndx_p1b=rmsd_cmp_array_c_ndxs[-1], 
                                                                                            ndx_p2a=rmsd_cmp_array_a_ndxs[-1],
                                                                                            ndx_p2b=rmsd_cmp_array_b_ndxs[0],
                                                                                            ndx_p3a=rmsd_cmp_array_d_ndxs[0], 
                                                                                            context_pose=pyrosetta.rosetta.core.pose.Pose(),  
                                                                                            sfx=self.scorefxn_ramaHyper,
                                                                                            harmonic_constraint_streght=1.0, #Smaller the stronger 
                                                                                            cart_constraint_weight_ramp=[1.0, 
                                                                                                                        0.8, 
                                                                                                                        0.5, 
                                                                                                                        0.1,
                                                                                                                        0.0],
                                                                                            aminoacid_letter_for_design=self.fake_aminoacid_letter_for_design )

                            
                            ##test_mini_pose.dump_pdb("/home/dadriano/tmp/testD%04d.pdb"%indx)
                            
                            #Check OMEGA angles
                            #ToDo, build something that allows better reconstruction of the loops. Idea is the snake-algorithm
                            low_bound=max(1, rmsd_cmp_array_c_ndxs[-1]-1)
                            win_len=loop_len+2
                            high_bound=min(test_mini_pose.total_residue(), (rmsd_cmp_array_c_ndxs[-1]+win_len))
                            #print "Omega L/H:", low_bound, high_bound, range(low_bound, high_bound+1)
                            for ires in xrange(low_bound, high_bound+1):
                                tmp_ome=abs(test_mini_pose.omega(ires))
                                #Check OMEGA after minimization
                                #Possible omega is essentially only 180 and 0 in proteins, don't allow more than 20deg change, or so...
                                if((tmp_ome < (180-self.allowed_omega_deviation)) and (tmp_ome > (0+self.allowed_omega_deviation)) ): #(tmp_ome > (0+self.allowed_omega_deviation)) and
                                        if self.b_debug_print:
                                            print "OME DEVIATION from planar is too large. Skipping due to this: !!!!: ", tmp_ome, low_bound, high_bound+1
                                        continue
                            

                            ##Check Clashes agains the Global context if b_use_global_context:
                            if b_use_global_context:
                                hasGC_BBclashes, GC_clash_map = self.check_pm_clashes_hash_vs_pose(global_context_pose_hash, 
                                                                                                    test_mini_pose,
                                                                                                    pick_mode=pyrosetta.rosetta.core.pose.PoseCoordPickMode_N_CA_C_CB)
                                if (hasGC_BBclashes):
                                    if (((len(list(GC_clash_map)))-num_global_context_clashes_initial) > 0 ):
                                        if self.b_debug_print:
                                            print "Has Global context clashes: ", list(target_clash_map)
                                        continue

                            
                            #Check that the loop and pre-post structures do exists in the clustered database!
                            ##print "TEST: ", max(1, low_bound-self.cluster_lr_fragment_size), min(test_mini_pose.total_residue(), high_bound+self.cluster_lr_fragment_size)
                            [b_is_good_fragment, 
                             [faulty_fragment_rmsd_avg,
                              faulty_fragment_rmsd_lim,
                              faulty_fragment_population_avg]] = self.check_span_for_ideal_fragments(test_mini_pose, 
                                                                                                            max(1, low_bound-self.cluster_lr_fragment_size), 
                                                                                                            min(test_mini_pose.total_residue(), high_bound+self.cluster_lr_fragment_size), 
                                                                                                            default_allowed_blind_raddi=self.fragment_blind_raddi_cut,
                                                                                                            rmsd_pert=self.rmsd_pert_to_clustered_raddi, 
                                                                                                            max_rmsd_limit=self.max_rmsd_limit_spans_check,
                                                                                                            spans_min_population_cut=self.check_spans_population_cut,
                                                                                                            cluster_fragment_len=self.cluster_lr_fragment_size,
                                                                                                            cluster_centers_coordinates=self.cluster_lr_centers_coordinates,
                                                                                                            cluster_centers_store_data=self.cluster_lr_centers_store_data,
                                                                                                            rmsd_calculator=fast_rmsd_calc)

                            ##print b_is_good_fragment, faulty_fragment_rmsd_avg, faulty_fragment_population_avg
                            ##assert 0==1
                            if(not b_is_good_fragment):
                                if self.b_debug_print:
                                    print "Bad fragment. Skipping this. Population after reassignment to clusters is too loow (Curr,Exp,Pop): ", faulty_fragment_rmsd_avg, faulty_fragment_rmsd_lim, faulty_fragment_population_avg
                                continue

                                

                            ##If you reach this point it should be fine
                            if self.b_debug_print:
                                print "Found ONE:", ("#%d"%len(closure_results_dic[ichndx,jchndx]),
                                                        "%d<->%d"%(closure_order[0],closure_order[2]),
                                                        "LL: %d"%loop_len,
                                                        "cNdx: %d"%indx,
                                                        "RMSD: %0.3f"%nc_points_rmsd,
                                                        "MaxAtmDist: %0.3f"%max(perAtom_dist_ABsides),
                                                        "SpanRMSD: %0.3f"%faulty_fragment_rmsd_avg,
                                                        "MaxPopSum: %0.3f"%faulty_fragment_population_avg )
                                print ("*SS (LvsT):", self.nScS_poses_ss[indx],
                                        "%d:%d"%((rmsd_cmp_array_a_ndxs-2)[-1],(rmsd_cmp_array_b_ndxs-2)[0]),
                                        "=",self.nScS_poses_ss[indx][(rmsd_cmp_array_a_ndxs-2)[-1]:(rmsd_cmp_array_b_ndxs-2)[0]+1], 
                                        "vs", pose_by_chain_ss[closure_order[0]-1][-4:],"-N--C-", pose_by_chain_ss[closure_order[2]-1][:4])

                            #print "D", get_labeled_residues(this_inpose=test_mini_pose)

                            #Add this RESULT to the dictionary of results
                            closure_results_dic[ichndx,jchndx].append([loop_len, 
                                                                        test_mini_pose.clone(), 
                                                                        self.nScS_NCpoints_size[indx],
                                                                        nc_points_rmsd ])
                            

                            binary_closure_matrix[ichndx][jchndx]=True
                            total_num_loop_solutions+=1

                            if( total_num_loop_solutions >= max_num_loop_solutions):
                                break
                        if( total_num_loop_solutions >= max_num_loop_solutions):
                            break
                    if( total_num_loop_solutions >= max_num_loop_solutions):
                        break
                        

                


        if self.b_debug_print:
            print "All loop finding done"
            print "Please come back soon to make some more loops, because if not I feel so alone and unuseful..."
        
        return closure_results_dic
    
    
    def clustered_aaProbToDic(self,
                    probability_table=[]): 
        #Minimum check
        assert (len(probability_table)==20)
            
        #Be carefull, this might mess the universe if the order is actually wrong
        result_aa_prob_dic={}
        result_aa_prob_dic["A"]=probability_table[0]
        result_aa_prob_dic["R"]=probability_table[1]
        result_aa_prob_dic["N"]=probability_table[2]
        result_aa_prob_dic["D"]=probability_table[3]
        result_aa_prob_dic["C"]=probability_table[4]
        result_aa_prob_dic["E"]=probability_table[5]
        result_aa_prob_dic["Q"]=probability_table[6]
        result_aa_prob_dic["G"]=probability_table[7]
        result_aa_prob_dic["H"]=probability_table[8]
        result_aa_prob_dic["I"]=probability_table[9]
        result_aa_prob_dic["L"]=probability_table[10]
        result_aa_prob_dic["K"]=probability_table[11]
        result_aa_prob_dic["M"]=probability_table[12]
        result_aa_prob_dic["F"]=probability_table[13]
        result_aa_prob_dic["P"]=probability_table[14]
        result_aa_prob_dic["S"]=probability_table[15]
        result_aa_prob_dic["T"]=probability_table[16]
        result_aa_prob_dic["W"]=probability_table[17]
        result_aa_prob_dic["Y"]=probability_table[18]
        result_aa_prob_dic["V"]=probability_table[19]
        
        return result_aa_prob_dic


    def __init__(self,
                 nc_points_rmsd_dist_threshold=0.9,
                 nc_points_individual_msd_dist_threshold=0.9,
                 loop_NCpoints_distance_max_threshold_to_try=1.0,
                 max_rmsd_limit_spans_check=1.0,
                 allowed_omega_deviation=30.0,
                 cluster_population_cut_off=100,
                 check_spans_population_cut=50,
                 fragment_blind_raddi_cut=0.5,
                 fake_aminoacid_letter_for_design="A",
                 min_allowed_loop_len=0,
                 max_allowed_loop_len=99,
                 min_aa_prob_for_label=0.05,
                 minTrimmedSSsize=4,
                 max_zero_CA_CA_distances=[4.0,3.0,2.0], #CA-CA max distance for considering it closable (has to contain ~max and min allowable/thinkable scenarios)
                 b_only_close_with_loop=False,
                 b_generate_loop_loop_closures=False,
                 b_debug_print=False,
                 b_debug_out=False,
                 clustered_fragments_path="./size_7/kcenters_stats.dat_r1",
                 clustered_fragments_assignments_path="./size_7/kcenters_assignments.dat_r1",
                 debug_out_path="./"):
        self.status=0
        print "Initializing Dloop sampling algorithm"
        
        #Debug print options
        self.b_debug_print=b_debug_print
        if self.b_debug_print:
            print "Enabled nasty debug prints"
        self.b_debug_out=b_debug_out
        self.debug_out_path=debug_out_path
        if self.b_debug_out: 
            print "Enabled debug out of structures to path: ", self.debug_out_path
            
        ##Sampler Reducer options
        ##self.keep_only_top_loops_num_max=200
        ###self.max_num_results_when_rebuilding=100
        
        #The clustered fragments Database:
        self.clustered_fragments_path=clustered_fragments_path
        self.clustered_fragments_assignments_path=clustered_fragments_assignments_path
        
        #ToDo: make this/next to be automatically determined from the database
        self.nc_points_rmsd_dist_threshold=nc_points_rmsd_dist_threshold  #1.5 permisive, 1.3 medium, 1.1 strict
        self.nc_points_individual_msd_dist_threshold=nc_points_individual_msd_dist_threshold #1.5 permisive, 1.3 medium, 1.1 strict
        self.loop_NCpoints_distance_max_threshold_to_try=loop_NCpoints_distance_max_threshold_to_try
        self.allowed_omega_deviation=allowed_omega_deviation  #Angle, 30.0 OK
        #self.min_allowed_score_per_res_before_des=0.0
        
        ##only_as_input_conectivity_order_mode=True  #True uses input order, False uses exhaustive ordering
        self.min_allowed_loop_len=min_allowed_loop_len #Values smaller than this will be ignored
        self.max_allowed_loop_len=max_allowed_loop_len
        self.cluster_population_cut_off=cluster_population_cut_off
        self.check_spans_population_cut=check_spans_population_cut
        #allowVariateSSsizeBy=1 #8 Deprecated to =0
        self.minTrimmedSSsize=minTrimmedSSsize
        
        #spans check
        self.max_rmsd_limit_spans_check=max_rmsd_limit_spans_check #Make this an option
        self.rmsd_pert_to_clustered_raddi=1.1 #check this
        self.fragment_blind_raddi_cut=fragment_blind_raddi_cut
        #break_trimming_at_hotspots=False
        
        #Misc
        self.fake_aminoacid_letter_for_design=fake_aminoacid_letter_for_design
        self.min_aa_prob_for_label=min_aa_prob_for_label
        
        #Database options
        self.b_only_close_with_loop=b_only_close_with_loop
        if self.b_only_close_with_loop:
                print "Closure only with loop types"
        else:
                print "Warning: I'll find closures not only with loops but using SSs types also."
                
        self.b_generate_loop_loop_closures=b_generate_loop_loop_closures
        if self.b_generate_loop_loop_closures:
                print "Closure only of SS types"
        else:
                print "!!!Warning: I'll also try to close Nloop<-->Cloop, this might increase your waiting time a loooot (database exploding)."
        
        #Scoring functions
        #self.scorefxn_soft_rep_name="beta_cart"
        #self.scorefxn_soft_rep = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function(self.scorefxn_soft_rep_name)
        self.scorefxn_vanilla_name="beta_cart"
        self.scorefxn_vanilla = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function(self.scorefxn_vanilla_name)
        
        
        self.scorefxn_ramaHyper_name="beta_cart"
        self.scorefxn_ramaHyper = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function(self.scorefxn_ramaHyper_name)
        self.scorefxn_ramaHyper.set_weight(pyrosetta.rosetta.core.scoring.rama , self.scorefxn_ramaHyper.get_weight(pyrosetta.rosetta.core.scoring.rama)*3.0)
        self.scorefxn_ramaHyper.set_weight(pyrosetta.rosetta.core.scoring.hbond_bb_sc , self.scorefxn_ramaHyper.get_weight(pyrosetta.rosetta.core.scoring.hbond_bb_sc)*3.0)
        self.scorefxn_ramaHyper.set_weight(pyrosetta.rosetta.core.scoring.hbond_lr_bb , self.scorefxn_ramaHyper.get_weight(pyrosetta.rosetta.core.scoring.hbond_lr_bb)*3.0)
        self.scorefxn_ramaHyper.set_weight(pyrosetta.rosetta.core.scoring.hbond_sc , self.scorefxn_ramaHyper.get_weight(pyrosetta.rosetta.core.scoring.hbond_sc)*3.0)
        #print scorefxn_ramaHyper
        
        
        #RMSD calculator
        #ultrafast RMSD calculator
        #self.num_threads_global=5
        #self.fast_rmsd_calculator= #pyrosetta.rosetta.core.numeric. #rmsd_calc(self.num_threads_global)
        
        #Read Database
        print "Reading clustered fragments database"
        #LR clustering
        
        [self.cluster_lr_centers_store_data,
            self.cluster_lr_fragment_size]= self.parse_dadriano_clustered_fragment_table(self.clustered_fragments_path, 
                                                                                                    atom_types_arr=["CA", "C", "O", "N"])
        self.cluster_lr_centers_coordinates = self.cluster_lr_centers_store_data["coordinates"]
        self.cluster_lr_centers_coordinates = self.cluster_lr_centers_coordinates.copy().view(
                                                dtype=float).reshape((self.cluster_lr_centers_coordinates.shape[0], -1, 3))
        self.num_clusters_lr=len(self.cluster_lr_centers_coordinates)
        print "Readed a total of %d clusters"%self.num_clusters_lr
        print "The fragment size is %d"%self.cluster_lr_fragment_size
        
        """
        #Now Read assignments (dadriano version)
        print "Reading assignments (LR and HR)"
        #LR assignments
        self.assignments_lr_store_data = self.parse_dadriano_clustered_fragment_assignments(
                                        clustered_fragments_assignments_path)
        print "Done"
        #Get contiguous indexes
        tmpIndex=0
        self.assignments_lr_contiguousIndex= np.zeros(len(self.assignments_lr_store_data["assignment"]),int)
        self.assignments_lr_contiguousIndex[0]=tmpIndex
        for i in range(0,len(self.assignments_lr_store_data["assignment"])-1):
                if ( (self.assignments_lr_store_data["span"][i,0]) == (self.assignments_lr_store_data["span"][i+1,0]-1) ):
                    self.assignments_lr_contiguousIndex[i+1]=tmpIndex
                else:
                    tmpIndex+=1
                    self.assignments_lr_contiguousIndex[i+1]=tmpIndex
        print "Number of contiguous LR fragments: ", self.assignments_lr_contiguousIndex.max()+1
        
        #Calculate TCM with jump step=self.cluster_lr_fragment_size-3
        self.tcm_jump_step=self.cluster_lr_fragment_size-3
        self.clusters_lr_TCM = np.zeros((self.num_clusters_lr,self.num_clusters_lr), int)
        for i in range(0,len(self.assignments_lr_store_data["assignment"])-self.tcm_jump_step):
                if ( self.assignments_lr_contiguousIndex[i] == self.assignments_lr_contiguousIndex[i+self.tcm_jump_step]):
                    self.clusters_lr_TCM[(self.assignments_lr_store_data["assignment"][i]-1)][(self.assignments_lr_store_data["assignment"][i+self.tcm_jump_step]-1)]+=1
        print "LR TCM created: ", self.clusters_lr_TCM.shape, ", with jump step: ", self.tcm_jump_step
        """

        self.nScS_cluster_ndxs=[]
        self.sstype_lib_dic={}
        #H= index 0; E= index 1; L= index 2
        for i in np.where(self.cluster_lr_centers_store_data["size"]>=self.cluster_population_cut_off)[0]:
                #To remove those starting with loops
                max_ss_prob=self.cluster_lr_centers_store_data["ss"][i].argmax(axis=1)
                #Hack convert ss=L to ss=E
                ##tmp_ndx_loop=np.where(max_ss_prob==2)
                ##max_ss_prob(tmp_ndx_loop)=1
                #End HAck to convert ss=L to ss=E
                if (self.b_only_close_with_loop) and (2 in max_ss_prob[1:-1]):
                    self.nScS_cluster_ndxs.append(i)
                    self.sstype_lib_dic[i]=max_ss_prob[:]
                else:
                    self.nScS_cluster_ndxs.append(i)
                    self.sstype_lib_dic[i]=max_ss_prob[:]
                    
        print "The number of potential loop closures from the Database that will be considered based in your options is: %d"% len(self.nScS_cluster_ndxs)
        print "ALL DONE"


        ###

        #DO NOT Create poses for all the clusters, just check omega angle
        #self.all_consider_clusters_poses=[]
        print "Removing possible closures with bad omega angles"
        self.all_consider_clusters_poses_omega_quality=[]
        for lr_fragment in xrange(self.num_clusters_lr):
                if lr_fragment in self.nScS_cluster_ndxs:
                    b_is_good_pose=True
                    for angle in self.cluster_lr_centers_store_data["bb"] [lr_fragment]:
                        if (abs(angle[2])<(180-self.allowed_omega_deviation)):
                                b_is_good_pose=False
                                break
                        
                    #self.all_consider_clusters_poses.append(tmp_pose.clone())
                    self.all_consider_clusters_poses_omega_quality.append(b_is_good_pose)
                else:
                    #self.all_consider_clusters_poses.append([])
                    self.all_consider_clusters_poses_omega_quality.append(False)
        self.all_consider_clusters_poses_omega_quality=np.asarray(self.all_consider_clusters_poses_omega_quality)
        #print "Not created poses for all the clusters. Total =", len(self.all_consider_clusters_poses), len(self.all_consider_clusters_poses_omega_quality)
        ##self.all_consider_clusters_poses=np.asarray(self.all_consider_clusters_poses)



        ###
        #Remove Faulty omega angles
        print "The number of potential loop closures before removing bad OMEGA angles is: %d"% len(self.nScS_cluster_ndxs)
        ##print len(where(self.all_consider_clusters_poses_omega_quality[self.nScS_cluster_ndxs]==False)[0])
        bad_omega_nScS_indexes=np.where(self.all_consider_clusters_poses_omega_quality[self.nScS_cluster_ndxs]==False)[0]
        for index in sorted(bad_omega_nScS_indexes, reverse=True):
                del self.nScS_cluster_ndxs[index]
        ##print len(where(self.all_consider_clusters_poses[self.nScS_cluster_ndxs]==False)[0])

        print "The number of potential loop closures after removing bad OMEGA angles is: %d"% len(self.nScS_cluster_ndxs)
        print "ALL DONE"


        print "Pregenerating the library of usable loops"
        ###
        #Pre-make a library of the loops
        self.nScS_poses=[]
        self.nScS_poses_bbxyz=[]
        self.nScS_poses_ss=[]
        self.nScS_NCpoints_distances=[]
        self.nScS_NCpoints_ndxs=[]
        self.nScS_NCpoints_size=[]
        
        for lr_frag_num,lr_fragment in enumerate(self.nScS_cluster_ndxs):
            #The pose container
            tmp_seq=""
            for i in xrange(self.cluster_lr_fragment_size+2):
                    tmp_seq+=self.fake_aminoacid_letter_for_design
            tmp_pose=pyrosetta.pose_from_sequence(tmp_seq)
            
            curr_res_ndx=2
            
            max_ss_prob=self.cluster_lr_centers_store_data["ss"][lr_fragment].argmax(axis=1)
            #Hack convert ss=L to ss=E
            ##tmp_ndx_loop=np.where(max_ss_prob==2)
            ##max_ss_prob(tmp_ndx_loop)=1
            #End HAck to convert ss=L to ss=E
            for angle in self.cluster_lr_centers_store_data["bb"] [lr_fragment]:
                tmp_pose.set_phi(curr_res_ndx, angle[0])
                tmp_pose.set_psi(curr_res_ndx, angle[1])
                tmp_pose.set_omega(curr_res_ndx, angle[2])
                curr_res_ndx+=1
                
            
            #Add 2 reslabels to all the residues
            #ToDo: Define the best way to use this, maybe separated NOTAA are easier to use under rosetta's "&&" logic
            for nndx in xrange(1, curr_res_ndx-1): #range(1, curr_res_ndx-1)
                ##array number is -1
                this_aa_prob_dic=self.clustered_aaProbToDic(self.cluster_lr_centers_store_data["aa"][lr_fragment][nndx-1])
                tmp_reslabel_aas=""
                for key in this_aa_prob_dic:
                    if (this_aa_prob_dic[key] >= self.min_aa_prob_for_label):
                        tmp_reslabel_aas=tmp_reslabel_aas+key
                assert(len(tmp_reslabel_aas)>0)
                tmp_reslabel_aas="PIKAA_"+tmp_reslabel_aas
                #Rosetta number is +1
                tmp_pose.pdb_info().add_reslabel((nndx+1), "DLOOPER")
                tmp_pose.pdb_info().add_reslabel((nndx+1), "PIKAA")
                tmp_pose.pdb_info().add_reslabel((nndx+1), tmp_reslabel_aas)
            
            ##is_loop_b=False
            for nndx in xrange(1, curr_res_ndx-1): #range(1, curr_res_ndx-1)
                #if is_loop_b:
                #    break
                #print lr_fragment, nndx-1
                #print self.sstype_lib_dic[lr_fragment][nndx-1]
                
                #This in principle is to consider only those posse that are loops that connect two SSs
                if ( (self.sstype_lib_dic[lr_fragment][nndx-1] == 2 ) and 
                     (not b_generate_loop_loop_closures) ):
                    ##is_loop_b=True
                    break
                for cndx in xrange(curr_res_ndx-2, nndx, -1): #xrange(curr_res_ndx-2, nndx, -1)
                    #This is, in principle, to consider only those closures that are loops that connect two SSs
                    if ( (self.sstype_lib_dic[lr_fragment][cndx-1] == 2 ) and 
                         (not b_generate_loop_loop_closures) ): #WT*!?? Double check if this is intended
                        break
                    tmp_dist=( (tmp_pose.residue(nndx+1).xyz( "CA" ) - 
                                                    tmp_pose.residue(cndx+1).xyz( "CA" )).norm() )
                        
                    #Possible ToDoS: 
                    #1. Use the next code to avoid loops shorter than what we want
                    # loop_len=cndx-nndx
                    # if(loop_len < max(1,self.min_allowed_loop_len)):
                    #      continue
                    #  elif(loop_len > self.max_allowed_loop_len):
                    #      continue
                    # print "loop size:", [nndx-1, cndx-1], cndx-nndx
                    self.nScS_poses.append(tmp_pose.clone())
                    self.nScS_poses_bbxyz.append(pose_atoms_to_nparray(tmp_pose,
                                                                       ["CA","C","O","N"])) #This is for speed
                    self.nScS_poses_ss.append(max_ss_prob[:])
                    self.nScS_NCpoints_distances.append( tmp_dist )
                    self.nScS_NCpoints_ndxs.append([nndx-1, cndx-1])
                    self.nScS_NCpoints_size.append(self.cluster_lr_centers_store_data["size"][lr_fragment])
                    
            """#Not needed, this comes for free from size-2 fragments
            #Add zero-size closure too:
            if ((lr_frag_num)==0) and ([0, 0] not in self.nScS_NCpoints_ndxs):
                for max_zero_CA_CA_dist in max_zero_CA_CA_distances:
                    print "Adding zero-size conection with limit of %0.2f A" %max_zero_CA_CA_dist
                    self.nScS_poses.append(tmp_pose.clone())
                    self.nScS_poses_bbxyz.append(pose_atoms_to_nparray(tmp_pose,
                                                                       ["CA","C","O","N"])) #This is for speed
                    self.nScS_poses_ss.append(max_ss_prob[:])
                    self.nScS_NCpoints_distances.append( max_zero_CA_CA_dist )
                    self.nScS_NCpoints_ndxs.append([0, 0])
                    self.nScS_NCpoints_size.append(9999999) #This connection is infinitely-true possible
            """
                    
        self.nScS_NCpoints_distances=np.asarray(self.nScS_NCpoints_distances)
        self.nc_points_max_dist_threshold=self.nScS_NCpoints_distances.max()
        print "Maximum closure distance using this database =",  self.nc_points_max_dist_threshold
        print "Number of closure-loops in database after expansion of all possible closure lengths=", len(self.nScS_NCpoints_distances)

        print "All init done!"
        self.status=1
        
        ###Init end

print "Done"

#RMSD and clustering routines!
def msd_bb_2_poses(pose1, 
                     init_res1, 
                     end_res1, 
                     pose2, 
                     init_res2, 
                     end_res2):
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

def msd_bb_2_poses_by_common_ndx(pose1, 
                                 pose2, 
                                 common_pose_ndxs):
    numRes=len(common_pose_ndxs)
    coorA=np.zeros(((4*numRes),3), float)
    coorB=np.zeros(((4*numRes),3), float)
    
    counter=0
    #print pose1
    #print pose2
    for indx in xrange(numRes):
        
        for dim in range(0,3): 
            coorA[counter,dim]=(pose1.residue(common_pose_ndxs[indx]).xyz("CA")[dim])
            coorB[counter,dim]=(pose2.residue(common_pose_ndxs[indx]).xyz("CA")[dim])
        counter+=1
        for dim in range(0,3): 
            coorA[counter,dim]=(pose1.residue(common_pose_ndxs[indx]).xyz("C")[dim])
            coorB[counter,dim]=(pose2.residue(common_pose_ndxs[indx]).xyz("C")[dim])
        counter+=1
        for dim in range(0,3): 
            coorA[counter,dim]=(pose1.residue(common_pose_ndxs[indx]).xyz("O")[dim])
            coorB[counter,dim]=(pose2.residue(common_pose_ndxs[indx]).xyz("O")[dim])
        counter+=1
        for dim in range(0,3): 
            coorA[counter,dim]=(pose1.residue(common_pose_ndxs[indx]).xyz("N")[dim])
            coorB[counter,dim]=(pose2.residue(common_pose_ndxs[indx]).xyz("N")[dim])
        counter+=1
        
    #Calculate the MSD
    msdVal=0.0
    for i in range(len(coorA)):
            tmpVal_buff= np.linalg.norm(coorA[i] - coorB[i])
            msdVal+=(tmpVal_buff*tmpVal_buff)
    msdVal=msdVal/len(coorA)
    msdVal=np.sqrt(msdVal)
    
    return msdVal

def msd_significantDifferences_bb_2_poses(pose1, 
                                         init_res1, 
                                         end_res1, 
                                         pose2, 
                                         init_res2, 
                                         end_res2,
                                         significant_val_threshold=0.01):
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
    tmp_atom_counter=0.0
    ##print "TEST", len(coorA), end_res1-init_res1
    for i in range(0,len(coorA),4):
        tmp_msdVal=0.0
        #Only for complete residues/bb
        for j in xrange(4):
            tmpVal_buff= np.linalg.norm(coorA[i+j] - coorB[i+j])
            tmp_msdVal+=(tmpVal_buff*tmpVal_buff)
        ###print tmp_msdVal
        if(tmp_msdVal>significant_val_threshold):
            msdVal+=tmp_msdVal
            tmp_atom_counter+=4.0
    ###print "TEST: ", msdVal, tmp_atom_counter
    msdVal=msdVal/tmp_atom_counter
    msdVal=np.sqrt(msdVal)
    return msdVal



def k_center_clusterer_poses(in_poses=[],
                             align_method="rmsd",
                             convergence_dist=0.5,
                             min_clusters_num=1,
                             max_clusters_num=9999,
                             pdb_info_label="DLOOPER",
                             b_find_true_cluster_center_after_kcenters=True):
    assert (convergence_dist >= 0.0)
    assert (len(in_poses)>1)
    curr_center=0
    data_size=len(in_poses)
    pose_size=in_poses[0][0].total_residue()
    print "Executing k-centers clustering in %d poses  with len=%d, with convergence criteria of %0.3f A :)"%(data_size,
                                                                                                              pose_size,
                                                                                                             convergence_dist)
    
    min_distances=np.zeros(data_size, float)
    assert (min_clusters_num >= 1)
    assert (min_clusters_num <= data_size)
    min_distances.fill(np.finfo(float).max)
    assignments=np.zeros(data_size, int)
    assignments.fill(curr_center)
    center_ndxs=[]
    if align_method=="rmsd":
        print "Using pose alignment (a.k.a RMSD)"
        while ( (min_distances.max() > convergence_dist) or
                (len(center_ndxs) < min_clusters_num) ):
            for indx in xrange(data_size):
                if(indx == curr_center):
                    min_distances[indx]=0.0
                    assignments[indx]=curr_center
                    continue
                tmp_rmsd= rmsd_atoms_by_ndxs(in_poses[curr_center][0],
                                          1,
                                          in_poses[curr_center][0].total_residue(),
                                          in_poses[indx][0],
                                          1,
                                          in_poses[indx][0].total_residue(),
                                            atoms=["CA","C","O","N"])
                if (tmp_rmsd < min_distances[indx]):
                    min_distances[indx]=tmp_rmsd
                    assignments[indx]=curr_center
            center_ndxs.append(curr_center)
            curr_center=min_distances.argmax()
            if ( len(center_ndxs) >= max_clusters_num):
                break
        ##
        ##START REASSIGN INSIDE CLUSTER
        if b_find_true_cluster_center_after_kcenters:
            print "Reassigning centers to true centers of the cluster"
            new_center_dic={}
            for ic in xrange(len(center_ndxs)):
                icenter=center_ndxs[ic]
                test_ndxs=np.where(assignments==icenter)[0]
                partial_dist_matrix=np.zeros((len(test_ndxs),len(test_ndxs)),
                                             float)
                new_center_internal_ndx=-1
                if (len(test_ndxs)>2):
                    for jn in xrange(len(test_ndxs)):
                        jndx=test_ndxs[jn]
                        for kn in xrange(jn, len(test_ndxs)):
                            kndx=test_ndxs[kn]
                            partial_dist_matrix[jn,kn]= rmsd_atoms_by_ndxs(in_poses[jndx][0],
                                                                      1,
                                                                      in_poses[jndx][0].total_residue(),
                                                                      in_poses[kndx][0],
                                                                      1,
                                                                      in_poses[kndx][0].total_residue(),
                                                                          atoms=["CA","C","O","N"])
                            partial_dist_matrix[kn,jn]=partial_dist_matrix[jn,kn]
                            
                    new_center_internal_ndx=partial_dist_matrix.sum(axis=0).argmin()
                    center_ndxs[ic]=test_ndxs[new_center_internal_ndx]
                    for jn in xrange(len(test_ndxs)):
                        jndx=test_ndxs[jn]
                        assignments[jndx]=center_ndxs[ic]
                        min_distances[jndx]=partial_dist_matrix[new_center_internal_ndx][jn]
                else:
                    new_center_dic[icenter]=icenter  
        ##END REASSIGN INSIDE CLUSTER 
        ##
    
    elif align_method=="msd":
        print "Not aligning poses!!! (a.k.a MSD)"
        while ( (min_distances.max() > convergence_dist) or
                (len(center_ndxs) < min_clusters_num) ):
            ##print "ASKJHDF"
            for indx in xrange(data_size):
                if(indx == curr_center):
                    min_distances[indx]=0.0
                    assignments[indx]=curr_center
                    continue
                tmp_msd= msd_bb_2_poses(in_poses[curr_center][0],
                                          1,
                                          in_poses[curr_center][0].total_residue(),
                                          in_poses[indx][0],
                                          1,
                                          in_poses[indx][0].total_residue())
                if (tmp_msd < min_distances[indx]):
                    min_distances[indx]=tmp_msd
                    assignments[indx]=curr_center
            center_ndxs.append(curr_center)
            curr_center=min_distances.argmax()
            if ( len(center_ndxs) >= max_clusters_num):
                break
        ##
        ##START REASSIGN INSIDE CLUSTER
        if b_find_true_cluster_center_after_kcenters:
            print "Reassigning centers to true centers of the cluster"
            new_center_dic={}
            for ic in xrange(len(center_ndxs)):
                icenter=center_ndxs[ic]
                test_ndxs=np.where(assignments==icenter)[0]
                partial_dist_matrix=np.zeros((len(test_ndxs),len(test_ndxs)),
                                             float)
                new_center_internal_ndx=-1
                if (len(test_ndxs)>2):
                    for jn in xrange(len(test_ndxs)):
                        jndx=test_ndxs[jn]
                        for kn in xrange(jn,len(test_ndxs)):
                            kndx=test_ndxs[kn]
                            
                            partial_dist_matrix[jn,kn]= msd_bb_2_poses(in_poses[jndx][0],
                                                                      1,
                                                                      in_poses[jndx][0].total_residue(),
                                                                      in_poses[kndx][0],
                                                                      1,
                                                                      in_poses[kndx][0].total_residue())
                            partial_dist_matrix[kn,jn]=partial_dist_matrix[jn,kn]
                            
                    new_center_internal_ndx=partial_dist_matrix.sum(axis=0).argmin()
                    center_ndxs[ic]=test_ndxs[new_center_internal_ndx]
                    for jn in xrange(len(test_ndxs)):
                        jndx=test_ndxs[jn]
                        assignments[jndx]=center_ndxs[ic]
                        min_distances[jndx]=partial_dist_matrix[new_center_internal_ndx][jn]
                else:
                    new_center_dic[icenter]=icenter  
        ##END REASSIGN INSIDE CLUSTER 
        ##
    
    elif align_method=="msd_label":
        print "Using only labeled residues. Not aligning poses!!! (a.k.a MSD of different parts)"
        to_align_positions=np.zeros((pose_size),bool)
        print "Detecting labeled residues as:", pdb_info_label
        for ipose in xrange(data_size):
            tmp_pdbinfo=in_poses[curr_center][0].pdb_info()
            for ires in xrange(pose_size):
                if tmp_pdbinfo.res_haslabel((ires+1), pdb_info_label):
                    to_align_positions[ires]=True
        to_align_res=np.where(to_align_positions==True)[0]
        if (len(to_align_res)> 0):
            print "Align residues ndx:", to_align_res
            while ( (min_distances.max() > convergence_dist) or
                    (len(center_ndxs) < min_clusters_num) ):
                for indx in xrange(data_size):
                    if(indx == curr_center):
                        min_distances[indx]=0.0
                        assignments[indx]=curr_center
                        continue
                    tmp_msd= msd_bb_2_poses_by_common_ndx(in_poses[curr_center][0], 
                                                         in_poses[indx][0], 
                                                         common_pose_ndxs=to_align_res+1)
                    if (tmp_msd < min_distances[indx]):
                        min_distances[indx]=tmp_msd
                        assignments[indx]=curr_center
                center_ndxs.append(curr_center)
                curr_center=min_distances.argmax()
                if ( len(center_ndxs) >= max_clusters_num):
                    break
            ##
            ##START REASSIGN INSIDE CLUSTER
            if b_find_true_cluster_center_after_kcenters:
                print "Reassigning centers to true centers of the cluster"
                new_center_dic={}
                for ic in xrange(len(center_ndxs)):
                    icenter=center_ndxs[ic]
                    test_ndxs=np.where(assignments==icenter)[0]
                    partial_dist_matrix=np.zeros((len(test_ndxs),len(test_ndxs)),
                                                 float)
                    new_center_internal_ndx=-1
                    if (len(test_ndxs)>2):
                        for jn in xrange(len(test_ndxs)):
                            jndx=test_ndxs[jn]
                            for kn in xrange(len(test_ndxs)):
                                kndx=test_ndxs[kn]
                                
                                partial_dist_matrix[jn,kn]= msd_bb_2_poses_by_common_ndx(in_poses[jndx][0], 
                                                                                     in_poses[kndx][0], 
                                                                                     common_pose_ndxs=to_align_res+1)
                        new_center_internal_ndx=partial_dist_matrix.sum(axis=0).argmin()
                        center_ndxs[ic]=test_ndxs[new_center_internal_ndx]
                        for jn in xrange(len(test_ndxs)):
                            jndx=test_ndxs[jn]
                            assignments[jndx]=center_ndxs[ic]
                            min_distances[jndx]=partial_dist_matrix[new_center_internal_ndx][jn]
                    else:
                        new_center_dic[icenter]=icenter  
            ##END REASSIGN INSIDE CLUSTER 
            ##
        elif (len(to_align_res) == 0):
            print "WARNING!! Non labeled residues (label==%s) found"%pdb_info_label
            print "Fallback: Switching automatically to full BB RMSD mode"
            print "Using pose alignment (a.k.a RMSD)"
            while ( (min_distances.max() > convergence_dist) or
                    (len(center_ndxs) < min_clusters_num) ):
                for indx in xrange(data_size):
                    if(indx == curr_center):
                        min_distances[indx]=0.0
                        assignments[indx]=curr_center
                        continue
                    tmp_rmsd= rmsd_atoms_by_ndxs(in_poses[curr_center][0],
                                              1,
                                              in_poses[curr_center][0].total_residue(),
                                              in_poses[indx][0],
                                              1,
                                              in_poses[indx][0].total_residue(),
                                                atoms=["CA","C","O","N"])
                    if (tmp_rmsd < min_distances[indx]):
                        min_distances[indx]=tmp_rmsd
                        assignments[indx]=curr_center
                center_ndxs.append(curr_center)
                curr_center=min_distances.argmax()
                if ( len(center_ndxs) >= max_clusters_num):
                    break
            ##
            ##START REASSIGN INSIDE CLUSTER
            if b_find_true_cluster_center_after_kcenters:
                print "Reassigning centers to true centers of the cluster"
                new_center_dic={}
                for ic in xrange(len(center_ndxs)):
                    icenter=center_ndxs[ic]
                    test_ndxs=np.where(assignments==icenter)[0]
                    partial_dist_matrix=np.zeros((len(test_ndxs),len(test_ndxs)),
                                                 float)
                    new_center_internal_ndx=-1
                    if (len(test_ndxs)>2):
                        for jn in xrange(len(test_ndxs)):
                            jndx=test_ndxs[jn]
                            for kn in xrange(len(test_ndxs)):
                                kndx=test_ndxs[kn]
                                
                                partial_dist_matrix[jn,kn]= rmsd_atoms_by_ndxs(in_poses[jndx][0],
                                                                              1,
                                                                              in_poses[jndx][0].total_residue(),
                                                                              in_poses[kndx][0],
                                                                              1,
                                                                              in_poses[kndx][0].total_residue(),
                                                                              atoms=["CA","C","O","N"])
                        new_center_internal_ndx=partial_dist_matrix.sum(axis=0).argmin()
                        center_ndxs[ic]=test_ndxs[new_center_internal_ndx]
                        for jn in xrange(len(test_ndxs)):
                            jndx=test_ndxs[jn]
                            assignments[jndx]=center_ndxs[ic]
                            min_distances[jndx]=partial_dist_matrix[new_center_internal_ndx][jn]
                    else:
                        new_center_dic[icenter]=icenter  
            ##END REASSIGN INSIDE CLUSTER 
            ##
            
        else:
            print "This is clearly a mistake sinnce the num of residues to align is <=0, stopping now!!!"
            assert(0==1)
            
        """
        elif align_method=="msd_significant_differences":
            print "Not aligning poses!!! (a.k.a MSD of different parts)"
            while ( (min_distances.max() > convergence_dist) or
                    (len(center_ndxs) < min_clusters_num) ):
                for indx in xrange(data_size):
                    if(indx == curr_center):
                        min_distances[indx]=0.0
                        assignments[indx]=curr_center
                        continue
                    tmp_msd= msd_significantDifferences_bb_2_poses(in_poses[curr_center][0],
                                              1,
                                              in_poses[curr_center][0].total_residue(),
                                              in_poses[indx][0],
                                              1,
                                              in_poses[indx][0].total_residue(),
                                              significant_val_threshold=0.1)
                    if (tmp_msd < min_distances[indx]):
                        min_distances[indx]=tmp_msd
                        assignments[indx]=curr_center
                center_ndxs.append(curr_center)
                curr_center=min_distances.argmax()
                if ( len(center_ndxs) >= max_clusters_num):
                    break
        """
    
    else:
        print "I don't know that distance metric!!!: ", align_method
        assert (0==1)
    print "All k-centers Done"
    return center_ndxs, assignments, min_distances

    
def k_center_clusterer_heterogeneusPoses(in_poses=[],
                                         align_method="rmsd",
                                         convergence_dist=0.5,
                                         min_clusters_num=1,
                                         max_clusters_num=9999,
                                         pdb_info_label="DLOOPER"):
    print "Number of elements to cluster=", len(in_poses)
    
    print "Separating by lenght"
    poses_by_len_dic={}
    for indx in xrange(len(in_poses)):
        #print "S: ", in_poses[indx].total_residue()
        if not in_poses[indx].total_residue() in poses_by_len_dic:
            poses_by_len_dic[in_poses[indx].total_residue()]=[[in_poses[indx].clone(),indx]]
        else:
            poses_by_len_dic[in_poses[indx].total_residue()].append([in_poses[indx].clone(),indx])
            
    print "Executing k-centers clustering by total-pose-lenght"
    clustered_poses_centers_ndx=[]
    clustered_poses_cluster_size=[]
    center_ndxs=[]
    for ikey in poses_by_len_dic:
        if (len(poses_by_len_dic[ikey]) > 1): #If there is something to actually cluster
            [center_ndxs,
             assignments, 
             min_distances]=k_center_clusterer_poses(poses_by_len_dic[ikey],
                                      align_method=align_method,
                                     convergence_dist=convergence_dist,
                                     min_clusters_num=min_clusters_num,
                                     max_clusters_num=max_clusters_num,
                                     pdb_info_label="DLOOPER")
            ##print "ALSKJL", center_ndxs, assignments

        else: #otherwise
            center_ndxs=[0]
            assignments=np.asarray([0])
            
        ##print assignments, center_ndxs
        for icenter in  center_ndxs:
            clustered_poses_centers_ndx.append(poses_by_len_dic[ikey][icenter][1])
            clustered_poses_cluster_size.append(len(np.where(assignments==icenter)[0]))
            #Add cluster size???
        #print "AS", clustered_poses_centers_ndx
    print "Finished Heterogeneus Pose Clustering "
    return [np.asarray(clustered_poses_centers_ndx), 
            np.asarray(clustered_poses_cluster_size)]


#Layer Design Definitions

#Make mono-polyAA pose
# replaces the residue at  <resid>  in  <pose>  with  <new_res>  with repacking
def mutate_residues_without_neighbors_consideration( inpose , 
                                                   mutant_positions , 
                                                   mutant_aas ):
    
    
    #if pose.is_fullatom() == False:
    #    IOError( 'mutate_residue only works with fullatom poses' )

    test_pose = inpose.clone()

    #Herein the intrincate method that Rosetta has to do something a simple as mutate a residue, unbelievable
    for mutant_aa_i in xrange(len(mutant_positions)):
        mutant_aa=mutant_aas[mutant_aa_i]
        mutant_position=mutant_positions[mutant_aa_i]
        new_aa_name=pyrosetta.rosetta.core.chemical.name_from_aa(pyrosetta.rosetta.core.chemical.aa_from_oneletter_code( mutant_aa ))
        restype_set=test_pose.residue_type_set_for_pose() #test_pose.residue( mutant_position ).residue_type_set()
        new_res = pyrosetta.rosetta.core.conformation.ResidueFactory.create_residue(rsd_type=restype_set.name_map(new_aa_name), 
                                                                         current_rsd=test_pose.residue(mutant_position),
                                                                         conformation=test_pose.conformation(),
                                                                         preserve_c_beta=False)
        #pyrosetta.rosetta.core.conformation.copy_residue_coordinates_and_rebuild_missing_atoms( source_rsd=test_pose.residue( mutant_position ),
        #                                                                              target_rsd=new_res, 
        #                                                                              conf=test_pose.conformation());
        test_pose.replace_residue( seqpos=mutant_position, 
                                  new_rsd_in=new_res, 
                                  orient_backbone=False );
        
        

    return test_pose

def make_polyAA_pose(in_pose, 
                     targetAA='V'):
    pos_for_mut=range(1, in_pose.total_residue()+1)
    #res_to_mut=[]
    #task_bools_all_true=[]
    #for i in range(in_pose.total_residue()):
    #    res_to_mut.append(targetAA)
    #    aa_bool_tmp=[]
    #    for j in range( 1 , 21 ):
    #        aa_bool_tmp.append(True)
    #    aa_bool = pyrosetta.rosetta.utility.vector1_bool()
    #    for j in range( 1 , 21 ):
    #        aa_bool.append(aa_bool_tmp[j-1]==True)
    #    ##layer_design_task.nonconst_residue_task(i).restrict_absent_canonical_aas( aa_bool )
    #    task_bools_all_true.append(aa_bool)
    connected_pose_mini_desig_polyAA=mutate_residues_without_neighbors_consideration( inpose=in_pose , 
                                                   mutant_positions=pos_for_mut , 
                                                   mutant_aas=(targetAA*len(pos_for_mut)) )
    return connected_pose_mini_desig_polyAA

"""
def mut_pos_to_aa_pose(in_pose, 
                       targetPositions=[],
                       targetAA='V'):
    assert (len(targetPositions)>0)
    assert (len(np.where(targetPositions<=0)[0])==0)
    pos_for_mut=targetPositions #range(1, in_pose.total_residue()+1)
    res_to_mut=[]
    task_bools_all_true=[]
    for i in targetPositions:
        res_to_mut.append(targetAA)
        aa_bool_tmp=[]
        for j in range( 1 , 21 ):
            aa_bool_tmp.append(True)
        aa_bool = pyrosetta.rosetta.utility.vector1_bool()
        for j in range( 1 , 21 ):
            aa_bool.append(aa_bool_tmp[j-1]==True)
        ##layer_design_task.nonconst_residue_task(i).restrict_absent_canonical_aas( aa_bool )
        #task_bools_all_true.append(aa_bool)
    connected_pose_mini_desig_polyAA=mutate_residues_without_neighbors_consideration( inpose=in_pose , 
                                                   mutant_positions=pos_for_mut , 
                                                   mutant_aas=res_to_mut )
    return connected_pose_mini_desig_polyAA
"""

def make_random_sequence_pose(in_pose, 
                              possibleAA=None):
    if ( possibleAA == None ):
        possibleAA=['R','H','K','D','E','S','T','N','Q','A','V','I','L','M','F','Y','W'] #,'P','C','G'
    possibleAA=np.asarray(possibleAA)
    pos_for_mut=range(1, in_pose.total_residue()+1)
    res_to_mut=[]
    #task_bools_all_true=[]
    for i in range(in_pose.total_residue()):
        targetAA=np.random.choice(possibleAA)
        res_to_mut.append(targetAA)
    #    aa_bool_tmp=[]
    #    for j in range( 1 , 21 ):
    #        aa_bool_tmp.append(True)
    #    aa_bool = pyrosetta.rosetta.utility.vector1_bool()
    #    for j in range( 1 , 21 ):
    #        aa_bool.append(aa_bool_tmp[j-1]==True)
    #    ##layer_design_task.nonconst_residue_task(i).restrict_absent_canonical_aas( aa_bool )
    #    task_bools_all_true.append(aa_bool)
    connected_pose_mini_desig_polyAA=mutate_residues_without_neighbors_consideration( inpose=in_pose , 
                                                   mutant_positions=pos_for_mut , 
                                                   mutant_aas=res_to_mut )
    return connected_pose_mini_desig_polyAA

def calc_pose_sasa_total_per_residue_and_atom(in_pose,
                                              probe_radii=2.0,
                                              extra_pertRot_repeats=0):
    
    pertRot_repeats=extra_pertRot_repeats #random_roll_pose always return+1 since the first result is the input itself
    #Sasa algorithm in rosetta is kinda broken, this ensures that we apply extra random-rotations to it in order to minimize the effect by calculating an average.
    if (pertRot_repeats>0):
        test_poses=random_roll_pose(inpose=in_pose,
                                     target_chains=[0], #0 == all the pose
                                     rot_mag=180.0,  #This is irrelevant, just need some numbers here
                                     trans_mag=10.0, #Again, this is irrelevant
                                     num_repeats=pertRot_repeats+1)  #
    else:
        test_poses=[in_pose.clone()]
        
    #print len(test_poses)
    #print test_poses
    #Calculate SASA
    #print "Init"
    total_SASA=0.0
    per_res_sasa=np.zeros((in_pose.total_residue()),float)
    for ipose in test_poses:
        atom_sasa = pyrosetta.rosetta.core.id.AtomID_Map_Real()
        rsd_sasa = pyrosetta.rosetta.utility.vector1_double()
        atom_map = pyrosetta.rosetta.core.id.AtomID_Map_bool()
        for jr in range( 1, ipose.total_residue()+1 ):
                        #Atoms 1 to 5
                        rsd = ipose.residue( jr );
                        for kandx in range (1, rsd.natoms()+1):
                                        atom = pyrosetta.rosetta.core.id.AtomID(kandx, jr)
                                        atom_map.set( atom, True )
                                        ###print atom, atom_map.get(atom)
        
        
        tmp_sasa=pyrosetta.rosetta.core.scoring.calc_per_atom_sasa(ipose,
                                                        atom_sasa,
                                                        rsd_sasa,
                                                        probe_radii,
                                                        False,
                                                        atom_map)
        #print total_SASA
        total_SASA+=tmp_sasa
        #print "+ ", total_SASA, tmp_sasa

        
        for jndx in xrange(len(rsd_sasa)):
            per_res_sasa[jndx]+=rsd_sasa[jndx+1]

    #print "End"
    #return [(total_SASA/(len(test_poses)*1.0)), 
    #        (np.asarray(per_res_sasa)/(len(test_poses)*1.0))]
    return [total_SASA, 
            per_res_sasa]

def res_nonBB_min_distance_to(target_pose=None,
                        target_resA=1,
                        target_resBs=[]):
    min_dist_dic={}
    tmp_resB=target_pose.residue(target_resA)
    for ires in target_resBs:
        min_dist_dic[ires]=np.finfo(float).max
        tmp_resA=target_pose.residue(ires)
        for iatom in xrange(5,tmp_resA.natoms()+1):
            #print tmp_resA
            for jatom in xrange(5,tmp_resB.natoms()+1):
                ##print iatom, jatom
                if ((not tmp_resA.atom_is_backbone(iatom)) and 
                    (not tmp_resB.atom_is_backbone(jatom))):
                    dist=(tmp_resA.atom(iatom).xyz()-tmp_resB.atom(jatom).xyz()).norm()
                    if (dist<min_dist_dic[ires]):
                        #print tmp_resA.atom_name(iatom), tmp_resB.atom_name(jatom), dist
                        min_dist_dic[ires]=dist
        ##print ires, min_dist_dic[ires]
    ##print min_dist_dic
    return min_dist_dic


def pack_pose_generic_and_get_packstat(in_pose,
                                       packaas=["A","L","I","V","F"],
                                       num_packstat_calcs=10,
                                       scorefxn=None):
    packstatfilter=pyrosetta.rosetta.protocols.simple_filters.PackStatFilter()
    pack_task={}
    for ires in  xrange(1,in_pose.total_residue()+1):
        pack_task[ires]=packaas
    packed_pose=pack_rotamers_fast(in_pose=in_pose,
                                   pack_task=pack_task,
                                   pack_scorefxn=scorefxn)
    pack_values=[]
    for i in xrange(num_packstat_calcs):
        pack_values.append(packstatfilter.compute(packed_pose))

    return np.average(pack_values), packed_pose

def pack_pose_generic(in_pose,
                       packaas=["A","L","I","V","F"],
                       scorefxn=None):
    packstatfilter=pyrosetta.rosetta.protocols.simple_filters.PackStatFilter()
    pack_task={}
    for ires in  xrange(1,in_pose.total_residue()+1):
        pack_task[ires]=packaas
    packed_pose=pack_rotamers_fast(in_pose=in_pose,
                                   pack_task=pack_task,
                                   pack_scorefxn=scorefxn)

    return packed_pose

def pack_rotamers_fast(in_pose,
                   pack_task={},
                   pack_scorefxn=None
                  ):
    test_pose = in_pose.clone() #Copy constructor
    pack_scorefxn(test_pose) 
    ##print("DEBUG", test_pose.total_residue())
    task = pyrosetta.standard_packer_task( test_pose )
    for ipos in xrange(1,in_pose.total_residue()+1):
        raa_bool = pyrosetta.rosetta.utility.vector1_bool(20)
        for jndx in range( 20 ):
                raa_bool[jndx+1]=False
        if ipos in pack_task:
            for jaa in pack_task[ipos]:
                mutant_aa = pyrosetta.rosetta.core.chemical.aa_from_oneletter_code( jaa )
                raa_bool[int(mutant_aa)]=True
        task.nonconst_residue_task( ipos
                ).restrict_absent_canonical_aas( raa_bool )
    
    #Pack
    pyrosetta.rosetta.core.pack.pack_rotamers(test_pose,
                                                pack_scorefxn,
                                                task)
    return test_pose

def rodrigues_rotation(a, 
                       angle, 
                       verbose=False):
    #from math import cos, sin, acos, asin, sqrt, pi, radians, degrees
    "use Rodrigues' rotation formula to get rotation matrix"
    a = np.array(a, dtype=float)
    a /= np.sqrt(np.inner(a, a)) # make unit vector
    #assert np.abs(np.sin(angle) - np.sin(np.arccos(np.cos(angle)))) < 1e-6
    #if verbose:
    #    print "rotation angle:", degrees(angle)
    #    print "rotation axis:", a
    omega = np.array([[   0., -a[2],  a[1]],
                   [ a[2],    0., -a[0]],
                   [-a[1],  a[0],    0.]], dtype=float)
    rm = (np.identity(3) + omega * np.sin(angle)
                            + np.dot(omega, omega) * (1 - np.cos(angle)))
    #if verbose:
    #    print "rotation matrix:\n", rm
    return rm

def random_roll_pose(inpose=None,
                     target_chains=[0], #Rosetta_numbers
                     rot_mag=1.0,
                     trans_mag=1.0,
                     num_repeats=3,
                     prepend_original_pose=True):
    
    assert(num_repeats>0)

    #By definition here chain = 0 == all pose
    start_res={0:1}
    stop_res={0:inpose.total_residue()}
    coor_array_original = {}
    com = {}
    for target_chain in target_chains:
        assert(target_chain<=inpose.conformation().num_chains())
        if (target_chain>0):
            start_res[target_chain] = inpose.conformation().chain_begin(target_chain);
            stop_res[target_chain] = inpose.conformation().chain_end(target_chain);


        #Extract original coordinates
        coor_array=[]
        for ires in xrange(start_res[target_chain], (stop_res[target_chain]+1)):
            for jatom in xrange(1, (inpose.residue(ires).natoms()+1)):
                tmp_xyz=np.zeros(3,float)
                for kdim in range(0,3): 
                    tmp_xyz[kdim] = inpose.residue(ires).xyz(jatom)[kdim]
                coor_array.append(tmp_xyz[:])
        coor_array_original[target_chain]=np.asarray(coor_array)
        com[target_chain]=coor_array_original[target_chain].mean(axis=0)

    #result_poses=[]
    #result_poses.append(inpose.clone())
    inpose_copy=inpose.clone()
    for irepeat in xrange(num_repeats-int(prepend_original_pose)):
        if ((irepeat<=0) and prepend_original_pose):
            yield inpose.clone()    
        else:
            for target_chain in target_chains:
                #Restore coordinates
                coor_array=np.copy(coor_array_original[target_chain])

                tVec=np.asarray([0.0,0.0,0.0])
                for iaxis in xrange(3):
                    tVec[iaxis] = pyrosetta.rosetta.numeric.random.gaussian() * trans_mag;

                #Calculate rotation
                angle=pyrosetta.rosetta.numeric.random.gaussian()*rot_mag
                rMat1 = rodrigues_rotation([pyrosetta.rosetta.numeric.random.gaussian(),
                                            pyrosetta.rosetta.numeric.random.gaussian(),
                                            pyrosetta.rosetta.numeric.random.gaussian()],
                                            np.radians(angle))

                ##print rMat1

                #print tVec
                #Apply -T,R,+T
                coor_array=coor_array-com[target_chain]
                coor_array=np.dot(coor_array,rMat1)
                coor_array=coor_array+com[target_chain]
                coor_array=coor_array+tVec

                #Apply_coordinates_to_pose
                tmp_cntr=0 #It is easier just to keep a counter
                for ires in xrange(start_res[target_chain], 
                                   (stop_res[target_chain]+1)):
                    for jatom in xrange(1, (inpose.residue(ires).natoms()+1)):
                        tmp_xyz=pyrosetta.rosetta.numeric.xyzVector_double_t(coor_array[tmp_cntr][0],
                                                                             coor_array[tmp_cntr][1],
                                                                             coor_array[tmp_cntr][2])
                        inpose_copy.set_xyz(pyrosetta.rosetta.core.id.AtomID(jatom, ires),
                                            tmp_xyz)
                        tmp_cntr+=1

            #Append to results
            yield inpose_copy.clone()
            #result_poses.append(inpose_copy.clone())
   
    #yield result_poses


def get_residues_per_layer(in_pose,
                           test_aa_large,
                           test_aa_muta_contact,
                           b_use_global_context,
                           in_global_context_pose,
                           layer_sasa_difinitions,
                           test_aa_large_ref_sasa,
                           layer_contact_distance_difinitions,
                           b_debug=False):
    if b_debug:
        print "--Calculating layers with definitions", layer_sasa_difinitions, "(nothing is your usual subway sandwich!)"
    tmp_mut_pose_largeAA=make_polyAA_pose(in_pose,
                                                test_aa_large)
    tmp_pose_size=tmp_mut_pose_largeAA.total_residue()
    [total_SASA, 
    per_res_sasa]=calc_pose_sasa_total_per_residue_and_atom(tmp_mut_pose_largeAA,
                                                                    probe_radii=1.7)
    #print per_res_sasa
    ##tmp_mut_pose_largeAA.dump_pdb("%s/test_layers/testF.pdb"%out_file_dir)
    by_layer_dic={}
    by_layer_dic["core"]=[]
    by_layer_dic["limbo"]=[]
    by_layer_dic["surface"]=[]
    by_layer_dic["interface"]=[]
    if b_use_global_context:
        #outpdbname+="_c%s"%in_context_pose_basename[0:-4]
        #result_pose=solution_poses[isol][jndx][0].clone()
        tmp_mut_pose_largeAA_wContext=tmp_mut_pose_largeAA.clone()
        tmp_mut_pose_largeAA_wContext.append_pose_by_jump(in_global_context_pose.clone(), 1)
        [total_wContext_SASA, 
        per_res_sasa_wContext]=calc_pose_sasa_total_per_residue_and_atom(tmp_mut_pose_largeAA_wContext,
                                                                         probe_radii=1.7)
        tmp_delta_sasa=np.asarray(per_res_sasa[:tmp_pose_size])-np.asarray(per_res_sasa_wContext[:tmp_pose_size])
        interface_ndxs=np.where(tmp_delta_sasa>=layer_sasa_difinitions["interface"])[0]
        #tmp_mut_pose_largeAA_wContext.dump_pdb("%s/test_layers/testF_wC.pdb"%out_file_dir)
        for lres in interface_ndxs:
            by_layer_dic["interface"].append(lres)

    ##print np.asarray(per_res_sasa).max()
    per_res_sasa=np.asarray(per_res_sasa)/test_aa_large_ref_sasa*100.0

    if b_debug:
        print "--Using pixie dust to detect core-neigbors"
    tmp_test_sasa=[]
    for lres in xrange(len(per_res_sasa)):
        ##if (lres not in by_layer_dic["interface"]):
            if (per_res_sasa[lres] <= layer_sasa_difinitions["core"]):
                by_layer_dic["core"].append(lres)
            elif (per_res_sasa[lres] <= layer_sasa_difinitions["limbo"]):
                ##print "LIMBO", lres+1
                by_layer_dic["limbo"].append(lres)

    tmp_mut_pose_smallAA=make_polyAA_pose(in_pose,
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
                #print "TEST:", lres+1, min_dist
                if (min_dist <= layer_contact_distance_difinitions["limbo"]):
                    if (per_res_sasa[lres] <layer_sasa_difinitions["surf"]):
                        by_layer_dic["limbo"].append(lres)
                    else:
                        by_layer_dic["surface"].append(lres)

                else: #Otherwise is a surface residue
                    #print "SURF",lres+1, min_dist
                    by_layer_dic["surface"].append(lres)
            else: #Otherwise there is no core and therefore no-limbo residues
                by_layer_dic["surface"].append(lres)
    
    return(by_layer_dic)


def get_packstat_nrepeats(in_pose,
                         nrepeats=5):
    result_array=[]
    for i in xrange(nrepeats):
        result_array.append(pyrosetta.rosetta.core.scoring.packstat.compute_packing_score(in_pose))
    return np.average(result_array)



