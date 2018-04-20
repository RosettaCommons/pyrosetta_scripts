#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   ~/enhance_enzyme_solubility.py
## @brief  calculating the distance to the active site and contact number for each residue in an enzyme 
## @author Raisa Noshin

import os, math, subprocess, csv

import argparse

from rosetta import *

from pyrosetta import *

init()

parser = argparse.ArgumentParser(prog='Enhance Enzyme Solubility', description="Use a residue's distance to active site, contact number, and PSSM score to isolate solubility-enhancing mutations that preserve catalytic activity.")
subparsers = parser.add_subparsers(help='Usage: python enhance_enzyme_solubility.py {RunMode} -flags', dest='path')

parser_one = subparsers.add_parser('GenerateMutations', help='Uses PSSM score, distance to active site, and contact number information to calculate solubilty-enhancing mutations that maintain catalytic activity.')
parser_one.add_argument('-p','--score_threshold',action='store',default=0,dest='pssm_threshold', help="Indicate the threshold above which PSSM scores will be favorable. Default: 0")
parser_one.add_argument('-d','--distance_threshold',action='store',default=15.0,dest='dist_threshold',help="Indicate the Calpha distance to the active site threshold above which will enhance favorablility. Default: 15.0")
parser_one.add_argument('-c','--contact_threshold',action='store',default=16,dest='contact_threshold',help="Indicate the contact number below which enhance favoribility. Default: 16")
parser_one.add_argument('-o','--outfile',action='store',dest='outfilename',default='',help="Specify the desired name for the final file containing all favorable mutations. Default: 'favorable_mutations.csv'")
parser_one.add_argument('-m','--pssm',action='store',required=True,dest='pssm_filename',help='Provide the file name for the csv file containing the PSSM scores formatted similar to the output of the script written to calculate PSSM score by J. Klesmith.')
parser_one.add_argument('-n','--distance_and_contact',required=True,action='store',dest='dist_and_contact_filename',help='Provide the file name for the csv file containing the distance to the active site and contact number for each residue formatted into four columns with the following information [Position, WT_Residue, Distance to Active Stie, Contact Number]')

parser_two = subparsers.add_parser('CalculateDistAndCN',help='Runs the script written to calculate the distance to the active site and the contact number of each residue in a protein.')
parser_two.add_argument('-o','--outfile',action='store',default='results.csv',dest='outfilename',help="Specify the desired name for the output file containing the distance to the active site and contact number (Default: 'results.csv'): ")
parser_two.add_argument('-p','--pdb',action='store',required=True,dest='pdb',help="Enter the path to the pdb file, noting that the last thing should have the form of 'filename.pdb' (REQUIRED)")
parser_two.add_argument('-c','--compiler',action='store',required=True,dest='oscompiler',help='Enter your operating system (i.e. Mac or Linux) (REQUIRED)')
parser_two.add_argument('-l','--lig',action='store_true',dest='lig',help='Use this flag to indicate that a ligand is present within the pdb structure.')
parser_two.add_argument('-h','--recordtype',action='store',default='H',dest='record_type',help='Enter the Atom ID of the ligand. If there are multiple, separate them by a space. Note: Must be specified if -l flag is used')
parser_two.add_argument('-r','--rosetta',action='store',dest='rpath',required=True,help='Enter the path to Rosetta. (REQUIRED)')
parswer_two.add_argument('-a','--asite',action='store',dest='asite',help='Enter the path to the file containing the active site coordinates for the pdb structure in question. Note: This is required if the -l flag is not used.')
options = parser.parse_args()

def enhance_enzyme_solubility():
	#Script for calculating the Calpha distance to active site as well as the contact number for each residue in the backbone of an enzyme
	
	#Check if the given files can be opened
	try:
		f = open(options.pdb)
		f.close()
		if options.asite:
			g = open(options.asite)
			g.close()
	except (OSError, IOError) as e:
		print e
		print "Unable to find input file or the file containing the active site. Please ensure that you specified the right path."
		return

	if options.oscompiler.lower() == 'mac':
		oscompiler = 'macosclang'
	elif options.oscompiler.lower() == 'linux':
		oscompiler = 'linuxgcc'
	else:
		print "Accepted operators are Mac and Linux, any other operator is not supported."
	if options.rpath[-1] != '/'
		options.rpath += '/'

	activesite = []
	backbone = []
	info = [['Position','WT_Residue','Min_Distance_to_Active_Site','Contact_Number']]
	contact_nums = []
	dist_asite = []

	#Extracting the ligand(s) from either the input pdb file or the provided coordinates

	if options.lig:
		all_ligands = options.record_type.split()


	#Cleaning pdb of all extraneous information and renumbering

	call_clean = options.rpath+'tools/protein_tools/scripts/clean_pdb.py '+options.pdb+' '+chain
	try:
		print "Cleaning and renumbering input pdb file..."
		status = os.system(call_clean)
		subprocess.call(call_clean,shell=True)
	except (OSError, ValueError) as e:
		print e
		print "Unable to renumber and clean the pdb file. Please check if your path to Rosetta is correct."
		return
	except subprocess.CalledProcessError as f:
		print f.output
		return 
	else:
		print "PDB successfully renumbered!"

	
def main():
	if options.path == "GenerateMutations":
		print "Running Generate Mutations..."
		generate_mutations()
	if options.path == 'CalculateDistAndCN':
		print "Running Calculate Distance and Contact Number..."
		enhance_enzyme_solubility()

if __name__ == '__main__':
	main()