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

	with open(options.pdb) as readfile:
		lines = readfile.readlines()
		first_atom = 0
		done = False
		for line in lines:
			first_atom += 1
			if done:
				break
			items = line.split()
			for key,value in amino_acids.iteritems():
				if ("ATOM" in items) & (value in items):
					done = True
		chain = lines[first_atom-2].split()[4]

	if options.lig:
		is_selected = False
		print "Extracting ligand coordinates..."
		all_ligands = options.record_type.split()
		with open(options.pdb) as readfile:
			lines = readfile.readlines()
			for line in lines: 
				items = line.split()
				if len(items[0]) > len(ligand):
					serial_num = items[0][len(ligand):]
					items[0] = items[0][:(len(items[0]) - len(ligand)+1)]
					items.insert(1,serial_num)
				if ligand not in items:
					continue
				format_errors = verify_line(items)
				for current_ligand in all_ligands:
					if (format_errors == 0) & (items[3] == current_ligand):
						temp = items[6:9]
						activesite.append(temp)
					else:
						if ((format_errors == 1) | (format_errors == -1)) & (items[3] == current_ligand):
							print "Unexpected line format in the following line: "
							print line
							print "Please reformat this line so that the XYZ coordinates of the ligand are in columns 6-8"
						if (items[2][len(items[2])-3:] == current_ligand) & ((format_errors == 2) | (format_errors == 3)):
							if format_errors == 2:
								temp = items[5:8]
							if format_errors == 3:
								temp = items[4:7]
							activesite.append(temp)
						if (format_errors == 4) & (items[3] == current_ligand):
							temp = items[5:8]
							activesite.append(temp)
		if len(activesite) > 0:
			is_selected = True
		if is_selected:
			print "Ligand coordinates successfully selected!"
		else:
			print "Unable to extract ligand coordinates from the input pdb file. Please check if the provided ligand ID('s) are correct."
			return
	else:
		message = False
		with open(options.asite) as readfile:
			for line in readfile:
				if not message:
					print "Opening file containing active site coordinates..."
					message = True
				items = line.split('\r')
				for item in items:
					to_append = item.split(',')
					temp = []
					for a in to_append:
						m = a[a.find('"')+1:a.rfind('"')]
						temp.append(m)
					if temp[-1] == '':
						temp = temp[:-1]
					if len(temp) != 0:
						activesite.append(temp)
			if message:
				print "Coordinates successfully obtained!"
			else:
				print "Unable to open and extract active site coordinates from provided file. Please ensure that the file containes 1 set of 3D coordinates in the first three columns of the sheet in the form 'X Y Z'"
				return

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

	cleaned_pdb = pdb[:len(pdb)-4]+'_'+chain+'.pdb'

	#Extracting Calpha backbone coordinates from cleaned pdb

	print "Extracting Calpha backbone..."

	with open(cleaned_pdb) as trunc_file:
		for line in trunc_file:
			items = line.split()
			if items[0] == 'ATOM':
				if (items[2] == 'CA'):
					position = [items[5]]
					position.append(items[3])
					position.extend(items[6:9])
					backbone.append(position)

	#Using the Average Degree Filter to calculate the Contact Number for each residue

	contact_nums = contact_number(cleaned_pdb,options.rpath, oscompiler, backbone)

	#Calculate the minimum distance to the active site for each residue

	dist_asite = distance_to_active_site(backbone,activesite)

	#Concatenate results to one list

	residue = 0
	match = True
	for res in dist_asite:
		temp = []
		current_contact_num_id = contact_nums[residue][0]
		current_contact_num = contact_nums[residue][1]
		if (res[1] == current_contact_num_id):
			temp.append(res[0])
			temp.append(res[1])
			temp.append(res[2])
			info.append(temp)
	if match == False:
		print "For missing contact numbers, refer to the file entitled 'contact_numbers.txt' generated in the current directory."

	#Write resulting list to a CSV file

	if len(info) > 4:
		print "Generating CSV file with results..."
		with open(options.outfilename, 'w') as newfile:
			write = csv.wrtier(newfile,quoting=csv.QUOTE_ALL)
			write.writerows(info)
		print "File generated and placed in current directory! Name: "+"'"+options.outfilename+"'"
	else:
		print "Unable to calculate the distance to active site. Please check if the columns in the input file are separated by spaces and no other punctuation."
	with open(options.outfilename) as thefile:
		if sum(1 for line in thefile) == 0:
			print "Unable to write results to a csv file. Printing resulting list instead."
			print info

	try:
		os.remove(cleaned_pdb)
		os.remove(cleaned_pdb[:cleaned_pdb.find('.')]+'.fasta')
		os.remove(cleaned_pdb[:len(cleaned_pdb)-4]+"_0001.pdb")
	except OSError as e:
		print e, "Can't locate cleaned file that was generated. Aborting..."
		return

	return options.outfilename

def verify_line(line):
	#searches for erros in formatting the pdb file
	value = 0
	length = len(line)
	if length < 12:
		if len(line[2]) > 3:
			value = 2
		elif len(line[3]) > 4:
			value += 1
		elif len(line[4]) > 4:
			value = 4
		elif len(line[9]) > 4:
			value = 0
		else:
			value = -1
		return value
	else:
		return value

def contact_number(cleaned_pdb, rpath,oscompiler,backbone):
	print "Calculating the contact number for each residue..."
	sub_call = rpath+'main/source/bin/rosetta_scripts.'+oscompiler+'release -database '+rpath+'main/database/'+' -parser:protocol '+'../contact_num.xml -s '+cleaned_pdb+' -ignore_unrecognized_res -ex1 -ex2 -use_input_sc -no_his_his'
def euclidean(list1, list2)
	
def main():
	if options.path == "GenerateMutations":
		print "Running Generate Mutations..."
		generate_mutations()
	if options.path == 'CalculateDistAndCN':
		print "Running Calculate Distance and Contact Number..."
		enhance_enzyme_solubility()

if __name__ == '__main__':
	main()