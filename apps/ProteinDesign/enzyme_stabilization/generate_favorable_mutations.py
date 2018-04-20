#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   ~/generate_favorable_mutations.py
## @brief  calculating the distance to the active site and contact number for each residue in an enzyme 
## @author Raisa Noshin

import os, math, subprocess, csv

import argparse

from rosetta import *

from pyrosetta import *

init()

parser = argparse.ArgumentParser(prog='Generate Favorable Mutations', description="Use a residue's distance to active site, contact number, and PSSM score to isolate favorable mutations that enhance solubility.")
subparsers = parser.add_subparsers(help='Usage: python generate_favorable_mutations.py {RunMode} -flags', dest='path')

parser_one = subparsers.add_parser('GenerateMutations', help='Generates favorable mutations in the event that the PSSM and distance to active site and contact number are already calculated and stored in specified files.')
parser_one.add_argument('-p','--score_threshold',action='store',default=0,dest='pssm_threshold',help="Indicate the threshold above which PSSM scores will be favorable. Default: 0")
parser_one.add_argument('-d','--distance_threshold',action='store',default=15.0,dest='dist_threshold',help="Indicate the Calpha distance to the active site threshold above which will enhance favorablility. Default: 15.0")
parser_one.add_argument('-c','--contact_threshold',action='store',default=16,dest='contact_threshold',help="Indicate the contact number below which enhance favoribility. Default: 16")
parser_one.add_argument('-o','--outfile',action='store',dest='outfilename',default='',help="Specify the desired name for the final file containing all favorable mutations. Default: 'favorable_mutations.csv'")
parser_one.add_argument('-m','--pssm',action='store',required=True,dest='pssm_filename',help='Provide the file name for the csv file containing the PSSM scores formatted similar to the output of the script written to calculate PSSM score by J. Klesmith.')
parser_one.add_argument('-n','--distance_and_contact',required=True,action='store',dest='dist_and_contact_filename',help='Provide the file name for the csv file containing the distance to the active site and contact number for each residue formatted into four columns with the following information [Position, WT_Residue, Distance to Active Stie, Contact Number]')

parser_two = subparsers.add_parser('CalculateDistAndCN',help='Runs the script written to calculate the distance to the active site and the contact number of each residue in a protein.')
parser_two.add_argument('-o','--outfile',action='store',default='results.csv',dest='outfilename',help="Specify the desired name for the output file containing the distance to the active site and contact number (Default: 'results.csv'): ")

parser_three = subparsers.add_parser('DistCNinBatch',help='Generate a batch of distances to active site and contact numbers.')
parser_three.add_argument('-n','--number_of_models',action='store',dest='num_models',help='Enter the number of models you seek to process.')
parser_three.add_argument('-e','--exact_struct',action='store',dest='exact',help='Enter the name of the pdb file containing the exact structure of the protein.')

parser_four = subparsers.add_parser('MutationsinBatch',help='Generate a batch of favorable mutations')
parser_four.add_argument('-s','--s_threshold',action='store',default=0,dest='score_threshold',help="Indicate the threshold above which PSSM scores will be favorable. Default: 0")
parser_four.add_argument('-w','--d_threshold',action='store',default=15.0,dest='distance_threshold',help="Indicate the Calpha distance to the active site threshold above which will enhance favorablility. Default: 15.0")
parser_four.add_argument('-l','--c_threshold',action='store',default=16,dest='c_threshold',help="Indicate the contact number below which enhance favoribility. Default: 16")
options = parser.parse_args()

amino_acids = {'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL','U':'SEL','O':'PYL'}

def generate_favorable_mutations(PSSM,dist_cn):
	if options.path == 'GenerateMutations':
		try:
			pssm_threshold = int(options.pssm_threshold)
			dist_threshold = float(options.dist_threshold)
			contact_threshold = int(options.contact_threshold)
		except (ValueError, IndexError) as e:
			print e, "One of the following three input thresholds is unrecognizable: ", options.pssm_threshold,"  ",options.dist_threshold,"  ",options.contact_threshold
			return
	if options.path == 'MutationsinBatch':
		try:
			pssm_threshold = int(options.score_threshold)
			dist_threshold = float(options.distance_threshold)
			contact_threshold = int(options.c_threshold)
		except (ValueError, IndexError) as e:
			print e, "One of the following three input thresholds is unrecognizable: ", options.score_threshold,"  ",options.distance_threshold,"  ",options.c_threshold
			return

	favorable_mutations = [['Position','WT','Favorable_Mutation','Distance_to_Active_Site','Contact_Number','PSSM_Score','Designation']] 
	all_mutations = [['Position','WT','Mutation','Distance_to_Active_Site','Contact_Number','PSSM_Score','Designation']]
	contact_and_distance = []
	rows = []
	columns = []

	#Extracting all rows from PSSM file

	print "Opening PSSM file..."
	with open(PSSM,'rU') as input:
		reader = csv.reader(input)
		data = list(list(rec) for rec in csv.reader(input,delimiter=','))
	for i in range(23):
		if i == 2:
			pass
		else:
			rows.append(data[i])
	if len(rows) == 0:
		print "Unable to open PSSM file. Please ensure that the file name provided includes the .csv extension."
		return

	#Transposing rows to get a list of columns from the PSSM file 

	unedited_columns = map(list,zip(*rows))

	#Editing columns list to remove items with 'None' listed instead of a score

	columns=[]
	for column in unedited_columns:
		counter = 0
		for row in column:
			row = row.rstrip()
			if row == 'None':
				counter += 1
		if counter == 0:
			columns.append(column)

	#Generating a list of just the PSSM scores of each residue

	scores = []
	for column in range(len(columns)):
		scores.append(columns[column][2:])

	if len(scores) == 0:
		print "Unable to extract scores from PSSM file. Please ensure that the file is a csv file and the .csv extension is included in the filename provided. "
		return
	else:
		print "PSSM scores successfully extracted! "

	#Extracting the distance to active site and contact number for each residue into a list of lists entitled values

	print "Opening distance to active site and contact number file..."
	with open(dist_cn) as input:
		rows = csv.reader(input)
		contact_and_distance = list(rows)
	mutations = scores[0]
	scores = scores[1:]
	values = contact_and_distance[1:]
	if len(values) == 0:
		print "Unable to extract distances and contact numbers from file. Please ensure that the file is a csv file and the .csv extension is included in the filename provided. "
		return
	else:
		print "Distance and contact number data successfully extracted!"
	if dist_cn.split('_')[0] != '2dfb':
		columns = columns[1:]
		for j in columns:
			for i in range(len(values)):
				for key, value in amino_acids.iteritems():
					if values[i][1] == value:
						WT = key
				if (int(values[i][0]) == int(j[0])) & (WT != j[1]):
					values.insert(i,[j[0],amino_acids[j[1]],0,50,"NOTE: Residue deleted from model"])
					for n in range(i+1,len(values)):
						p = int(values[n][0])
						values[n][0] = p + 1
						values[n][0] = str(values[n][0])
				else:
					continue
	else:
		pass

	#Filtering favorable mutations

	print "Detecting favorable mutations..."
	for position,score in enumerate(scores):
		if position > len(values):
			print "Your models are missing one or more residues. Only residues that are present have been evaluated."
			break
		current_position = values[position][0]
		current_wildtype = values[position][1]
		current_distance = float(values[position][2])
		current_contactnum = int(values[position][3])
		for pointer,current_score in enumerate(score):
			current_mutation = mutations[pointer]
			for key, value in amino_acids.iteritems():
				if current_wildtype == value:
					current_wildtype_small = key
			temp1 = [current_position,current_wildtype,current_mutation,current_distance,current_contactnum,current_score,current_wildtype_small+current_position+current_mutation]
			all_mutations.append(temp1)
			if (int(current_score) >= pssm_threshold) & (current_distance >= dist_threshold) & (current_contactnum <= contact_threshold):
				temp = [current_position,current_wildtype,current_mutation,current_distance,current_contactnum,current_score,current_wildtype_small+current_position+current_mutation]
				favorable_mutations.append(temp)
				print "Favorable mutation detected at position "+current_position+" from "+current_wildtype+" --> "+current_mutation

	#Writing all possible mutations to a csv file

	if len(all_mutations) > 7:
		print "Generating csv file with all possible mutations..."
		filename = dist_cn[:dist_cn.find('.')]+"_all_mutations.csv"
		with open(filename,'w') as newfile:
			write = csv.writer(newfile,quoting=csv.QUOTE_ALL)
			write.writerows(all_mutations)
		print "File containing all possible mutations generated and placed in project folder! Name: all_mutations.csv"
		with open(filename) as thefile:
			if sum(1 for line in thefile) == 0:
				print "Unable to write results to a csv file. Printing all possible mutations instead."
				print all_mutations

	#Writing favorable mutations to a csv file 

	outfilename = dist_cn[:dist_cn.find('.')]+"_favorable_mutations.csv"
	
	if len(favorable_mutations) > 7:
		print "All favorable mutations generated!"
		print "Generating csv file with results..."
		with open(outfilename, 'w') as newfile:
			write = csv.writer(newfile,quoting=csv.QUOTE_ALL)
			write.writerows(favorable_mutations)
		print "File containing all favorable mutations generated and placed in current directory! Name: "+outfilename
		with open(outfilename) as thefile:
			if sum(1 for line in thefile) == 0:
				print "Unable to write results to a csv file. Printing favorable mutations instead."
				print favorable_mutations
	else:
		print "No favorable mutations detected with the current provided thresholds. Please check if these values are accurate, and that the columns in the input file are separated by spaces and no other punctuation."
		return

	return outfilename

def generate_favorable_mutations_batch():
	
	pssm = raw_input("Enter the name of the file containing the PSSM scores for the decoy set in the form 'filename.csv': ")
	print ""
	true = raw_input("Enter the name of the file containing the distance to active site and contact number of the exact structure in the form 'filename.csv': ")
	dist_cn = [true]
	print ""
	with open(true[:true.find("_")]+"_decoyset_dist_cn_names.txt") as thefile:
		for line in thefile:
			temp = line.splitlines()
			dist_cn.extend(temp)	
	all_files = []
	for model in dist_cn:
		outputname = generate_favorable_mutations(pssm,model)
		all_files.append(outputname)
	all_files = all_files[1:]
	with open(true[:true.find('_')]+"_decoyset_favorable_mutation_names.txt","w") as thefile:
		for item in all_files:
			thefile.write("%s\n" % item)

def enzyme_solubility_factors():
	#Script for calculating the Calpha distance to active site as well as the contact number for each residue in the backbone of an enzyme
	marker = 0
	while marker == 0:
		pdb = raw_input("Enter the name of the file containing the desired pdb structure in the form of 'filename.pdb': ")
		try:
			f = open("../"+pdb)
			f.close()
			marker = 1
		except (OSError, IOError) as e:
			print e, "Unable to find input file. Please ensure that you spelled the name correctly and that it is present in the current directory."
			marker = 0
	oscompiler = raw_input("Are you using a Mac derivative or Linux operating system? ")
	lig = raw_input("Enter whether the input pdb structure contains information for one or more ligands present in the active site (True or False): ")
	if lig == 'False':
		asite = raw_input("Enter the name of the file containing the activesite coordinates in a csv file with each 3D coordinate on a new line: ")
	else:
		asite = ''
		record_type = raw_input("Enter the record name in the first column designating the ligand, or press enter to accept the default ('HETATM'): ")
	if record_type == '':
		record_type = 'HETATM'
	if (oscompiler == 'Mac') | (oscompiler == 'MAC') | (oscompiler == 'mac'):
		oscompiler = 'macosclang'
	elif (oscompiler == 'Linux') | (oscompiler == 'linux') | (oscompiler == 'LINUX'):
		oscompiler = 'linuxgcc'
	else:
		print "Invalid operator chosen. Please check your spelling. Accepted answers: Mac, MAC, mac, Linux, linux, LINUX"
	rpath = '../../Rosetta/'
	activesite = []
	backbone = []
	info = [['Position','WT_Residue','Min_Distance_to_Active_Site','Contact_Number']]
	contact_nums = []
	dist_asite = []
	all_ligands = []

	#Extracting the ligand(s) from either the input pdb file or the provided coordinates 

	if lig == 'True':
		number_of_ligands = input("Enter the number of ligands present in the active site of this input structure: ")

		for ligand in range(1,number_of_ligands+1):
			current_ligand = raw_input("Enter the ligand ID of ligand number " + str(ligand) + ": ")
			all_ligands.append(current_ligand)
		is_selected = False
		print "Extracting ligand coordinates..."
		with open('../'+pdb) as readfile:
			lines = readfile.readlines()
			first_atom = 0
			done = False
			for line in lines:
				first_atom += 1
				if done:
					break
				items = line.split()
				for key, value in amino_acids.iteritems():
					if ("ATOM" in items) & (value in items):
						done = True
			chain = lines[first_atom-2].split()[4]
			for line in lines:
				items = line.split()
				if len(items[0]) > len(record_type):
					serial_num = items[0][len(record_type):]
					items[0] = items[0][:(len(items[0]) - len(record_type)+1)]
					items.insert(1,serial_num)
				if record_type not in items:
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
			print "Unable to extract ligand coordinates from input pdb file. Please check if the provided ligand ID('s) are correct."
			return

	if asite != '':
		with open('../'+pdb) as readfile:
			lines = readfile.readlines()
			first_atom = 0
			done = False
			for line in lines:
				first_atom += 1
				if done:
					break
				items = line.split()
				for key, value in amino_acids.iteritems():
					if ("ATOM" in items) & (value in items):
						done = True
			chain = lines[first_atom-2].split()[4]
		message = False
		with open('../'+asite) as readfile:
			for line in readfile:
				if not message:
					print "Opening file containing active site coordinates..."
					message = True
				items = line.split('\r')
				for item in items:
					to_append = item.split(',')
					activesite.append(to_append)
			if message:
				print "Coordinates successfully obtained!"
			else:
				print "Unable to open and extract active site coordinates form file provided. Please ensure that the file contains 1 set of 3D coordinates separated by spaces on each line in the form of 'X Y Z'"
				return

	#Cleaning pdb of all extraneous information and renumbering 

	call_clean = rpath+'tools/protein_tools/scripts/clean_pdb.py ../'+pdb+' '+chain
	try:
		print "Cleaning and renumbering input pdb file..."
		status = os.system(call_clean)
		subprocess.call(call_clean,shell=True)
	except (OSError, ValueError) as e:
		print e, ": Unable to renumber and clean the pdb file. Please check if your path to Rosetta is correct, or that your pdb file is in the current directory."
		return
	except subprocess.CalledProcessError as f:
		print f.output
		return
	else:
		print "Pdb successfully renumbered!"

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

	contact_nums = contact_number(cleaned_pdb,rpath,oscompiler,backbone)

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
			temp.append(current_contact_num)
			info.append(temp)
			residue += 1
		else:
			match = False
			print "Unable to display contact number for " + res[1] + " at position " + res[0]
			temp.append(res[0])
			temp.append(res[1])
			temp.append(res[2])
			info.append(temp)
	if match == False:
		print "For missing contact numbers, refer to the file entitled 'contact_numbers.txt' generated in the current directory."
	
	#Write resulting list to a CSV file

	if len(info) > 4:
		print "Generating csv file with results..."
		with open(options.outfilename, 'w') as newfile:
			write = csv.writer(newfile,quoting=csv.QUOTE_ALL)
			write.writerows(info)
		print "File generated and placed in current directory! Name: "+"'"+options.outfilename+"'"
	else:
		print "Unable to calculate the distance to active site. Please check if the columns in the input file are separated by spaces and no other punctuation."
	with open(options.outfilename) as thefile:
		if sum(1 for line in thefile) == 0:
			print "Unable to write results to a csv file. Printing resulting list instead."
			print info

def enzyme_solubility_factors_no_input(pdb,oscompiler,lig,numlig,ligIDs,asite,record_type,outfilename):
	#Script for calculating the Calpha distance to active site as well as the contact number for each residue in the backbone of an enzyme

	rpath = '../../Rosetta/'
	activesite = []
	backbone = []
	info = [['Position','WT_Residue','Min_Distance_to_Active_Site','Contact_Number']]
	contact_nums = []
	dist_asite = []
	all_ligands = ligIDs
	position_before_renumbering = []
	#Extracting the ligand(s) from either the input pdb file or the provided coordinates 

	if lig == 'True':
		is_selected = False
		print "Extracting ligand coordinates..."
		with open('../'+pdb) as readfile:
			lines = readfile.readlines()
			first_atom = 0
			done = False
			for line in lines:
				first_atom += 1
				if done:
					break
				items = line.split()
				for key, value in amino_acids.iteritems():
					if ("ATOM" in items) & (value in items):
						done = True
			chain = lines[first_atom-2].split()[4]
			for line in lines:
				items = line.split()
				if len(items[0]) > len(record_type):
					serial_num = items[0][len(record_type):]
					items[0] = items[0][:(len(items[0]) - len(record_type)+1)]
					items.insert(1,serial_num)
				if record_type not in items:
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
			print "Unable to extract ligand coordinates from input pdb file. Please check if the provided ligand ID('s) are correct."
			return

	if asite != '':
		with open('../'+pdb) as readfile:
			lines = readfile.readlines()
			first_atom = 0
			done = False
			for line in lines:
				first_atom += 1
				if done:
					break
				items = line.split()
				for key, value in amino_acids.iteritems():
					if ("ATOM" in items) & (value in items):
						done = True
			chain = lines[first_atom-2].split()[4]
		message = False
		with open('../'+asite) as readfile:
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
				print "Unable to open and extract active site coordinates form file provided. Please ensure that the file contains 1 set of 3D coordinates separated by spaces on each line in the form of 'X Y Z'"
				return

	#Cleaning pdb of all extraneous information and renumbering 

	call_clean = rpath+'tools/protein_tools/scripts/clean_pdb.py ../'+pdb+' '+chain
	try:
		print "Cleaning and renumbering input pdb file..."
		status = os.system(call_clean)
		subprocess.call(call_clean,shell=True)
		print "Pdb successfully renumbered!"
	except (OSError, ValueError) as e:
		print e, ": Unable to renumber and clean the pdb file. Please check if your path to Rosetta is correct, or that your pdb file is in the current directory."
		return
	except subprocess.CalledProcessError as f:
		print f.output
		

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

	contact_nums = contact_number(cleaned_pdb,rpath,oscompiler,backbone)

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
			temp.append(current_contact_num)
			info.append(temp)
			residue += 1
		else:
			match = False
			print "Unable to display contact number for " + res[1] + " at position " + res[0]
			temp.append(res[0])
			temp.append(res[1])
			temp.append(res[2])
			info.append(temp)
	if match == False:
		print "For missing contact numbers, refer to the file entitled 'contact_numbers.txt' generated in the current directory."
	
	#Write resulting list to a CSV file

	if len(info) > 4:
		print "Generating csv file with results..."
		with open(outfilename, 'w') as newfile:
			write = csv.writer(newfile,quoting=csv.QUOTE_ALL)
			write.writerows(info)
		print "File generated and placed in current directory! Name: "+"'"+outfilename+"'"
	else:
		print "Unable to calculate the distance to active site. Please check if the columns in the input file are separated by spaces and no other punctuation."
	with open(outfilename) as thefile:
		if sum(1 for line in thefile) == 0:
			print "Unable to write results to a csv file. Printing resulting list instead."
			print info

	try:
		os.remove(cleaned_pdb)
		os.remove(cleaned_pdb[:cleaned_pdb.find('.')]+".fasta")
		os.remove(cleaned_pdb[:len(cleaned_pdb)-4]+"_0001.pdb")
	except OSError as e:
		print e, "Can't locate cleaned file that was just generated. Aborting..."
		return

	return outfilename

def euclidean(list1,list2):
	#This function calculates the euclidean distance from two lists containing 3-dimensional coordinates
	point1 = []
	point2 = []
	
	for coordinate in list1:
		if coordinate == '':
			continue
		else:
			to_append = coordinate.rstrip()
			point1.append(to_append)
	for coordinate in list2:
		if coordinate == '':
			continue
		else:
			to_append = coordinate.rstrip()
			point2.append(to_append)
	print point1,point2
	print ""
	try:
		x1,y1,z1 = float(point1[0]),float(point1[1]),float(point1[2])
		x2,y2,z2 = float(point2[0]),float(point2[1]),float(point2[2])
	except (ValueError, IndexError) as e:
		print e, "One of the following two coordinates is unrecognizable: ", point1,"    ",point2
		return

	return math.sqrt(math.pow((x2-x1),2)+math.pow((y2-y1),2)+math.pow((z2-z1),2)) 

def distance_to_active_site(backbone,activesite):
	#takes each residue in the backbone and calculates the euclidean distance 
	#between all the points of the ligand, finally selecting the minimum distance as the official distance to active site
	list = []
	for res in backbone:
		temp = res[2:]
		min_dist = 100000000.0
		for coordinate in activesite:
			dist = euclidean(temp,coordinate)
			if dist < min_dist:
				min_dist = dist
		temp1 = [res[0],res[1],min_dist]
		list.append(temp1)
	return list

def contact_number(cleaned_pdb,rpath,oscompiler,backbone):
	#calls the contact_num.xml script written to use the AverageDegree Filter from Rosetta to calculate the contact number
	print "Calculating the contact number for each residue..."
	sub_call = rpath+'main/source/bin/rosetta_scripts.'+oscompiler+'release -database '+rpath+'main/database/'+' -parser:protocol '+'../contact_num.xml -s '+cleaned_pdb+' -ignore_unrecognized_res -ex1 -ex2 -use_input_sc -no_his_his_pairE'+' -flip_HNQ -overwrite > contact_numbers.txt'
	try:
		status = os.system(sub_call)
		subprocess.call(sub_call,shell=True)
		print "Contact Numbers successfully generated!"
	except (OSError, ValueError):
		print e, ": Unable to calculate the contact number. Please check if your path to Rosetta is correct, and if you specified the correct OS and compiler. Also make sure your pdb file, this script, and the file entitled 'contact_num.xml' are within the current directory."
		return
	except subprocess.CalledProcessError as f:
		print f.output
		return
	with open('contact_numbers.txt') as contact:
		lines = []
		output = []
		for line in contact:
			lines.append(line)
		begin = -1
		end = -1
		for element in lines:
			if end == -1:
				end = element.find("END FILTER")
				if end != -1:
					break
			else:
				break
			if begin == -1:
				begin = element.find("BEGIN FILTER")
			else:
				if(isinstance(element.find("Connectivity"),(int,float))):
					line = element.split()
					amino_acid = line[-3][:3]
					temp = [amino_acid,line[-1]]
					output.append(temp)
	if len(output) != len(backbone):
		print "WARNING in processing the contact number: Numbers generated do not fully match the number of backbone residues."
	return output

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

def main():
	if options.path == "GenerateMutations":
		print "Running Generate Mutations..."
		generate_favorable_mutations(options.pssm_filename,options.dist_and_contact_filename)
	if options.path == "CalculateDistAndCN":
		print "Running Calculate Distance and Contact Number..."
		enzyme_solubility_factors()
	if options.path == "DistCNinBatch":
		print "Starting Batch..."
		enzyme_solubility_factors_info()
	if options.path == "MutationsinBatch":
		print "Generating Favorable Mutations in Batch..."
		generate_favorable_mutations_batch()

if __name__ == '__main__':
	main()