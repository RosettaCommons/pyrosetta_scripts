#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   ~/enhance_enzyme_solubility.py
## @brief  calculating the distance to the active site and contact number for each residue in an enzyme, and also using filters on PSSM, contact nubmer and distance to active site to generate all favorable mutations for a given enzyme. 
## @author Raisa Noshin

import os, math, subprocess, csv

import argparse

from rosetta import *

from pyrosetta import *

init()

parser = argparse.ArgumentParser(prog='Enhance Enzyme Solubility', description="Use a residue's distance to active site, contact number, and PSSM score to isolate solubility-enhancing mutations that preserve catalytic activity.")
subparsers = parser.add_subparsers(help='Usage: python enhance_enzyme_solubility.py {RunMode} -flags', dest='path')

parser_one = subparsers.add_parser('GenerateMutations', help='Uses PSSM score, distance to active site, and contact number information to calculate solubilty-enhancing mutations that maintain catalytic activity.')
parser_one.add_argument('-s','--score_threshold',action='store',default=0,dest='pssm_threshold', help="Indicate the threshold above which PSSM scores will be favorable. Default: 0")
parser_one.add_argument('-d','--distance_threshold',action='store',default=15.0,dest='dist_threshold',help="Indicate the Calpha distance to the active site threshold above which will enhance favorablility. Default: 15.0")
parser_one.add_argument('-c','--contact_threshold',action='store',default=16,dest='contact_threshold',help="Indicate the contact number below which enhance favoribility. Default: 16")
parser_one.add_argument('-o','--outfile',action='store',dest='outfilename',default='enzyme_0001',help="Specify the desired name for the project. Default: 'enzyme_0001'")
parser_one.add_argument('-m','--pssm',action='store',required=True,dest='pssm_filename',help='Provide the file name for the csv file containing the PSSM scores formatted similar to the output of the script written to calculate PSSM score by J. Klesmith.')
parser_one.add_argument('-n','--distance_and_contact',required=True,action='store',dest='dist_and_contact_filename',help='Provide the file name for the csv file containing the distance to the active site and contact number for each residue formatted into four columns with the following information [Position, WT_Residue, Distance to Active Stie, Contact Number]')
parser_one.add_argument('-p','--possibilities',action='store_true',dest='possibilites',help='Use this flag to indicate that you would like a file generated with all possible mutations and their respective information.')

parser_two = subparsers.add_parser('CalculateDistAndCN',help='Runs the script written to calculate the distance to the active site and the contact number of each residue in a protein.')
parser_two.add_argument('-o','--outfile',action='store',default='results.csv',dest='outfilename',help="Specify the desired name for the output file containing the distance to the active site and contact number (Default: 'results.csv'): ")
parser_two.add_argument('-p','--pdb',action='store',required=True,dest='pdb',help="Enter the path to the pdb file, noting that the last thing should have the form of 'filename.pdb' (REQUIRED)")
parser_two.add_argument('-c','--compiler',action='store',required=True,dest='oscompiler',help='Enter your operating system (i.e. Mac or Linux) (REQUIRED)')
parser_two.add_argument('-l','--lig',action='store_true',dest='lig',help='Use this flag to indicate that a ligand is present within the pdb structure.')
parser_two.add_argument('-h','--recordtype',action='store',default='H',dest='record_type',help='Enter the Atom ID of the ligand. If there are multiple, separate them by a space. Note: Must be specified if -l flag is used')
parser_two.add_argument('-r','--rosetta',action='store',dest='rpath',required=True,help='Enter the path to Rosetta. (REQUIRED)')
parswer_two.add_argument('-a','--asite',action='store',dest='asite',help='Enter the path to the file containing the active site coordinates for the pdb structure in question. Note: This is required if the -l flag is not used.')
options = parser.parse_args()

amino_acids = {'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL','U':'SEL','O':'PYL'}

#Distance to active site: Defined as the minimum Calpha distance of a specified residue to the active site. All active site coordinates provided are used to evaluate the minimum. In optimal conditions, this value is high (>15). 
#Contact Number: Defined as the number of adjacent residues to a specified residue. In optimal conditions, this value is low (<16)
#PSSM: Defined as the evolutionary conservation of a specifed residue. This script runs with the assumption that the format of the file is in accordance with J. Klesmith's PSSM-calculating script (GitHub username: Jklesmith). In optimal conditions, this value is high (>0)
#Above filters obtained from http://www.pnas.org/content/114/9/2265


def dist_contact_filters():
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
	if options.rpath[-1] != '/':
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

def euclidean(list1,list2):
	#calculates the euclidean distance from two lists containing 3-dimensional coordinates
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
	
def distance_to_active_site(backbone, activesite):
	#takes each residue in the backbone and calculates the euclidean distance between all the points of the ligand, finally selecting the minimum distance as the official distance to active site 
	list = []
	for res in backbone:
		temp = res[2:]
		min_dist = 100000000.0
		for coordinate in activesite:
			dist = euclidean(temp, coordinate)
			if dist < min_dist:
				min_dist = dist
		temp1 = [res[0],res[1],min_dist]
		list.append(temp1)
	return list

def enhance_enzyme_solubility():
	PSSM = options.pssm_filename
	dist_cn = options.dist_and_contact_filename
	try:
		pssm_threshold = int(options.pssm_threshold)
		dist_threshold = float(options.dist_threshold)
		contact_threshold = int(options.contact_threshold)
	except (ValueError, IndexError) as e:
		print e, "One of the following three input thresholds is unrecognizable: ", options.pssm_threshold,"  ",options.dist_threshold,"  ",options.contact_threshold
		return

	favorable_mutations = [['Position','WT','Favorable_Mutation','Distance_to_Active_Site','Contact_Number','PSSM_Score','Designation']] 
	all_mutations = [['Position','WT','Mutation','Distance_to_Active_Site','Contact_Number','PSSM_Score','Designation']]
	contact_and_distance = []
	rows = []
	columns = []

	#Extracting all rows from PSSM file, with the exception of a blank 3rd row

	print "Opening PSSM file..."
	with open(PSSM, 'rU') as input:
		reader = csv.reader(input)
		data = list(list(rec) for rec in csv.reader(input,delimiter=','))
	for i in range(23):
		if i == 2:
			pass
		else:
			rows.append(data[i])
	if len(rows) == 0:
		print "Unable to open PSSM file. Please ensure that the file name provided includes the .csv extension. Example: 'filename.csv'"
		return

	#Transposing rows to get a list of columns from the PSSM file

	unedited_columns = map(list,zip(*rows))

	#Editing columns list to remove items with 'None' listed instead of a score

	for column in unedited_columns:
		counter = 0
		for row in column:
			row = row.rstrip()
			if row == "None":
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

	#Extracting the distance to active site and contact number for each residue into a list of lists intitled values

	print "Opening distance to active site and contact number file..."
	with open(dist_cn) as input:
		rows = csv.reader(input)
		contact_and_distance = list(rows)
	mutations = scores[0]
	scores = scores[1:]
	values = contact_and_distance[1:]
	if len(values) == 0:
		print "Unable to extract distances and contact numbers from file. Please ensure that the file is a csv file provided in the format 'filename.csv'"
		return
	else:
		print "Distance and contact number data successfully extracted!"
	columns = columns[1:]
	#below makes sure the residues in the values list (contact numbers and distance to active site) match those in the columns list (from PSSM)
	for j in columns:
		for i in range(len(values)):
			for key, value in amino_acids.iteritems():
				if values[i][1] == value:
					WT = key
				if (int(values[i][0]) == int(j[0])) & (WT != j[1]):
					values.insert(i,[j[0],amino_acids[j[1]],0,50,"Note: Residue deleted from model"])
					for n in range(i+1,len(values)):
						p = int(values[n][0])
						values[n][0] = p + 1
						values[n][0] = str(values[n][0])
				else:
					continue	

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

	if options.possibilities:
		if len(all_mutations) > 7:
			print "Generating csv file with all possible mutations..."
			filename = options.outfilename+"_all_mutations.csv"
			with open(filename, 'w') as newfile:
				write = csv.writer(newfile, quoting=csv.QUOTE_ALL)
				write.writerows(all_mutations)
			with open(filename) as thefile:
				if sum(1 for line in thefile) <= 1:
					print "Unable to write results to a csv file. Printing all possible mutations instead..."
					print all_mutations
				else:
					print "File containing all possible mutations generated and placed in project folder! Name: "+filename
		else:
			print "Error in generating mutations. Please check that your files are formatted correctly and that the paths given for them are correct."
			return

	#Writing favorable mutations to a csv file

	outfilename = options.outfilename
	if len(favorable_mutations) > 7:
		print "All favorable mutations generated!" 
		print "Generating csv file with results..."
		with open(outfilename, 'w') as newfile:
			write = csv.writer(newfile,quoting=csv.QUOTE_ALL)
			write.writerows(favorable_mutations)
		with open(outfilename) as thefile:
			if sum(1 for line in thefile) <= 1:
				print "Unable to write results to a csv file. Printing favorable mutations instead..."
				print favorable_mutations
			else:
				print "File containing favorable mutations generated and placed in current dirrectory! Name: "+outfilename
	else:
		print "No favorable mutations detected with the current provided threshholds. Please check if these values are accurate, and that the files are of csv type, with no extraneous spaces or punctuation. Thresholds used were: PSSM - "+pssm_threshold+" Distance to Active Site - "+dist_threshold+" Contact Number - "+contact_threshold
		return

	return outfilename

def main():
	if options.path == "GenerateMutations":
		print "Running Generate Mutations..."
		enhance_enzyme_solubility()
	if options.path == 'CalculateDistAndCN':
		print "Running Calculate Distance and Contact Number..."
		dist_contact_filters()

if __name__ == '__main__':
	main()