# @file: calculate_flex_ddG.py
# @author: Rebecca Alford (ralford3#@jhu.edu)
# @brief: Calculate the ddG of mutation from an ensemble of backrub models
# Uses both ensemble averages and the generalized additive model for calculation

import sys, os
import numpy as np
import sqlite3

from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

#################################################################################
## Glaobl Data - Generalized Additive Model Constants per Score Type
## From Barlow et al. BioRxiv 2017
a = { "fa_sol" : 6.940, "hbond_sc" : 1.902, "hbond_bb_sc" : 0.063, "fa_rep" : 1.659, "fa_elec" : 0.697, "hbond_lr_bb" : 2.738, "fa_atr" : 2.313 }
b = { "fa_sol" : 6.722, "hbond_sc" : 1.999, "hbond_bb_sc" : 0.452, "fa_rep" : 1.836, "fa_elec" : 0.122, "hbond_lr_bb" : 1.179, "fa_atr" : 1.649 }
#################################################################################

def main( args ): 

	parser = OptionParser( usage="usage: %prog --dbfile ddG.db3" )
	parser.set_description(main.__doc__)

	parser.add_option('--dbfile', '-d', 
		action="store",
		help="Path to sqlite3 datbase file"
		)

	(options, args) = parser.parse_args(args=args[1:])
	global Options
	Options = options 

	if ( not Options.dbfile ): 
		sys.exit( "Must provide flags --dbfile! Exiting..." )

	# Connect to the sqlite database
	conn = sqlite3.connect( Options.dbfile )
	c = conn.cursor()

	# Calculate ensemble statistics for the native and mutant poses with normal weighting
	native_model_scores, mutant_model_scores = calculate_weighted_ensembles( c )
	print mutant_model_scores.mean() - native_model_scores.mean()
	print mutant_model_scores.std()
	print native_model_scores.std()

	# Calculate ensemble statistics for the native and mutant poses with GAM weighting
	native_GAM_model_scores, mutant_GAM_model_scores = calculate_GAM_weighted_ensembles( c )
	print mutant_GAM_model_scores.mean() - native_GAM_model_scores.mean()
	print mutant_GAM_model_scores.std()
	print native_GAM_model_scores.std()

#################################################################################
# @brief: Get energy function weights from the datbaase table
def get_score_function_weights( c ): 

	# Get the score type IDs from the score_function_weights table
	c.execute('SELECT score_type_id FROM score_function_weights')
	score_type_ids = [ x[0] for x in c.fetchall() ]

	# Get the score type weights from the score_function_weights table
	c.execute('SELECT weight FROM score_function_weights')
	score_type_weights = [ x[0] for x in c.fetchall() ]

	# Create and return a dictionary mapping type ids to type weights
	return dict(zip(score_type_ids, score_type_weights))

#################################################################################
# @brief: Get a mapping between the score type name and the score type id (according
# to the sqlite3 database)
def map_score_type_name_to_id( c ): 

	# Make an empty dictionary for the mapping
	score_type_name_to_id = {}

	# Get score function weights
	sfxn_weights = get_score_function_weights( c )

	for st_id in sfxn_weights: 

		c.execute('SELECT (score_type_name) FROM score_types WHERE score_type_id={st_id}'.format(st_id=st_id))
		score_type_name = str(c.fetchall()[0][0])
		score_type_name_to_id[ score_type_name ] = int(st_id)

	return score_type_name_to_id

#################################################################################
# @brief: Calculate the ensemble average ddG given an sqlite database of weighted scores
# Also calculate ensemble statistics based on individual score terms
def calculate_weighted_ensembles( c ):

	# Figure out how many strucutres are in the ensemble
	c.execute('SELECT ({coi}) FROM {tn}'.format(coi='struct_id',tn='structure_scores'))
	number_of_structures = len(list(set(c.fetchall())))
	if ( number_of_structures % 3 != 0 ): 
		sys.exit( "Incomplete backrub ensemble database: Number of entries should be a multiple of three!" )
	ensemble_size = number_of_structures/3

	# Make a list of indices that coontain the wild type and mutant models
	mutant_indices = []
	wt_indices = []
	for i in range(1, ensemble_size+1):
		mutant_indices.append( i*3 ) 
		wt_indices.append( (i*3) - 1 )

	# For each structure in the ensemble, calculate the total score
	native_total_scores = []
	mutant_total_scores = []
	for i in range(1,number_of_structures+1): 

		# Get all of the non-zero scores for struct_id = i
		c.execute('SELECT ({coi1}) FROM {tn} WHERE {coi2}={struct_id}'.\
			format(coi1='score_value', tn='structure_scores', coi2='struct_id', struct_id=str(i)))
		raw_scores = [ x[0] for x in c.fetchall() ]

		# Store in the appropriate score array
		if ( i in wt_indices ):
			native_total_scores.append( sum(raw_scores) )
		elif ( i in mutant_indices ): 
			mutant_total_scores.append( sum(raw_scores) )

	# Currently calculating stats here, but will move later
	native_arr = np.array( native_total_scores )
	mutant_arr = np.array( mutant_total_scores )

	return native_arr, mutant_arr
	
#################################################################################
# @brief: A logistic function that transforms a score type into a symmetrically 
# weighted value given fitted constants a and b 
# Formula: y = -e^a + (2*e^a)/(1 + e^(-score * e^-b ) )
def y_of_x( a, b, score ): 

	const = np.exp( a )
	y = -const + (2*const)/( 1 + np.exp( -score * np.exp( -b ) ) )
	return y

#################################################################################
# @brief: Calculate the GAM ensemble average ddG given an sqlite database of weighted scores
# Also calculate ensemble statistics based on individual score terms
def calculate_GAM_weighted_ensembles( c ):

	# Figure out how many strucutres are in the ensemble
	c.execute('SELECT ({coi}) FROM {tn}'.format(coi='struct_id',tn='structure_scores'))
	number_of_structures = len(list(set(c.fetchall())))
	if ( number_of_structures % 3 != 0 ): 
		sys.exit( "Incomplete backrub ensemble database: Number of entries should be a multiple of three!" )
	ensemble_size = number_of_structures/3

	# Make a list of indices that coontain the wild type and mutant models
	mutant_indices = []
	wt_indices = []
	for i in range(1, ensemble_size+1):
		mutant_indices.append( i*3 ) 
		wt_indices.append( (i*3) - 1 )

	# Get a mapping between score type IDs and names
	score_type_name_to_id =	map_score_type_name_to_id( c )

	# For each structure in the ensemble, calculate the total score
	native_total_scores = []
	mutant_total_scores = []
	for i in range(1,number_of_structures+1): 

		# Get all of the score_type_ids for all non-zero scores when struct_id = i
		c.execute('SELECT ({coi1}) FROM {tn} WHERE {coi2}={struct_id}'.\
			format(coi1="score_type_id", tn="structure_scores", coi2="struct_id", struct_id=str(i)))
		score_type_ids = [ x[0] for x in c.fetchall() ]

		# Get all of the non-zero scores for struct_id = i
		c.execute('SELECT ({coi1}) FROM {tn} WHERE {coi2}={struct_id}'.\
			format(coi1='score_value', tn='structure_scores', coi2='struct_id', struct_id=str(i)))
		raw_scores = [ x[0] for x in c.fetchall() ]

		# Make a compact dictionary for easy search
		unweighted_scores = dict(zip(score_type_ids, raw_scores))

		# Calculate the GAM weighted total score from the individual scores
		gam_score_types = [ 'fa_sol', 'hbond_sc', 'hbond_bb_sc', 'fa_rep', 'fa_elec', 'hbond_lr_bb', 'fa_atr' ]

		total_model_score = 0
		for score_type in gam_score_types: 

			score_type_id = score_type_name_to_id[ score_type ]
			unwt_score = unweighted_scores[ score_type_id ]
			const_a = a[ score_type ]
			const_b = b[ score_type ]
			gam_weighted_score = y_of_x( const_a, const_b, unwt_score )
			total_model_score += gam_weighted_score

		# Store in the appropriate score array
		if ( i in wt_indices ):
			native_total_scores.append( total_model_score )
		elif ( i in mutant_indices ): 
			mutant_total_scores.append( total_model_score )

	# Currently calculating stats here, but will move later
	native_arr = np.array( native_total_scores )
	mutant_arr = np.array( mutant_total_scores )

	return native_arr, mutant_arr

if __name__ == "__main__" : main(sys.argv)
