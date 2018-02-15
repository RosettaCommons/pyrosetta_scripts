#!/home/koepnick/local/anaconda2/bin/python

import os, sys
import argparse    
import re
import numpy as np
from rosetta import *
from pyrosetta import *
from sklearn import manifold

RE_ATOM_RECORD = re.compile( 'ATOM' )
RE_IRDATA_SCORE = re.compile( 'IRDATA SCORE ' )
RE_IRDATA_ENERGY = re.compile( 'IRDATA ENERGY ' )
RE_IRDATA_PDL = re.compile( 'IRDATA PDL ' )

# convert tm_score to distance measure
def dist( tm_score ):
	dist = 1.0 - pow( tm_score, 2 )
	#dist = 1.0 - tm_score

	return dist


# get PDL and energy from pdb files
def get_irdata( pdb_list ):
	print 'Retrieving IRDATA...'

	irdata_dict = {}

	for pdb in pdb_list:
		pdb_irdata = {}
		pdb_irdata[ 'irdata_pdl' ] = []

		with open( pdb, 'r' ) as f:
			for line in f.readlines():
				if RE_ATOM_RECORD.match( line ):
					continue
				elif RE_IRDATA_PDL.match( line ):
					pdb_irdata[ 'irdata_pdl' ].append( line.lstrip( 'IRDATA_PDL . ' ).split( ',' )[0] )
				elif RE_IRDATA_SCORE.match( line ):
					pdb_irdata[ 'irdata_score' ] = float( line.strip().split()[2] )
				elif RE_IRDATA_ENERGY.match( line ):
					pdb_irdata[ 'irdata_energy' ] = float( line.strip().split()[2] )

		irdata_dict[ pdb ] = pdb_irdata

	return irdata_dict


# calculate TMscore
# NB: this Rosetta implementation of TMalign yields different numbers than the original FORTRAN program
def tm_score( decoy_pose, ref_pose ):
	tm_align = rosetta.protocols.hybridization.TMalign()

	# check return value of tm_align() for failure
	if not tm_align.apply( decoy_pose, ref_pose ) == 0:
		raise ValueError( 'TMalign failed...' )

	return tm_align.TMscore( ref_pose.total_residue() )


# lookup distance in dictionary
def get_dist( d, key1, key2 ):
	try:
		if key1 in d and key2 in d[ key1 ]:
			return d[ key1 ][ key2 ]
		else:
			return d[ key2 ][ key1 ]
	except:
		raise ValueError( 'Unable to find ({0}, {1}) in dictionary!'.format( key1, key2 ) )


# convert distance dictionary to ndarray matrix
# returns a distance matrix and an index
def dict2mtx( dictionary ):
	# make an index of dictionary keys
	idx = dictionary.keys()

	# now make a 2D ndarray of size NxN, populate with zeros
	shape = ( len( idx ), ) * 2
	mtx = np.zeros( shape )

	# fill in southwest corner of distance matrix
	for i, key1 in enumerate( idx ):
		for j, key2 in list( enumerate( idx ) )[:i]:
			mtx[ i ][ j ] = get_dist( dictionary, key1, key2 )
		
	# symmetrize distance matrix	
	dist_mtx = np.maximum( mtx, mtx.transpose() )
	
	return dist_mtx, idx


# write distance matrix to file
def write_distance_matrix( idx, dist_mtx, filename ):
	
	# write index as first line
	outlines = [ '\t'.join( idx ) ]

	for i in range( len( idx ) ):
		mtx_row = []
		for j in range( i ):
			mtx_row.append( '{0:0.4f}'.format( dist_mtx[ i ][ j ] ) )

		outlines.append( '\t'.join( mtx_row ) )		
	
	with open( filename, 'w' ) as f:
		f.write( '\n'.join( outlines ) )
		

# try to place pose in a cluster
def cluster_pose( decoy_id, decoy_pose, dicts, threshold ):
	decoy_dict, dist_dict, energy_dict = dicts

	dists = {}
	for ref_id, ref_pose in decoy_dict.iteritems():
		tm = tm_score( decoy_pose, ref_pose )

		# if decoy belongs in existing cluster, return
		if tm >= threshold:
			print 'Rejected as duplicate of {0}\nTMscore: {1}'.format( ref_id, tm )
			return
	
		# store distance
		dists[ ref_id ] = dist( tm )

	# if no parent cluster found, add decoy to list of cluster centers
	decoy_dict[ decoy_id ] = decoy_pose
	print 'Accepting cluster #{0}: {1}'.format( len( decoy_dict ), decoy_id )

	# store distances and energy
	dist_dict[ decoy_id ] = dists


# calculate pairwise distances between a set of pdbs
def calc_distances( decoy_pdbs, dicts ):	
	decoy_dict, dist_dict, energy_dict = dicts
	print 'Begin distance calculation...'

	for decoy_id in decoy_pdbs:
			
		# if for some reason we've already evaluated this decoy, skip
		if decoy_id in decoy_dict:
			continue

		print 'Evaluating {0}...'.format( decoy_id )

		dists = {}

		# calculate distances to all known decoys
		decoy_pose = rosetta.core.import_pose.pose_from_file( decoy_id )
		for ref_id, ref_pose in decoy_dict.iteritems():
			dists[ ref_id ] = dist( tm_score( decoy_pose, ref_pose ) )

		# store data
		decoy_dict[ decoy_id ] = decoy_pose
		dist_dict[ decoy_id ] = dists


# calculate pairwise distances between a set of pdbs, with clustering
def calc_cluster_distances( sorted_decoy_pdbs, dicts, threshold, energy_depth, max_clust ):	
	decoy_dict, dist_dict, energy_dict = dicts
	print 'Begin clustering and distance calculation...'
	
	# set energy baseline
	energy_baseline = energy_dict[ sorted_decoy_pdbs[0] ]

	# compare each new structure to known cluster centers
	for decoy_id in sorted_decoy_pdbs:
		
		# if we've surpassed the energy depth, break
		if abs( energy_dict[ decoy_id ] - energy_baseline ) > energy_depth:
			break

		# if for some reason we've already accepted this decoy, skip
		if decoy_id in decoy_dict:
			continue

		# initialize pose
		decoy_pose = rosetta.core.import_pose.pose_from_file( decoy_id )
		
		# try to cluster pose	
		print 'Evaluating {0}...'.format( decoy_id )		
		cluster_pose( decoy_id, decoy_pose, dicts, threshold )		

		# break as soon as we've accumulated the maximum number of clusters
		if len( decoy_dict ) >= max_clust:
			break


# for rosetta input
def rosetta_mds( args, dicts ):
	decoy_dict, dist_dict, energy_dict = dicts

	# first initialize a key for the native
	if args.native:
		print 'Evaluating native: {0}'.format( args.native )
		native_pose = rosetta.core.import_pose.pose_from_file( args.native )
		
		dist_dict[ 'NATIVE' ] = {}
		decoy_dict[ 'NATIVE' ] = native_pose

		scorefxn = rosetta.core.scoring.get_score_function()
		energy_dict[ 'NATIVE' ] = scorefxn( native_pose )

	# do any clustering on provided silent file
	if args.silent:
		# read in silent file
		sfd = rosetta.core.io.silent.SilentFileData()
		sfd.read_file( args.silent )

		# get a list of tags and sort by energy
		sf_tags = list( sfd.tags() )
		sf_tags.sort( key=lambda x: sfd.get_structure( x ).get_energy( 'score' ) )
		baseline = sfd.get_structure( sf_tags[0] ).get_energy( 'score' )

		# iterate over silent structures
		for tag in sf_tags:
			silent_struct = sfd.get_structure( tag )
			sc = silent_struct.get_energy( 'score' )

			# break as soon as we've breached the specified energy depth
			if args.energy_depth < abs( baseline - sc ):
				break
			energy_dict[ tag ] = sc

			# initialize pose
			decoy_pose = rosetta.core.pose.Pose()
			silent_struct.fill_pose( decoy_pose )

			# try to cluster pose
			print 'Evaluating {0}...'.format( tag )
			cluster_pose( tag, decoy_pose, dicts, args.clust_thresh )

			# break as soon as we've accumulated the maximum number of clusters
			if len( decoy_dict ) >= args.max_clust:
				break

	# finally handle provided decoys
	if args.l:
		with open( args.l ) as f:
			decoy_pdbs = [ line.strip() for line in f.readlines() ]

			scorefxn = rosetta.core.scoring.get_score_function()
			for decoy_id in decoy_pdbs:
				
				# if for some reason we've already evaluated this decoy, skip
				if decoy_id in decoy_dict:
					continue
				
				decoy_pose = rosetta.core.import_pose.pose_from_file( decoy_id )
				energy_dict[ decoy_id ] = scorefxn( decoy_pose )

				print 'Evaluating {0}...'.format( decoy_id )
				dists = {}

				# calculate distances to all known decoys
				for ref_id, ref_pose in decoy_dict.iteritems():
					dists[ ref_id ] = dist( tm_score( decoy_pose, ref_pose ) )

				# store data
				decoy_dict[ decoy_id ] = decoy_pose
				dist_dict[ decoy_id ] = dists
				

# for foldit input
def foldit_mds( args, dicts ):
	decoy_dict, dist_dict, energy_dict = dicts

	# first initialize a key for the native
	if args.native:
		print 'Evaluating native: {0}'.format( args.native )
	
		dist_dict[ 'NATIVE' ] = {}
		decoy_dict[ 'NATIVE' ] = rosetta.core.import_pose.pose_from_file( args.native )

		irdata_dict = get_irdata( [ args.native ] )
		energy_dict[ 'NATIVE' ] = irdata_dict[ args.native ][ args.energy_type ]
	
	# then do any clustering
	if args.cluster:
		with open( args.cluster ) as f:
			decoy_pdbs = [ line.strip() for line in f.readlines() ]
			
			# collect energy data for decoys
			irdata_dict = get_irdata( decoy_pdbs )
			for pdb_id, irdata in irdata_dict.iteritems():
				energy_dict[ pdb_id ] = irdata[ args.energy_type ]
			
			# sort decoys
			sorted_decoy_pdbs = sorted( decoy_pdbs, key=energy_dict.get, reverse=( args.energy_type=='irdata_score' ) )
			
			calc_cluster_distances( sorted_decoy_pdbs, dicts, args.clust_thresh, args.energy_depth, args.max_clust )

	# finally handle provided decoys
	if args.solutions:
		with open( args.solutions ) as f:
			decoy_pdbs = [ line.strip() for line in f.readlines() ]
			
			# collect score data for decoys
			irdata_dict = get_irdata( decoy_pdbs )
			for pdb_id, irdata in irdata_dict.iteritems():
				energy_dict[ pdb_id ] = irdata[ args.energy_type ]

			calc_distances( decoy_pdbs, dicts )

	
# master function
def run( args ):
	# initialize rosetta
	#init( '-chemical:include_patches patches/NtermProteinFull.txt patches/CtermProteinFull.txt' )	
	init( '-in:file:silent_read_through_errors' )

	# create dictionaries for decoy poses, energies, and distances
	decoy_dict = {}
	dist_dict = {}
	energy_dict = {}
	dicts = ( decoy_dict, dist_dict, energy_dict )	

	# support silent file input?
	if args.silent or args.l:
		rosetta_mds( args, dicts )
	
	if args.solutions or args.cluster: 
		foldit_mds( args, dicts )

	# convert dictionary to distance matrix and index
	dist_mtx, idx = dict2mtx( dist_dict )

	# dump distance matrix?
	if args.dump_mtx:
		write_distance_matrix( idx, dist_mtx, args.dump_mtx )

	# now run MDS
	mds = manifold.MDS( n_components=2, n_jobs=1, dissimilarity='precomputed', random_state=11 )
	coords = mds.fit_transform( dist_mtx )
		
	print 'Distance Matrix:'
	print dist_mtx

	print 'MDS coordinates:'
	print coords

	# add id and score component to ndarray
	# if available, write PDL username for pdb_id
	pdb_id = np.array( [[ os.path.basename( key ) for key in idx ]] )
	sc = np.array( [[ energy_dict[ key ] for key in idx ]] )

	coords = np.concatenate( ( pdb_id.T, coords, sc.T ), axis=1 )
	print coords
	
	# write tab-delimited coordinates to file
	with open( args.out, 'w' ) as f:
		outlines = [ '\t'.join( coord ) for coord in coords ]
		f.write( '\n'.join( outlines ) )


if __name__ =='__main__':
	parser = argparse.ArgumentParser(description='')

	# arguments
	parser.add_argument( '-l', type=str, help='input list of pdbs' )
	parser.add_argument( '-solutions', type=str, help='input list of Foldit solution pdbs' )
	parser.add_argument( '-cluster', type=str, help='input list of Foldit solution pdbs to be clustered' )
	parser.add_argument( '-silent', type=str, help='input silent file to be clustered' )
	parser.add_argument( '-native', type=str, help='native pdb file' )
	parser.add_argument( '-energy_type', choices=[ 'irdata_energy', 'irdata_score', 'score' ], default='irdata_score', help='value to use for decoy energy; also used in clustering' )
	parser.add_argument( '-out', type=str, default='mds.dat', help='output file for MDS coordinates' )
	parser.add_argument( '-dump_mtx', type=str, help='dump the raw distance matrix to file' )
	parser.add_argument( '-energy_depth', type=float, default=500, help='maximum energy depth for clustering' )
	parser.add_argument( '-clust_thresh', type=float, default=0.90, help='TM-score threshold, above which two pdbs are considered members of the same cluster' ) 
	parser.add_argument( '-max_clust', type=int, default=200, help='the maximum number of cluster centers to plot' )
	parser.add_argument( '-verbose', action='store_true', help='verbose output' )

	args = parser.parse_args()

	# check for valid arguments
	if not ( args.l or args.solutions or args.cluster or args.silent ):
		parser.print_help()
		raise ValueError( 'Invalid input arguments' )

	# unless using verbose option, redirect stdout to /dev/null
	if not args.verbose:
		fnull = open( os.devnull, 'w' )
		sys.stdout = fnull

	run( args )	
