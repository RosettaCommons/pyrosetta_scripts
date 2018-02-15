#!/home/koepnick/local/anaconda2/bin/python

#===============================================================================
# Brian Koepnick
# September 3, 2016
#
# 	Used to generate x0 and SD values for backbone torsions of canonical ABEGO
# loops.
#
#===============================================================================

import os, sys, argparse
import re
import math, numpy
import rosetta
import rosetta.core.import_pose


def mean_angle( angles ):
	# converting each angle to unit vector, v is the sum of all vectors
	v = ( sum( [ math.sin( x ) for x in angles ] ), sum( [ math.cos( x ) for x in angles ] ) )

	# r is the magnitude of vector v
	r = math.sqrt( pow( v[0], 2 ) + pow( v[1], 2 ) )
	mean_r = r / len( angles )

	mean = math.atan2( v[0], v[1] )
	stdev = math.sqrt( -2 * math.log( mean_r ) )
	
	return mean, stdev


if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description='')
	parser.add_argument( 'dbfile', type=str, nargs='*', help='i.e. aa_abego_allfrequencies_nonorm' )
	parser.add_argument( '-raw_file', type=str, help='raw data from a previous run' )
	parser.add_argument( '-radians', action='store_true', help='output dihedral angles in radians' )
	parser.add_argument( '-no_raw', action='store_true', help='do not print raw data file' )
	parser.add_argument( '-verbose', action='store_true', help='verbose output' )
	args = parser.parse_args()

	if not args.verbose:
		FNULL = open( os.devnull, 'w' )
		sys.stdout = FNULL
	
	results = {}
	if args.dbfile:
		with open( args.dbfile[ 0 ], 'r' ) as dbfile:
			in_records = [ x.strip().split() for x in dbfile.readlines() ] 

		rosetta.init()

		raw_data = []
		for record in in_records:
			assert( len( record ) == 7 )
			pdb_filename, context, start, end, loopseq, extraseq, abego = tuple( record )
			
			# check for this context/abego combination in results dictionary
			if not context in results:
				results[ context ] = {}

			if not abego in results[ context ]:
				results[ context ][ abego ] = []

			pdb_id = os.path.splitext( os.path.basename( pdb_filename ) )[0]

			# read in pdb
			assert( os.path.isfile( pdb_filename ) )
			pose = rosetta.core.import_pose.pose_from_pdb( pdb_filename )

			# get one residue of context on either side of loop
			# why are we converting everything to radians here if the default is just to convert back to degrees???
			loop_bb = []
			for ii in range( int( start ) - 1, int( end ) + 2 ):
				omega = pose.omega( ii ) * math.pi / 180.0
				phi = pose.phi( ii ) * math.pi / 180.0
				psi = pose.psi( ii ) * math.pi / 180.0
				loop_bb.append( ( phi, psi, omega ) )

			# add this loop backbone to results dictionary at appropriate key
			results[ context ][ abego ].append( loop_bb )
			if args.radians:
				torsions = '\t'.join( [ '{{{0}}}'.format( ' '.join( [ '{0}'.format( torsion ) for torsion in bb ] ) ) for bb in loop_bb ] ) #radians
			else:
				torsions = '\t'.join( [ '{{{0}}}'.format( ' '.join( [ '{0}'.format( torsion * 180.0 / math.pi ) for torsion in bb ] ) ) for bb in loop_bb ] ) #degrees
			raw_data.append( ( context, abego, pdb_id, torsions ) )
	
	elif args.raw_file:
		print 'Reading in raw data file: ' + args.raw_file + '...'

		with open( args.raw_file, 'r' ) as raw_file:
			in_records = [ x.strip().split( '\t' ) for x in raw_file.readlines() ]

		re_brace = re.compile( '{|}' )
		for record in in_records:
			context = record.pop( 0 )
			abego = record.pop( 0 )
			pdb_id = record.pop( 0 )

			torsions = [ re.sub( re_brace, '', x ).split() for x in record ]

			if not context in results:
				results[ context ] = {}
			
			if not abego in results[ context ]:
				results[ context ][ abego ] = []

			loop_bb = []
			for res in torsions:
				loop_bb.append( [ float( x ) for x in res ] )

			results[ context ][ abego ].append( loop_bb )


	print 'Writing individual loop files...'
	out_records = []
	for context, abegos in results.iteritems():
		for abego, loops in abegos.iteritems():
			print 'Writing ' + context + ' ' + abego
			
			# get length of loop from first entry
			loop_mean_std = []
			for ii in range( len( loops[0] ) ):
				phi_ii = [ loop[ ii ][0] for loop in loops ]
				#phi_mean = numpy.mean( phi_ii )
				#phi_std  = numpy.std( phi_ii, ddof=1 )
				phi_mean, phi_std = mean_angle( phi_ii )

				psi_ii = [ loop[ ii ][1] for loop in loops ]
				#psi_mean = numpy.mean( psi_ii )
				#psi_std  = numpy.std( psi_ii, ddof=1 )
				psi_mean, psi_std = mean_angle( psi_ii )
				
				omega_ii = [ loop[ ii ][2] for loop in loops ]
				#omega_mean = numpy.mean( omega_ii )
				#omega_std  = numpy.std( omega_ii, ddof=1 )
				omega_mean, omega_std = mean_angle( omega_ii )

				loop_mean_std.append( ( phi_mean, phi_std, psi_mean, psi_std, omega_mean, omega_std ) )

			prefix, suffix = ( 'A', 'A' )
			if context[0] == 'E':
				prefix = 'B'
			if context[1] == 'E':
				suffix = 'B'
			outfile = prefix + abego + suffix + '.dat'

			with open( outfile, 'w' ) as f:
				if args.radians:
					f.write( '\n'.join( [ '\t'.join( [ '{0:7.6f}'.format( value ) for value in ii ] ) for ii in loop_mean_std ] ) )	
				else:
					f.write( '\n'.join( [ '\t'.join( [ '{0:8.3f}'.format( value * 180.0 / math.pi ) for value in ii ] ) for ii in loop_mean_std ] ) )	
	if not args.no_raw:		
		#outlines = [ 'context\tabego\tpdb\ttorsions' ]
		with open( 'raw.dat', 'w' ) as f:
			f.write( '\n'.join( [ '\t'.join( datum ) for datum in raw_data ] ) )
