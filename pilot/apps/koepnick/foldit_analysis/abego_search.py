#!/home/koepnick/local/anaconda2/bin/python

#===============================================================================
# Brian Koepnick
# February 2, 2016
#
# Search pdb structures for particular ABEGO string
#
#===============================================================================

import rosetta
import rosetta.core.import_pose
import re, numpy
import os, sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument( 'input_pdbs', type=str, nargs='*', help='input pdb files' )
parser.add_argument( '-l', type=str, help='list of pdb files' )
parser.add_argument( '-abego', type=str, required=True, help='ABEGO search string' )
args = parser.parse_args()

if args.l:
	with open( args.l, 'r' ) as f:
		args.input_pdbs.extend( [ line.strip() for line in f.readlines() ] )

if not args.input_pdbs:
	print 'Input pdb files required'
	quit()


RE_ABEGO = re.compile( args.abego.upper() )

def calc_abego( omega, phi, psi ):
	if numpy.fabs( omega ) < 90.0:
		return 'O'
	elif phi >= 0.0:
		if ( -100.0 <= psi and psi < 100.0 ):
			return 'G'
		else:
			return 'E'
	else:
		if ( -125.0 <= psi and psi < 50.0 ):
			return 'A'
		else: 
			return 'B'

rosetta.init()

outlines = []

# Get a test structure
for pdb in args.input_pdbs:
	pdb_id = os.path.splitext( os.path.basename( pdb ) )[0]

	# read in pdb
	assert( os.path.isfile( pdb ) )
	pose = rosetta.core.import_pose.pose_from_pdb( pdb )

	total_res = pose.total_residue()

	abego = []
	for ii in range( 1, total_res ):
		omega = pose.omega( ii )
		phi = pose.phi( ii )
		psi = pose.psi( ii )
		abego.append( calc_abego( omega, phi, psi ) )
	abego_str = ''.join( abego )
	#match = RE_ABEGO.search( abego_str )
	it = RE_ABEGO.finditer( abego_str )
	for match in it:
		outline = '{0}[{1}]: {2}'.format( pdb_id, match.start(), match.group() )
		print outline
	
#print '\n'.join( outlines )
