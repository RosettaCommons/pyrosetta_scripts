#!/home/koepnick/local/anaconda2/bin/python

#===============================================================================
# Brian Koepnick
# February 2, 2016
#
# Print ABEGO string for pdb structures
#
#===============================================================================

import rosetta
import rosetta.core.import_pose
import numpy
import os

import argparse

parser = argparse.ArgumentParser()
parser.add_argument( 'input_pdbs', type=str, nargs='+', help='input pdb files' )
args = parser.parse_args()

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
	outlines.append( pdb_id )

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
	outlines.append( abego_str )
	
print '\n'.join( outlines )

'''
# write abego data
summary_lines = []
for pdb_id in result.keys():
	summary_lines.append( '{0}\n{1}'.format( pdb_id, result[pdb_id]['abego'] ) )
outfile = args.out + '_abego.out'		
with open ( outfile, 'w' ) as f:
	f.write( '\n'.join( summary_lines ) )
'''
