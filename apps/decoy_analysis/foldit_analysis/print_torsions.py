#!/home/koepnick/local/anaconda2/bin/python

#===============================================================================
# Brian Koepnick
# September 27, 2016
#
# 	Print phi, psi, omega for each residue of input PDB.
#
#===============================================================================

import os, sys, argparse
import math, numpy
import rosetta
import rosetta.core.import_pose


if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description='')
	parser.add_argument( 'pdbs', type=str, nargs='*', help='input pdb files' )
	parser.add_argument( '-out', type=str, default='torsions.out', help='output filename' )
	#parser.add_argument( '-radians', action='store_true', help='output dihedral angles in radians' )
	args = parser.parse_args()

	rosetta.init()
	outlines = []
	for pdb_filename in args.pdbs:
		print pdb_filename
		assert( os.path.isfile( pdb_filename ) )
		pose = rosetta.core.import_pose.pose_from_pdb( pdb_filename )

		for ii in range( 1, pose.total_residue()+1 ):		
			phi = pose.phi( ii )
			psi = pose.psi( ii )
			omega = pose.omega( ii )

			outline = '{0:.6f}\t{1:.6f}\t{2:.6f}'.format( phi, psi, omega )
			outlines.append( outline )

	with open( args.out, 'w' ) as f:
		f.write( '\n'.join( outlines ) )
