#!/home/koepnick/local/anaconda2/bin/python

#===============================================================================
# Brian Koepnick
# December 5, 2015
#
# Edited May 17, 2016:
#	Use Rosetta to calculate phi/psi instead of BioPython
#
# Scan a structure for regions that violate ideal ABEGO rules
#
# Output: The ABEGO string of the input structure, and for each 5mer that 
# 	violates ideal ABEGO rules: 
#	1. the sequence position of the 5mer
#	2. the ABEGO string of the 5mer
#
#===============================================================================

import os, sys, argparse
import re, numpy
#import Bio
#from Bio.PDB.PDBParser import PDBParser
#from Bio.PDB.Polypeptide import PPBuilder
import rosetta
import rosetta.core.import_pose

parser = argparse.ArgumentParser(description='')
parser.add_argument( 'pdbs', type=str, nargs='*', help='a pdb file' )
parser.add_argument( '-l', type=str, help='a list of pdb filenames' )
parser.add_argument( '-out', type=str, default='abego_violations.dat', help='output file' )
parser.add_argument( '-thresh', type=int, default=0, help='number of ABEGO violations to allow' )
parser.add_argument( '-verbose', action='store_true', help='output violation info' )
args = parser.parse_args()

if not args.verbose:
	FNULL = open( os.devnull, 'w' )
	sys.stdout = FNULL

pdb_list = []
if args.l:
	with open( args.l, 'r' ) as listfile:
		pdb_list = [ x.strip() for x in listfile.readlines() ]
 
pdb_list.extend( args.pdbs )
assert pdb_list

hh_loops = [ 'B','E','G','BB','GB','GBB','BAB','GABB','BBBB','GBBBB' ]
he_loops = [ 'GB','BA','GB','GBA','BAA' ]
eh_loops = [ 'AB','BBB','AB','GBB' ]
ee_loops = [ 'GG','EA','AA','BG','AAAG','AA','BAAGB' ]

extended_loops = []
for loop in hh_loops:
	extended_loops.append( 'AAAA' + loop + 'AAAA' )
for loop in he_loops:
	extended_loops.append( 'AAAA' + loop + 'BBBB' )
for loop in eh_loops:
	extended_loops.append( 'BBBB' + loop + 'AAAA' )
for loop in ee_loops:
	extended_loops.append( 'BBBB' + loop + 'BBBB' )


ideal_abegos = set([ 'AAAAA','BBBBB' ])
for loop in extended_loops:
	for i in range( len( loop ) - 4 ):
		ideal_abegos.add( loop[ i:i+5 ] )

RE_IDEAL = re.compile( '|'.join( ideal_abegos ) )

def calc_abego( omega, phi, psi ):
    if numpy.fabs( omega ) < 90.0:
        return 'O'
    elif phi >= 0.0:
        if ( -100.0 <= psi and psi < 100.0 ):
            return 'G'
        else:
            return 'E'
    else:
        if ( -75.0 <= psi and psi < 50.0 ):
            return 'A'
        else:
            return 'B'

#ppb = PPBuilder()
rosetta.init()

pdb_violations = []
pass_pdbs = []
for pdb_filename in pdb_list:
	print pdb_filename
	violations = 0

#	pdb_dir = os.path.dirname( pdb_filename )
#	pdb_id = os.path.basename( pdb_filename ).split('.')[0]
#	print pdb_id

	pdb_id = os.path.splitext( os.path.basename( pdb_filename ) )[0]

	# read in pdb
	assert( os.path.isfile( pdb_filename ) )
	pose = rosetta.core.import_pose.pose_from_pdb( pdb_filename )

#	structure = Bio.PDB.PDBParser().get_structure( 'structure', pdb_filename )
#	polypeptide = ppb.build_peptides( structure )[0]
#	phipsi_list = polypeptide.get_phi_psi_list()
	
	total_res = pose.total_residue()

	abego = []
	for ii in range( 1, total_res ):
		omega = pose.omega( ii )
		phi = pose.phi( ii )
		psi = pose.psi( ii )
		abego.append( calc_abego( omega, phi, psi ) )
	abego_str = ''.join( abego )
	
#	abegos = []
#	for phipsi in phipsi_list:
#		if not ( phipsi[0] and phipsi[1] ):
#			continue
#		abegos.append( calc_abego( 180.0, phipsi[0]*180.0, phipsi[1]*180 ) )
#	print ''.join( abegos )

	# ignore terminal residues
	#abegos = abegos[1:-1]

	for i in range( 1, len( abego ) - 4 ):
		substr = ''.join( abego[ i:i+5 ] )
		if not RE_IDEAL.match( substr ):
			violations += 1
			print i+1, substr

	pdb_violations.append( ( pdb_id, abego_str, violations ) )

#pass_pdbs = filter( lambda x: x[1] <= args.thresh, pdb_violations )
outlines = [ 'pdb_id\tabego_str\tviolations' ]
with open( args.out, 'w' ) as f:
	outlines.extend( [ '{0[0]}\t{0[1]}\t{0[2]}'.format( data ) for data in pdb_violations ] )
	f.write( '\n'.join( outlines ) )
