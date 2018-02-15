#!/home/koepnick/local/anaconda2/bin/python

#===============================================================================
# Brian Koepnick
# October 7, 2016
# Edited: November 8, 2017
#
# Tabulate the DSSP assigment and backbone torsions for each residue of input 
# structures.
#
#===============================================================================

from pyrosetta import *
import os,sys,argparse    

parser = argparse.ArgumentParser(description='')

# basic argument types
parser.add_argument('pdbs', type=str, nargs='*', help='input pdbs')
parser.add_argument('-l', type=str, help='list of input pdbs')
parser.add_argument('-multi_out', action='store_true', help='output a file for each SS type')
parser.add_argument('-prefix' , type=str, help='prefix for output file')
parser.add_argument('-verbose', action='store_true', help='verbose output')

args = parser.parse_args()

if not args.prefix:
	args.prefix = 'dssp'

if args.l:
	with open( args.l, 'r' ) as f:
		args.pdbs.extend( [ line.strip() for line in f.readlines() ] )

if not args.pdbs:
	parser.print_help()
	sys.exit()

if not args.verbose:
	FNULL = open( os.devnull, 'w' )
	sys.stdout = FNULL

pyrosetta.init()

data = { 'H': [], 'E': [], 'L': [] }
for pdb in args.pdbs:
	assert os.path.exists( pdb )

	pose = rosetta.core.import_pose.pose_from_file( pdb )
	dssp = rosetta.core.scoring.dssp.Dssp( pose )
	dssp.dssp_reduced()
	secstruct = dssp.get_dssp_secstruct()

	try:
		#assert( pose.total_residue() == len( secstruct ) ) # waters are considered as residues, but are not acknowledged by DSSP
		
		for ii in range( len( secstruct) ):
			resnum = ii+1
			ss = secstruct[ ii ]
			phi = pose.phi( resnum )
			psi = pose.psi( resnum )
			omega = pose.omega( resnum )
			while( omega < 0 ):
				omega += 360.0

			record = ( pdb, resnum, ss, phi, psi, omega )
			data[ ss ].append( record )
	except:
		print 'ERROR: {0} has {1} residues, but secstruct of length {2}'.format( pdb, pose.total_residue(), len( secstruct ) )

all_outlines = []	
header = [ '#' + '\t'.join( ( 'resid', 'ss', 'phi', 'psi', 'omega' ) ) ]

for key, value in data.iteritems():
	outfile = args.prefix + '_' + key + '.dat'
	outlines=[]
	for record in value:
		outlines.append( '{0[0]}:{0[1]}\t{0[2]}\t{1}'.format( record, '\t'.join( [ '{0:8.3f}'.format( val ) for val in record[-3:] ] ) ) )
	
	if args.multi_out:
		with open( outfile, 'w' ) as f:
			f.write( '\n'.join( header + outlines ) )
	
	all_outlines.extend( outlines )

outfile = args.prefix + '_torsions.dat'
with open( outfile, 'w' ) as f:
	f.write( '\n'.join( header + all_outlines ) )
