#!/work/fordas/virtualenvs/fragment_profiling_v0.1/bin/python

#===============================================================================
# Before running, must source Alex Ford's virtual environment:
#	source /work/fordas/virtualenvs/fragment_profiling_v0.1/bin/activate
#
# Usage:
# 	profile_vall.py file1.pdb [ file2.pdb ... ]
#
#===============================================================================

import interface_fragment_matching.fragment_fitting
import interface_fragment_matching.structure_database
import interface_fragment_matching.fragment_fitting.rmsd_calc
import tables
import rosetta
import rosetta.core.import_pose
from interface_fragment_matching.structure_database import StructureDatabase
from interface_fragment_matching.fragment_fitting import FragmentSpecification
from interface_fragment_matching.fragment_fitting.rmsd_calc import atom_array_broadcast_rmsd
import numpy
import os
#import matplotlib.pyplot as plt

import argparse

parser = argparse.ArgumentParser()
parser.add_argument( 'input_pdbs', type=str, nargs='+', help='input pdb files' )
parser.add_argument( '-rmsd_thresh', type=float, default=0.5, help='maximum tolerated RMSD for fragments' )
parser.add_argument( '-frag_size', type=int, help='fragment length' )
parser.add_argument( '-out', type=str, default='fragments', help='id to prepend to output filename' )
args = parser.parse_args()

if not args.frag_size:
	args.frag_size = 9

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

#Load vall structure database
sdb = StructureDatabase("/work/fordas/databases/structure/vall.h5")

vall_residues = sdb.residues.read()
rosetta.init()

summary_out = []
result = {}

test_fragment_type = FragmentSpecification( args.frag_size, "CA" )
all_vall_fragment_start_residues, all_vall_fragments = test_fragment_type.fragments_from_source_residues( vall_residues )

#test_fragment_type_4mer = FragmentSpecification(4, "CA")
#all_vall_fragment_start_residues_4mer, all_vall_fragments_4mer = test_fragment_type_4mer.fragments_from_source_residues( vall_residues )

# Get a test structure
for pdb in args.input_pdbs:
	pdb_id = os.path.splitext( os.path.basename( pdb ) )[0]

	# Vall contains ~3800000 9mers
	# 0.001 percentile = top 38 9mers
	# high_percentile_alignment = numpy.percentile(test_to_vall_alignment, .001, axis=0)

	if pdb_id in result:
		test_structure_pose = result[pdb_id]['pose']
	else:
		# read in pdb
		assert( os.path.isfile( pdb ) )
		test_structure_pose = rosetta.core.import_pose.pose_from_pdb( pdb )

		# initialize dict entry
		pose = test_structure_pose
		abego = []
		result[pdb_id] = { 'pose' : pose, 'abego' : abego }

	total_res = test_structure_pose.total_residue()
	test_structure_residues = StructureDatabase.extract_residue_entries_from_pose( test_structure_pose )

	test_fragment_start_indices, test_fragments = test_fragment_type.fragments_from_source_residues( test_structure_residues )

	test_to_vall_alignment = atom_array_broadcast_rmsd( all_vall_fragments["coordinates"], test_fragments[:total_res]["coordinates"] )

	min_alignment = test_to_vall_alignment.min(axis=0)
	#worst_frag_rmsd = max( min_alignment )
	#mean_frag_rmsd = numpy.mean( min_alignment )

	for frame, rmsd in enumerate( min_alignment, start=1 ):
		if rmsd > args.rmsd_thresh:
			abego = []
			for ii in range( frame, frame + args.frag_size ):
				omega = test_structure_pose.omega( ii )
				phi = test_structure_pose.phi( ii )
				psi = test_structure_pose.psi( ii )
				abego.append( calc_abego( omega, phi, psi ) )
			abego_string = ''.join( abego )
			#print "Violation at residue {0}, RMSD {1:5.3f}, {2}".format( frame, rmsd, abego_string )
			result[pdb_id]['abego'].append( ( frame, rmsd, abego_string ) )

# write abego data
summary_lines = []
for pdb_id in result.keys():
	for data in result[pdb_id]['abego']:
		summary_lines.append( '{0}\t{1[0]}\t{1[1]:5.3f}\t{1[2]}'.format( pdb_id, data ) )
outfile = args.out + '_abego.out'		
with open ( outfile, 'w' ) as f:
	f.write( 'pdb_id\tres\trmsd\tabego\n' )
	f.write( '\n'.join( summary_lines ) )
