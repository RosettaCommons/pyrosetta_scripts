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
import matplotlib.pyplot as plt

import argparse

parser = argparse.ArgumentParser()
parser.add_argument( 'input_pdbs', type=str, nargs='+', help='input pdb files' )
parser.add_argument( '-plot_off', action='store_true', default=False, help='suppress PNG output files' )
args = parser.parse_args()

#Load vall structure database
sdb = StructureDatabase("/work/fordas/databases/structure/vall.h5")

vall_residues = sdb.residues.read()

# Specify fragment type for extraction
test_fragment_type = FragmentSpecification(4, "CA")

all_vall_fragment_start_residues, all_vall_fragments = test_fragment_type.fragments_from_source_residues( vall_residues )

# Get a test structure
rosetta.init()

summary_out = [ 'pdb_id\tworst_fragment_rmsd' ]
for pdb in args.input_pdbs:
	assert( os.path.isfile( pdb ) )
	pdb_id = os.path.splitext( os.path.basename( pdb ) )[0]
	outfile = pdb_id + '_4mer.png'
	
	#test_structure_pose = rosetta.core.import_pose.pose_from_pdb("/work/koepnick/projects/ideal_fragments/2lta_rsmn3x1_NMR.pdb")
	test_structure_pose = rosetta.core.import_pose.pose_from_pdb( pdb )
	total_res = test_structure_pose.total_residue()
	test_structure_residues = StructureDatabase.extract_residue_entries_from_pose( test_structure_pose )

	test_fragment_start_indicies, test_fragments = test_fragment_type.fragments_from_source_residues( test_structure_residues )

	test_to_vall_alignment = atom_array_broadcast_rmsd( all_vall_fragments["coordinates"], test_fragments[:total_res]["coordinates"] )

	min_alignment = test_to_vall_alignment.min(axis=0)

	# Vall contains ~3800000 9mers
	# 0.001 percentile = top 38 9mers
	high_percentile_alignment = numpy.percentile(test_to_vall_alignment, .001, axis=0)

	summary_out.append( '{0}\t{1}'.format( pdb_id, max( high_percentile_alignment ) ) )

	if not args.plot_off:
		plt.clf()
		plt.plot( high_percentile_alignment, color='blue', label='.001 percentile' )
		plt.plot( min_alignment, color='red', label='best fragment' )
		plt.ylim(0,0.2)
		plt.legend()
		plt.savefig( outfile )

with open ( 'summary.out', 'w' ) as f:
	f.write( '\n'.join( summary_out ) )
