#!/work/fordas/virtualenvs/fragment_profiling_v0.1/bin/python

#===============================================================================
# Before running, source Alex Ford's virtual environment:
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
parser.add_argument( 'input_pdbs', type=str, nargs='*', help='input pdb files' )
parser.add_argument( '-l', type=str, help='list of pdb files' )
parser.add_argument( '-frag_size', type=int, nargs='*', help='fragment length' )
parser.add_argument( '-out', type=str, default='fragments', help='id to prepend to output filename' )
parser.add_argument( '-plot_off', action='store_true', default=False, help='suppress PNG output files' )
args = parser.parse_args()

if args.l:
	with open( args.l, 'r' ) as f:
		args.input_pdbs.extend( [ line.strip() for line in f.readlines() ] )

if not args.input_pdbs:
	print 'Input pdb files required'
	quit()

if not args.frag_size:
	args.frag_size = [9]

#Load vall structure database
sdb = StructureDatabase("/work/fordas/databases/structure/vall.h5")

vall_residues = sdb.residues.read()
rosetta.init()

lines = []
plt.clf()
summary_header = [ 'pdb_id' ]
summary_out = []
result = {}

# Specify fragment type for extraction
for frag_size in sorted( args.frag_size ):
	nmer = '{0}mer'.format( frag_size )
	summary_header.append( nmer )
	
	test_fragment_type = FragmentSpecification( frag_size, "CA" )
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
			fig = result[pdb_id]['fig'] 
			lines = result[pdb_id]['lines']
			test_structure_pose = result[pdb_id]['pose']
		else:
			# read in pdb
			assert( os.path.isfile( pdb ) )
			test_structure_pose = rosetta.core.import_pose.pose_from_pdb( pdb )

			# initialize dict entry
			pose = test_structure_pose
			fig = plt.figure()
			lines = []
			worst_frag = []
			mean_frag = []
			result[pdb_id] = { 'worst_frag' : worst_frag, 'mean_frag' : mean_frag, 'fig' : fig, 'lines' : lines, 'pose' : pose }

		total_res = test_structure_pose.total_residue()
		test_structure_residues = StructureDatabase.extract_residue_entries_from_pose( test_structure_pose )
	
		test_fragment_start_indices, test_fragments = test_fragment_type.fragments_from_source_residues( test_structure_residues )

		test_to_vall_alignment = atom_array_broadcast_rmsd( all_vall_fragments["coordinates"], test_fragments[:total_res]["coordinates"] )

		min_alignment = test_to_vall_alignment.min(axis=0)
		worst_frag_rmsd = max( min_alignment )
		mean_frag_rmsd = numpy.mean( min_alignment )

		result[pdb_id]['worst_frag'].append( worst_frag_rmsd )
		result[pdb_id]['mean_frag'].append( mean_frag_rmsd )

		if not args.plot_off:
			ax1 = fig.add_subplot( 111 )
			color = next( ax1._get_lines.color_cycle )
			new_ax = ax1.twinx()
			new_ln, = new_ax.plot( min_alignment, color=color, label=nmer )
			#new_ax.set_ylim( 0, 2 * worst_frag_rmsd )
			new_ax.set_ylim( 0, 1.0 )
			new_ax.tick_params( right='off', labelright='off' )
			lines.append( new_ln )

if not args.plot_off:
	for pdb_id in result.keys():
		fig = result[pdb_id]['fig']
		lines = result[pdb_id]['lines']

		labels = [ line.get_label() for line in lines ]
		#plt.figlegend( ( l1, l2 ), ( '9mer', '4mer' ), loc=0 )
		fig.legend( lines, labels, loc=1 )

		#fig.set_xlabel( 'Residue position' )
		#fig.set_ylabel( 'Fragment RMSD' )
		fig.savefig( '{0}_frag.png'.format( pdb_id ) )


# write worst_fragment data
summary_lines = []
for pdb_id in result.keys():
	data = map( lambda x: '%.4f' % x, result[pdb_id]['worst_frag'] )
	data[:0] = [ pdb_id ]
	summary_lines.append( '\t'.join( data ) )

outfile = args.out + '_worst.out'
with open ( outfile, 'w' ) as f:
		f.write( '\t'.join( summary_header ) + '\n' )
		f.write( '\n'.join( summary_lines ) )

# write mean_fragment data
summary_lines = []
for pdb_id in result.keys():
	data = map( lambda x: '%.4f' % x, result[pdb_id]['mean_frag'] )
	data[:0] = [ pdb_id ]
	summary_lines.append( '\t'.join( data ) )

outfile = args.out + '_mean.out'
with open ( outfile, 'w' ) as f:
		f.write( '\t'.join( summary_header ) + '\n' )
		f.write( '\n'.join( summary_lines ) )
