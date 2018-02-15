#!/work/fordas/virtualenvs/fragment_profiling_v0.1/bin/python

#===============================================================================
# Before running, source Alex Ford's virtual environment:
#   source /work/fordas/virtualenvs/fragment_profiling_v0.1/bin/activate
#
# Usage:
#   profile_vall.py file1.pdb [ file2.pdb ... ]
#
#===============================================================================

import os,sys

# apparently my version of numpy is incompatible with Alex's libraries, 
# ... so edit the PYTHONPATH before loading those libraries
#os.environ[ 'PYTHONPATH' ] = '/work/fordas/virtualenvs/close_structure_next/lib/python2.7/site-packages/:' + os.environ[ 'PYTHONPATH' ]
#sys.path.insert( 1, '/work/fordas/virtualenvs/close_structure_next/lib/python2.7/site-packages/' )
sys.path.insert( 1, '/work/fordas/virtualenvs/close_structure_next/src/interface-fragment-matching/' )

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
import matplotlib.pyplot as plt
import argparse

def main( arguments ):

	parser = argparse.ArgumentParser()
	parser.add_argument( 'input_pdbs', type=str, nargs='*', help='input pdb files' )
	parser.add_argument( '-l', type=str, help='list of pdb files' )
	parser.add_argument( '-frag_size', type=int, nargs='*', help='fragment_length' )
	parser.add_argument( '-out', type=str, default='score.dat', help='output filename' )
	parser.add_argument( '-update', type=str, help='update an existing file with new data' )

	args = parser.parse_args( arguments )

	# parse arguments
	if args.l:
		with open( args.l, 'r' ) as f:
			args.input_pdbs.extend( [ line.strip() for line in f.readlines() ] )

	if not args.input_pdbs:
		raise ValueError( 'input pdb files required' )

	if not args.frag_size:
		args.frag_size = [9]

	#Load vall structure database
	sdb = StructureDatabase("/work/fordas/databases/structure/vall.h5")

	vall_residues = sdb.residues.read()
	rosetta.init()

	# if we're updating an existing file, go ahead and read it into memory
	data = {}
	if( args.update ):
		with open( args.update, 'r' ) as f:
			records = [ line.strip().split() for line in f.readlines() ]

			header = records.pop( 0 )
			assert( header[0] == 'pdb_id' )
			header.pop( 0 )

			# make sure the old data table contains the appropriate data
			for frag_size in args.frag_size:
				assert( 'profile_worst_{0}mer'.format( frag_size ) in header )
				assert( 'profile_mean_{0}mer'.format( frag_size ) in header )
		
			# read file into dict
			for record in records:
				pdb_id = record.pop( 0 )
				data[ pdb_id ] = {}
				for field, value in zip( header, record ):
					data[ pdb_id ][ field ] = float( value )

	poses = {}
	for frag_size in sorted( args.frag_size ):
		field_worst = 'profile_worst_{0}mer'.format( frag_size )
		field_mean	= 'profile_mean_{0}mer'.format( frag_size )
		
		test_fragment_type = FragmentSpecification( frag_size, "CA" )
		all_vall_fragment_start_residues, all_vall_fragments = test_fragment_type.fragments_from_source_residues( vall_residues )

		#test_fragment_type_4mer = FragmentSpecification(4, "CA")
		#all_vall_fragment_start_residues_4mer, all_vall_fragments_4mer = test_fragment_type_4mer.fragments_from_source_residues( vall_residues )

		for pdb in args.input_pdbs:
			pdb_id = os.path.splitext( os.path.basename( pdb ) )[0]
				
			# skip for which data already exists
			if not pdb_id in data:
				data[ pdb_id ] = {}

			elif field_worst in data[ pdb_id ] and field_mean in data[ pdb_id ]:
				continue

			# Vall contains ~3800000 9mers
			# 0.001 percentile = top 38 9mers
			# high_percentile_alignment = numpy.percentile(test_to_vall_alignment, .001, axis=0)

			# if we've already loaded this pdb, just get the loaded pose
			if pdb_id in poses:
				test_structure_pose = poses[pdb_id]
	
			else:
				# read in pdb
				assert( os.path.isfile( pdb ) )
				test_structure_pose = rosetta.core.import_pose.pose_from_pdb( pdb )

				# initialize dict entry
				pose = test_structure_pose
				poses[pdb_id] = pose

			total_res = test_structure_pose.total_residue()
			test_structure_residues = StructureDatabase.extract_residue_entries_from_pose( test_structure_pose )

			test_fragment_start_indices, test_fragments = test_fragment_type.fragments_from_source_residues( test_structure_residues )

			test_to_vall_alignment = atom_array_broadcast_rmsd( all_vall_fragments["coordinates"], test_fragments[:total_res]["coordinates"] )

			min_alignment = test_to_vall_alignment.min(axis=0)
			worst_frag_rmsd = max( min_alignment )
			mean_frag_rmsd = numpy.mean( min_alignment )

			data[pdb_id][ field_worst ] = worst_frag_rmsd
			data[pdb_id][ field_mean ] = mean_frag_rmsd


	# incorporate new results into any data that already existed in update file
	#for pdb_id, record in result.iteritems():
	#	if not pdb_id in data:
	#		data[ pdb_id ] = {}
	#	data[ pdb_id ].update( record )
		
	# formulate header
	summary_header = [ 'pdb_id' ]
	for frag_size in sorted( args.frag_size ):
		summary_header.extend( [ '{0}_{1}mer'.format( key, frag_size ) for key in ( 'profile_worst', 'profile_mean' ) ] )

	# formulate body
	summary_lines = [ '\t'.join( summary_header ) ]
	for pdb_id in sorted( data.keys() ):
		out_record = [ pdb_id ]
		for frag_size in sorted( args.frag_size ):
			out_record.extend( [ '{0:5.3f}'.format( data[ pdb_id ][ '{0}_{1}mer'.format( key, frag_size ) ] ) for key in ( 'profile_worst', 'profile_mean' ) ] )
		summary_lines.append( '\t'.join( out_record ) )

	# write output
	with open ( args.out, 'w' ) as f:
		f.write( '\n'.join( summary_lines ) )

if __name__ == '__main__':
	main( sys.argv[1:] ) # omit argv[0], which should just be this script's name
