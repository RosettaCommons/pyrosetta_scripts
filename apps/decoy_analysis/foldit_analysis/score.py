#!/home/koepnick/local/anaconda2/bin/python

#from rosetta import *
from pyrosetta import *
import os,sys
import argparse

def main( arguments ):

	parser = argparse.ArgumentParser()
	parser.add_argument( 'input_pdbs', type=str, nargs='*', help='input pdb files' )
	parser.add_argument( '-l', type=str, help='list of pdb files' )
	parser.add_argument( '-out', type=str, default='score.dat', help='output filename' )
	parser.add_argument( '-scorefxn', type=str, default='talaris2013_cart', help='score function name' )
	parser.add_argument( '-scorefxn_foldit', action='store_true', help='use foldit_design scorefxn and cmdline options (overrides -scorefxn)' )
	parser.add_argument( '-resmax', action='store_true', help='also report maximum single-residue score components' )
	parser.add_argument( '-update', type=str, help='update an existing file with new data' )

	args = parser.parse_args( arguments )

	# parse arguments
	if args.l:
		with open( args.l, 'r' ) as f:
			args.input_pdbs.extend( [ line.strip() for line in f.readlines() ] )

	if not args.input_pdbs:
		raise ValueError( 'input pdb files required' )

	if( args.scorefxn_foldit ):
		init( '-envsmooth_zero_negatives', '-rama_power 3' )
		args.scorefxn = 'talaris2013_cart'
		scorefxn = create_score_function( args.scorefxn )
		scorefxn.set_weight( rosetta.core.scoring.score_type_from_name( 'cart_bonded' ), 2.0 ) 
		scorefxn.set_weight( rosetta.core.scoring.score_type_from_name( 'envsmooth' ), 2.0 ) 
	else:
		init()
		scorefxn = create_score_function( args.scorefxn )

	# if we're updating an existing data table, go ahead and read it into memory
	data = {}
	if( args.update ):
		with open( args.update, 'r' ) as f:
			records = [ line.strip().split() for line in f.readlines() ]

			header = records.pop( 0 )
			assert( header[0] == 'pdb_id' )
			header.pop( 0 )
		
			# read file into dict
			for record in records:
				pdb_id = record.pop( 0 )
				new_record = {}

				# try/catch block in case there is a problem reading in previous data
				try:
					#data[ pdb_id ] = {}
					for field, value in zip( header, record ):
						#data[ pdb_id ][ field ] = float( value )
						new_record[ field ] = float( value )
					data[ pdb_id ] = new_record
				except Exception:
					pass

	# iterate over input pdbs
	scoretypes = scorefxn.get_nonzero_weighted_scoretypes()
	for pdb in args.input_pdbs:
		pdb_id = os.path.splitext( os.path.basename( pdb ) )[0]

		# if we already have data for this pdb_id, skip
		if pdb_id in data:
			continue
		data[ pdb_id ] = {}

		pose = rosetta.core.import_pose.pose_from_file( pdb )
		nres = pose.total_residue()
		total_score = scorefxn( pose )

		data[ pdb_id ][ 'nres' ] = nres
		for scoretype in scoretypes:
			score_name = rosetta.core.scoring.name_from_score_type( scoretype )
			data[ pdb_id ][ score_name + '_res' ] = pose.energies().total_energies().get( scoretype ) / float( nres )
			data[ pdb_id ][ score_name + '_resmax' ] = max( [ pose.energies().residue_total_energies( ii )[ scoretype ] for ii in range( 1, nres+1 ) ] )

		# also get total score
		data[ pdb_id ][ 'score_res' ] = total_score / float( nres )
		data[ pdb_id ][ 'score_resmax' ] = max( [ pose.energies().residue_total_energy( ii ) for ii in range( 1, nres+1 ) ] )

	# formulate header
	fields = [ 'nres', 'score_res', 'score_resmax' ]
	fields.extend( [ rosetta.core.scoring.name_from_score_type( s ) + '_res' for s in scoretypes ] )
	if args.resmax:
		fields.extend( [ rosetta.core.scoring.name_from_score_type( s ) + '_resmax' for s in scoretypes ] )
	header = [ 'pdb_id' ] + fields
	outlines = [ '\t'.join( header ) ]

	# formulate body
	for pdb_id, record in sorted( data.iteritems() ):
		out_record = [ pdb_id ]	
		for field in fields:
			if field in record:
				out_record.append( '{0:5.3f}'.format( record[ field ] ) )
				#out_record.extend( [ '{0:5.3f}'.format( record[ stype ] ) for stype in fields ] )
			else:
				out_record.append( 'NA' )
		outlines.append( '\t'.join( out_record ) )

	# write output
	with open ( args.out, 'w' ) as f:
		f.write( '\n'.join( outlines ) )

if __name__ == '__main__':
	main( sys.argv[1:] ) # omit argv[0], which should just be this script's name
