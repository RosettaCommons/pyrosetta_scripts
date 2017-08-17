import pyrosetta
from pyrosetta.rosetta import protocols 
import argparse 
from Bio.Data import IUPACData 

# parse command line arguments 
parser = argparse.ArgumentParser()
parser.add_argument('--input_pose', help='PDB file of input pose')
parser.add_argument('--mutation', help='The mutation you would like to make in the form A123B')
parser.add_argument('--output_filename', help='Name of the output PDB') 
args = parser.parse_args()

# initialize PyRosetta
pyrosetta.init() 

# create a pose
pose = pyrosetta.pose_from_file(args.input_pose)

# get a scorefxn
score_function = pyrosetta.create_score_function('ref2015') 

# mutate the residue 
target = int(args.mutation[1:-1]) 
new_res = IUPACData.protein_letters_1to3[args.mutation[-1]].upper()
mutate = protocols.simple_moves.MutateResidue(target, new_res)
mutate.apply(pose) 

# output a PDB of the mutation 
pose.dump_pdb(args.output_filename)
