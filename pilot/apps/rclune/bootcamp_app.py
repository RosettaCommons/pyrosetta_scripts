import sys
import argparse
from pyrosetta import*
import numpy as np

init(extra_options="-ignore_unrecognized_res")

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--structure", required=True, help="PDB file to load")
args = parser.parse_args()

# 3. Load the Pose from File
print(f'Input file: {args.structure}')
mypose = pose_from_pdb(args.structure)
print(f"Loaded pose with {mypose.total_residue()} residues from: {args.structure}")

# 4. Score the Pose
sfxn = rosetta.core.scoring.get_score_function()
myscore = sfxn(mypose)
print(f"Original score: {myscore}")

# 5. Create Monte Carlo Protocol to Optimize Pose

# C. Initialize a MonteCarloObject and loop over residues
mymc = rosetta.protocols.moves.MonteCarlo(mypose, sfxn, 0.5)

for ii in range(1000):

    # A. get your random residue and random perturbations
    myres = np.random.randint(1, mypose.total_residue()+1)
    my_phi = np.random.normal()
    my_psi = np.random.normal()
    print(sfxn(mypose))

    # B. perturb the pose, paying attention to if the residue is not an amino acid
    try:
        orig_phi = mypose.phi(myres)
        orig_psi = mypose.psi(myres)
        mypose.set_phi(myres, orig_phi + my_phi)
        mypose.set_psi(myres, orig_psi + my_psi)

    except:
        print(f"Residue {myres} is not an amino acid")
        pass
    
    mymc.boltzmann(mypose)


print(mymc.last_accept())
final_pose = mymc.lowest_score_pose()
final_score = mymc.lowest_score()
print(f"Final Score: {final_score}")