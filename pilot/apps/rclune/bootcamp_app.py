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
mymc = rosetta.protocols.moves.MonteCarlo(mypose, sfxn, 1.0)

# 7. MoveMap stuff 
mm = rosetta.core.kinematics.MoveMap()
mm.set_chi(True)
mm.set_bb(True)

min_opts = rosetta.core.optimization.MinimizerOptions("lbfgs_armijo_atol", 0.01, True)
minimizer = rosetta.core.optimization.AtomTreeMinimizer()

# 7 Packing and Minimizing
# Packing
tf = rosetta.core.pack.task.TaskFactory()
task = tf.create_task_and_apply_taskoperations(mypose)
task.restrict_to_repacking()
rosetta.core.pack.pack_rotamers(mypose, sfxn, task)

minimizer.run(mypose, mm, sfxn, min_opts)

# 6 PyMOL
the_observer = rosetta.protocols.moves.PyMOLObserver()
the_observer.pymol().apply(mypose)

acceptance_rate = 0
cumulative_acceptance_rate = 0
print_rate = 100

for ii in range(1000):
    print(f"Step {ii}: ")

    # A. get your random residue and random perturbations
    myres = np.random.randint(1, mypose.total_residue()+1)
    my_phi = np.random.normal()
    my_psi = np.random.normal()

    # B. perturb the pose, paying attention to if the residue is not an amino acid
    try:
        orig_phi = mypose.phi(myres)
        orig_psi = mypose.psi(myres)
        mypose.set_phi(myres, orig_phi + my_phi)
        mypose.set_psi(myres, orig_psi + my_psi)

    except:
        print(f"Residue {myres} is not an amino acid")
        pass
    
    rosetta.core.pack.pack_rotamers(mypose, sfxn, task)
    minimizer.run(mypose, mm, sfxn, min_opts)

    is_accepted = mymc.boltzmann(mypose)
    print(is_accepted)
    if is_accepted:
        acceptance_rate += 1
        cumulative_acceptance_rate += 1
    
    the_observer.pymol().apply(mypose)

    if ii % print_rate == 0 and ii != 0:
        print(f"Acceptance rate: {acceptance_rate/print_rate}")
        acceptance_rate = 0 # reset acceptance rate to 0
        print(f"Cumulative acceptance rate: {cumulative_acceptance_rate/ii}")
        print(mymc.show_counters())
        print(f"Average score: {mypose.energies().total_energy()}")



print(mymc.last_accept())
final_pose = mymc.lowest_score_pose()
final_score = mymc.lowest_score()
print(f"Final Score: {final_score}")
final_pose.dump_pdb(f"optimized_{args.structure}")