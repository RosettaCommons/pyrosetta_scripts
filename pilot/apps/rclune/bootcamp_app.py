import sys
import argparse
from pyrosetta import*
import numpy as np

# Lab 4
def identify_secondary_structure_spans(ss):
    """
    Takes a string of H's E's and spaces and returns a list describing
    how many secondary structure elements were found and the first and 
    last residues that define each element. 

    :param ss: A string of H's, E's and spaces defining the secondary
    structure of your pose
    :returns: Returns a list where the length of the list is the same
    as the number of secondary structure elements in the function. 
    The list will contain tuples where the first residue in the tuple 
    is the first residue of the SS element, and the second is the
    last residue in the SS element. 
    """

    elements = []
    start = None

    # Rosetta starts counting at 1, sigh
    for ii in range(1, len(ss)+1):
        if len(ss) == 0:
            print("Empty string given.")
            return elements
        
        #if len(ss) == 1:
        #    print("String of length 1 given, no secondary structure.")
        #    return elements
        
        current_char = ss[ii-1]
        if current_char in "EH":
            if start is None:
                start = ii
            if ii == len(ss) or ss[ii] != current_char:
                elements.append((start, ii))
                start = None

    return elements

def get_edges(ss_string):

    ss_elements = identify_secondary_structure_spans(ss_string)

    edges = []
    start = 1

    midpoint0 = (ss_elements[0][0] + ss_elements[0][1])//2
    midpoint1 = None
    midpoint2 = None
    jump_num = 1
    
    #all_vals = []
    #for element in ss_elements: 
    #    all_vals.append(element[0])
    #    all_vals.append(element[1])

    for ii in range(len(ss_elements)-1): 
        midpoint1 = (ss_elements[ii][0] + ss_elements[ii][1]) // 2
        edges.append((midpoint1, start, -1))
        edges.append((midpoint1, ss_elements[ii][1], -1)) 
        
        start = ss_elements[ii][1] + 1
        midpoint2 = (ss_elements[ii][1] + ss_elements[ii+1][0])//2
        # jump_edge
        edges.append((midpoint0, midpoint2, jump_num))
        jump_num += 1
        edges.append((midpoint2, start, -1))
        edges.append((midpoint2, ss_elements[ii+1][0]-1, -1))

        # jump edge
        next_midpoint = (ss_elements[ii+1][0] + ss_elements[ii+1][1]) // 2
        edges.append((midpoint0, next_midpoint, jump_num))
        jump_num += 1
        start = ss_elements[ii+1][0]
        

    # the last element: 
    midpoint = (ss_elements[-1][0] + ss_elements[-1][1])//2
    edges.append((midpoint, start, -1))
    edges.append((midpoint, len(ss_string), -1))
    
    return edges

def fold_tree_from_ss(mypose):
    """

    :param mypose: 
    :returns: A FoldTree
    """
    mydsspmv = rosetta.protocols.moves.DsspMover()
    mydsspmv.apply(mypose)
    ss_string = mypose.secstruct()

    return fold_tree_from_dssp_string(ss_string)


def fold_tree_from_dssp_string(ss_string):
    """
    Takes the string returned by DSSP and creates a FoldTree
    :param ss_string: This string comes from Rosetta/PyRosetta's DSSP
    code. It is a string of H's and E's (and maybe some other letters)
    that defines the secondary structure elements of a given pose.
    :returns: A FoldTree
    """
    myft = rosetta.core.kinematics.FoldTree()

    edges = get_edges(ss_string)

    for edge in edges:
        myft.add_edge(edge[0], edge[1], edge[2])
    
    return myft



def main():
    init(extra_options="-ignore_unrecognized_res")

    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--structure", required=True, help="PDB file to load")
    args = parser.parse_args()

    # 3. Load the Pose from File
    print(f'Input file: {args.structure}')
    mypose = pose_from_pdb(args.structure)
    print(f"Loaded pose with {mypose.total_residue()} residues from: {args.structure}")
    
    myft = fold_tree_from_ss(mypose)

    # 4. Score the Pose
    sfxn = rosetta.core.scoring.get_score_function()

    # lab 4 enabling a new score term linear_chainbreak
    sfxn.set_weight(rosetta.core.scoring.ScoreType.linear_chainbreak, 1)
    
    # add the cutpoint variants
    rosetta.core.pose.correctly_add_cutpoint_variants(mypose)

    myscore = sfxn(mypose)
    print(f"Original score: {myscore}")
    
    # 5. Create Monte Carlo Protocol to Optimize Pose
    
    # C. Initialize a MonteCarloObject and loop over residues
    temperature = 1.0
    mymc = rosetta.protocols.moves.MonteCarlo(mypose, sfxn, temperature)
    
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
    print(f"Temperature = {temperature}")

if __name__ == "__main__":
    main()
