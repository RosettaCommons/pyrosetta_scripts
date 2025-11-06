"""
2025 Rosetta Bootcamp++ Lab 6
Creating a Mover based off of the code in Lab 4 (FoldTree + MonteCarlo)

(c) Copyright Rosetta Commons Member Institutions. 
(c) This file is part of the Rosetta software suite and is made available under license. 
(c) The Rosetta software is developed by the contributing members of the Rosetta Commons. 
(c) For more information, see http://www.rosettacommons.org. Questions about this can be 
(c) addressed to University of Washington CoMotion, email: license@uw.edu.

Author:  Rachel Clune
"""

from pyrosetta import *
import numpy as np

class BootCampMover(rosetta.protocols.moves.Mover):
    def __init__(self, sfxn=rosetta.core.scoring.get_score_function(), num_iterations=1000):
        super().__init__(self) 
        self._sfxn = sfxn
        self._num_iterations = num_iterations

    def apply(self, mypose):
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
                
                # taken from the string_splitter homeworks
                # Checks if the character is E or H and then if it's different
                # from the previous character, if there is one
                current_char = ss[ii-1]
                if current_char in "EH":
                    if start is None:
                        start = ii
                    if ii == len(ss) or ss[ii] != current_char:
                        elements.append((start, ii))
                        start = None
            
            return elements
        
        def get_edges(ss_string):
            """
            Function to turn a string of H's and E's, like what will be returned
            from the DSSP function, into edges for a FoldTree. 
            Edges are the form (start, end, n), the last value should be -1 if
            it is a peptide edge and an integer from 1 to the number of jumps
            if it is a jump point. 
            :param ss: the secondary structure string output by DSSP
            :return: a list of tuples storing three values each
            """
        
            ss_elements = identify_secondary_structure_spans(ss_string)
        
            edges = []
            start = 1
        
            # this is the midpoint that all the jump points will be starting
            # from
            midpoint0 = (ss_elements[0][0] + ss_elements[0][1])//2
            midpoint1 = None
            midpoint2 = None
            jump_num = 1 # used to index the jump points
        
        
            for ii in range(len(ss_elements)-1): 
                # NOTE: This code may have an issue if two ss structure
                # elements are directly next to each other (no spaces in
                # between)
                # Yeah, it tries to recreate an edge that already exists. 
                # Not worth fixing in this lab IMO, just note its an issue
                # probably need to add a condition about whether or not it's
                # in a ss element or not
        
                # Or I could change everything so that the ss_elements is just
                # an array, instead of an array of tuples, add 1 to the front
                # and the length of the structure to the end and then loop
                # over each pair of numbers instead of each ss element. 
                # Would need to figure out how to treat the jump edges 
                # appropriately
        
                # I treat the last secondary structure element separately 
                # because it has a different end than the rest and no jump
                # points based on how I set everything up
        
                # the start/end for each edge is the midpoint of the element
                # or the midpoint of the region between the elements, unless
                # its a starting edge, so start needs to be 1, or the last 
                # edge, in which case end needs to be the same as the number
                # of residues, or the length of the string from DSSP
        
                # See image in Lab 4, it will help clarify this
        
                # peptide edge pointing towards the N terminus (<-) 
                # midpoint is in a region that has a ss element
                midpoint1 = (ss_elements[ii][0] + ss_elements[ii][1]) // 2
                edges.append((midpoint1, start, -1))
                # peptide_edge pointing towards the C terminus (->)
                edges.append((midpoint1, ss_elements[ii][1], -1)) 
                
        
                start = ss_elements[ii][1] + 1 # update the start 
                # this midpoint is in a region that does not have an ss element
                midpoint2 = (ss_elements[ii][1] + ss_elements[ii+1][0])//2
        
                # jump_edge
                edges.append((midpoint0, midpoint2, jump_num))
                jump_num += 1 # update jump numbering
        
                # peptide edge pointing towards the N terminus (<-)
                edges.append((midpoint2, start, -1))
                # peptide edge pointing towards the C terminus (->)
                edges.append((midpoint2, ss_elements[ii+1][0]-1, -1))
        
                # jump edge (had to do this here instead of waiting for the next
                # loop because it otherwise would have added a jump edge as the
                # first edge, which I don't want. There is probably a cleaner 
                # way to write this algorithm. Will feed into an LLM if I have
                # time/motivation later.)
                next_midpoint = (ss_elements[ii+1][0] + ss_elements[ii+1][1]) // 2
                edges.append((midpoint0, next_midpoint, jump_num))
                jump_num += 1
                start = ss_elements[ii+1][0]
                
        
            # the last element: (treated differently because the end needs to be
            # the length of the ss string)
            midpoint = (ss_elements[-1][0] + ss_elements[-1][1])//2
            edges.append((midpoint, start, -1))
            edges.append((midpoint, len(ss_string), -1))
            
            return edges

        def fold_tree_from_ss(pose):
            """
            Takes a pose and returns a fold tree. 
            :param pose: A rosetta.core.pose.Pose (https://graylab.jhu.edu/PyRosetta.documentation/pyrosetta.rosetta.core.pose.html#pyrosetta.rosetta.core.pose.Pose)
            object
            :returns: A FoldTree (https://graylab.jhu.edu/PyRosetta.documentation/pyrosetta.rosetta.core.kinematics.html#pyrosetta.rosetta.core.kinematics.FoldTree)
            """
            mydsspmv = rosetta.protocols.moves.DsspMover()
            mydsspmv.apply(pose)
            ss_string = mypose.secstruct()
        
            return fold_tree_from_dssp_string(ss_string)

        def fold_tree_from_dssp_string(ss_string):
            """
            Takes the string returned by DSSP and creates a FoldTree
            :param ss_string: This string comes from Rosetta/PyRosetta's DSSP
            code. It is a string of H's and E's (and maybe some other letters)
            that defines the secondary structure elements of a given pose.
            :returns: A FoldTree object (https://graylab.jhu.edu/PyRosetta.documentation/pyrosetta.rosetta.core.kinematics.html#pyrosetta.rosetta.core.kinematics.FoldTree)
            """
            myft = rosetta.core.kinematics.FoldTree()
        
            edges = get_edges(ss_string)
        
            for edge in edges:
                myft.add_edge(edge[0], edge[1], edge[2])
            
            return myft

        myft = fold_tree_from_ss(mypose)
        
        # Lab 2, 4. Score the Pose
        sfxn = self.get_sfxn()
        
        # lab 4 enabling a new score term linear_chainbreak
        sfxn.set_weight(rosetta.core.scoring.ScoreType.linear_chainbreak, 1)
        
        # Lab 4, add the cutpoint variants
        rosetta.core.pose.correctly_add_cutpoint_variants(mypose)
        
        myscore = sfxn(mypose)
        print(f"Original score: {myscore}")
        
        # Lab 2, 5. Create Monte Carlo Protocol to Optimize Pose
        
        # Lab 2, 5C. Initialize a MonteCarloObject and loop over residues
        temperature = 1.0
        mymc = rosetta.protocols.moves.MonteCarlo(mypose, sfxn, temperature)
        
        # Lab 2, 7. MoveMap stuff 
        mm = rosetta.core.kinematics.MoveMap()
        mm.set_chi(True)
        mm.set_bb(True)
        
        min_opts = rosetta.core.optimization.MinimizerOptions("lbfgs_armijo_atol", 0.01, True)
        minimizer = rosetta.core.optimization.AtomTreeMinimizer()
        
        # Lab2, 7 Packing and Minimizing
        # Packing
        tf = rosetta.core.pack.task.TaskFactory()
        task = tf.create_task_and_apply_taskoperations(mypose)
        task.restrict_to_repacking()
        rosetta.core.pack.pack_rotamers(mypose, sfxn, task)
        
        minimizer.run(mypose, mm, sfxn, min_opts)
        
        acceptance_rate = 0
        cumulative_acceptance_rate = 0
        
        # how often to print information about the MC acceptance rates
        print_rate = 100
        
        for ii in range(self.get_num_iterations()):
            print(f"Step {ii}: ")
        
            # Lab2, 5A. get your random residue and random perturbations
            myres = np.random.randint(1, mypose.total_residue()+1)
            my_phi = np.random.normal()
            my_psi = np.random.normal()
        
            # Lab 2, 5B. perturb the pose, paying attention to if the residue is not an amino acid
            try:
                orig_phi = mypose.phi(myres)
                orig_psi = mypose.psi(myres)
                mypose.set_phi(myres, orig_phi + my_phi)
                mypose.set_psi(myres, orig_psi + my_psi)
        
            except:
                print(f"Residue {myres} is not an amino acid")
                pass
            
            # pack and minimize the rotamers after the MC step
            # Probably only need to do this if the MC step was accepted?
            rosetta.core.pack.pack_rotamers(mypose, sfxn, task)
            minimizer.run(mypose, mm, sfxn, min_opts)
        
            # Run the metropolis algorithm for acceptance
            is_accepted = mymc.boltzmann(mypose)
            #print(is_accepted)
        
            # collect information on acceptance rates
            if is_accepted:
                acceptance_rate += 1
                cumulative_acceptance_rate += 1
        
            # print the acceptance criteria
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

    def get_name():
        return self.__class__.__name__
    
    @staticmethod
    def mover_name(self):
        return self.__class__.__name__
    
    def provide_xml_schema(self, xsd):
        pass
    
    def set_sfxn(self, new_sfxn):
        self._sfxn = new_sfxn
    
    def get_sfxn(self):
        return self._sfxn
    
    def set_num_iterations(self, new_num_iterations):
        self._num_iterations = new_num_iterations
    
    def get_num_iterations(self):
        return self._num_iterations
    
    def parse_my_tag(self, tag, datamap):
        if tag.hasOption("num_iterations"):
            iters = tag.get_option_int("num_iterations", 1)
            self.set_num_iterations(iters)

        if tag.hasOption("sfxn"):
            sfxn_choice = tag.get_option_string("sfxn", "commandline")
            mysfxn = rosetta.core.scoring.parse_score_function(tag, sfxn_choice, datamap)
            self.set_sfxn(mysfxn)

    
