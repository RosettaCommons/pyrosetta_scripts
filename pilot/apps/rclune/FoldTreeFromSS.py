from pyrosetta import *
import numpy as np

from dataclasses import dataclass

@dataclass
class _Built:
    ft: rosetta.core.kinematics.FoldTree
    loops: List[rosetta.protocols.loops.Loop]
    loop_for_residue: List[int] 

class FoldTreeFromSS: 
    def __init__(self, 
                 pose: rosetta.core.pose.Pose, 
                 loop_left: int = 2,
                 loop_right: int = 3):
        """
        :param pose: The pose object for your input structure
        :param loop_left: How many residues you perturb to the left
        of your cutpoint to close it
        :param loop_right: How many residues you perturb to the right
        of your cutpoint to close it
        """
        self._pose = pose
        self._loop_left = loop_left
        self._loop_right = loop_right
        
        mydsspmv = rosetta.protocols.moves.DsspMover()
        mydsspmv.apply(self._pose)
        self._ss_string = self._pose.secstruct()

        self._ss_elements = self.identify_secondary_structure_spans()
        self._cutpoints = self.calculate_cutpoints()
        self._loops = []
        self._loop_for_residue=[]
        self._ft = rosetta.core.kinematics.FoldTree()
        #self._loop_data = _Built(ft=None, loops=[], loop_for_residue=[])

        self.create_loop_list()
        self.create_loop_for_residue_list()
        self.fold_tree_from_ss()
        
    def create_loop_list(self):
        # Just need to count the number of loops
        # I am assuming there is a loop between each ss element
        
        cutpoints = self.calculate_cutpoints()

        # index of 0 means that no loop closure is needed
        # so I'm having the 0th element not be a loop!
        self._loops.append(None)

        for cutpoint in cutpoints:
            self._loops.append(rosetta.protocols.loops.Loop(cutpoint-self._loop_left, cutpoint+self._loop_right, cutpoint))
            
    def create_loop_for_residue_list(self):
        """
        If there is no loop to close (it's in the first or last peptide
        edge) then it should be 0. Otherwise it should follow the loop
        to the nearest cutpoint.
        """

        midpoints = []

        for ii in range(len(self._ss_elements)-1):
            midpoints.append((self._ss_elements[ii][0] + self._ss_elements[ii][1])//2)
            midpoints.append((self._ss_elements[ii][1] + self._ss_elements[ii+1][0])//2)
        midpoints.append((self._ss_elements[-1][0] + self._ss_elements[-1][1])//2)

        # if it's before the first midpoint then it should return 0
        # if it's after the last midpoint then it should return 0

        start = 0
        for i, end in enumerate(midpoints):
            self._loop_for_residue.extend([i] * (end - start-1))
            start = end-1 

        # add zeros to the remaining portion of the list: 
        self._loop_for_residue.extend([0] * (self._pose.total_residue()-start))


    
    def get_fold_tree(self) -> rosetta.core.kinematics.FoldTree:
        """
        Getter for the FoldTree
        """
        #self.fold_tree_from_ss()
        return self._ft
    
    def get_loop(self, index):
        #self.create_loop_list()
        return self._loops[index]
    
    def get_index_from_loop_for_residue(self, index):
        #self.create_loop_for_residue_list()
        return self._loop_for_residue[index]
    
    #def loop(self, index: int) -> rosetta.protocols.loops.Loop:
    #    """
    #    This function takes the loop index and returns the actual loops
    #    """
    #
    #    return self._residues_for_loop[index]
    
    def calculate_cutpoints(self):
        """
        Don't know if I'll need this but I feel like it'll be useful.
        Create a function that just creates a list of the indices of
        the residues that are to the left of the cutpoint.
        """
        cutpoints = []

        # first cutpoint is the end of the first ss element
        cutpoints.append(self._ss_elements[0][1])

        for ii in range(1, len(self._ss_elements)-1):
            cutpoints.append(self._ss_elements[ii][0]-1)
            cutpoints.append(self._ss_elements[ii][1])

        # last cutpoint is at the beginning of the last element
        cutpoints.append(self._ss_elements[-1][0]-1)

        print("CUTPOINTS")
        print(cutpoints)

        return cutpoints



    def identify_secondary_structure_spans(self):
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
        print(self._ss_string)

        start = None
    
        # Rosetta starts counting at 1, sigh
        for ii in range(1, len(self._ss_string)+1):
            if len(self._ss_string) == 0:
                print("Empty string given.")
                return elements
            
            # taken from the string_splitter homeworks
            # Checks if the character is E or H and then if it's different
            # from the previous character, if there is one
            current_char = self._ss_string[ii-1]
            if current_char in "EH":
                if start is None:
                    start = ii
                if ii == len(self._ss_string) or self._ss_string[ii] != current_char:
                    elements.append((start, ii))
                    start = None
        return elements
    
    def get_edges(self):
        """
        Function to turn a string of H's and E's, like what will be returned
        from the DSSP function, into edges for a FoldTree. 
        Edges are the form (start, end, n), the last value should be -1 if
        it is a peptide edge and an integer from 1 to the number of jumps
        if it is a jump point. 
        :param ss: the secondary structure string output by DSSP
        :return: a list of tuples storing three values each
        """
    
        #self._ss_elements = identify_secondary_structure_spans(self._ss_string)
    
        edges = []
        start = 1
    
        # this is the midpoint that all the jump points will be starting
        # from
        midpoint0 = (self._ss_elements[0][0] + self._ss_elements[0][1])//2
        midpoint1 = None
        midpoint2 = None
        jump_num = 1 # used to index the jump points
    
    
        for ii in range(len(self._ss_elements)-1): 
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
            midpoint1 = (self._ss_elements[ii][0] + self._ss_elements[ii][1]) // 2
            edges.append((midpoint1, start, -1))
            # peptide_edge pointing towards the C terminus (->)
            edges.append((midpoint1, self._ss_elements[ii][1], -1)) 
            
    
            start = self._ss_elements[ii][1] + 1 # update the start 
            # this midpoint is in a region that does not have an ss element
            midpoint2 = (self._ss_elements[ii][1] + self._ss_elements[ii+1][0])//2
    
            # jump_edge
            edges.append((midpoint0, midpoint2, jump_num))
            jump_num += 1 # update jump numbering
    
            # peptide edge pointing towards the N terminus (<-)
            edges.append((midpoint2, start, -1))
            # peptide edge pointing towards the C terminus (->)
            edges.append((midpoint2, self._ss_elements[ii+1][0]-1, -1))
    
            # jump edge (had to do this here instead of waiting for the next
            # loop because it otherwise would have added a jump edge as the
            # first edge, which I don't want. There is probably a cleaner 
            # way to write this algorithm. Will feed into an LLM if I have
            # time/motivation later.)
            next_midpoint = (self._ss_elements[ii+1][0] + self._ss_elements[ii+1][1]) // 2
            edges.append((midpoint0, next_midpoint, jump_num))
            jump_num += 1
            start = self._ss_elements[ii+1][0]
            
    
        # the last element: (treated differently because the end needs to be
        # the length of the ss string)
        midpoint = (self._ss_elements[-1][0] + self._ss_elements[-1][1])//2
        edges.append((midpoint, start, -1))
        edges.append((midpoint, len(self._ss_string), -1))
        
        return edges
    
    def fold_tree_from_ss(self):
        """
        Takes a pose and returns a fold tree. 
        :param mypose: A rosetta.core.pose.Pose (https://graylab.jhu.edu/PyRosetta.documentation/pyrosetta.rosetta.core.pose.html#pyrosetta.rosetta.core.pose.Pose)
        object
        :returns: A FoldTree (https://graylab.jhu.edu/PyRosetta.documentation/pyrosetta.rosetta.core.kinematics.html#pyrosetta.rosetta.core.kinematics.FoldTree)
        """

        mydsspmv = rosetta.protocols.moves.DsspMover()
        mydsspmv.apply(self._pose)
        #self._ss_string = self._pose.secstruct()
    
        return self.fold_tree_from_dssp_string()
    
    def fold_tree_from_dssp_string(self):
        """
        Takes the string returned by DSSP and creates a FoldTree
        :param ss_string: This string comes from Rosetta/PyRosetta's DSSP
        code. It is a string of H's and E's (and maybe some other letters)
        that defines the secondary structure elements of a given pose.
        :returns: A FoldTree object (https://graylab.jhu.edu/PyRosetta.documentation/pyrosetta.rosetta.core.kinematics.html#pyrosetta.rosetta.core.kinematics.FoldTree)
        """
    
        edges = self.get_edges()
    
        for edge in edges:
            self._ft.add_edge(edge[0], edge[1], edge[2])
        
        