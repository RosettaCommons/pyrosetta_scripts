import pytest

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
    print(len(ss))

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

ss1 = "   EEEEE   HHHHHHHH  EEEEE   IGNOR EEEEEE   HHHHHHHHHHH  EEEEE  HHHH   "
expected1 = [(4, 8), (12, 19), (22, 26), (36, 41), (45, 55), (58, 62), (65, 68)]

ss2 = "HHHHHHH   HHHHHHHHHHHH      HHHHHHHHHHHHEEEEEEEEEEHHHHHHH EEEEHHH "
expected2 = [(1, 7), (11, 22), (29, 40), (41, 50), (51, 57), (59, 62), (63, 65)]

ss3 = "EEEEEEEEE EEEEEEEE EEEEEEEEE H EEEEE H H H EEEEEEEE"
expected3 = [(1,9), (11, 18), (20, 28), (30, 30), (32, 36), (38, 38), (40, 40), (42, 42), (44, 51)]

ss4 = ""
expected4 = []

ss5 = "E"
expected5 = [(1,1)]

ss6 = "EH"
expected6 = [(1,1),(2,2)]

@pytest.mark.parametrize("value,expected", 
                          [
                              (ss1, expected1),
                              (ss2, expected2),
                              (ss3, expected3),
                              (ss4, expected4),
                              (ss5, expected5),
                              (ss6, expected6)
                          ])
def test_identify_secondary_structure_spans(value, expected):
    assert identify_secondary_structure_spans(value) == expected
    

