import pytest
from pyrosetta import *

from bootcamp_app import identify_secondary_structure_spans, get_edges, fold_tree_from_dssp_string

init(extra_options="-ignore_unrecognized_res")

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

ss7 = "EEEEEESSSDFSDHH    HHHEEEEISUDHHHHHHEE"
expected7 = [(1,6),(14,15),(20,22),(23,26),(31,36),(37,38)]

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
    
ss = "   EEEEEEE    EEEEEEE         EEEEEEEEE    EEEEEEEEEE   HHHHHH         EEEEEEEEE         EEEEE     "
expected_edges = [(7, 1, -1), (7, 10, -1), (7, 12, 1), (12, 11, -1), (12, 14, -1), (7, 18, 2), (18, 15, -1), (18, 21, -1), (7, 26, 3), (26, 22, -1), (26, 30, -1), (7, 35, 4), (35, 31, -1), (35, 39, -1), (7, 41, 5), (41, 40, -1), (41, 43, -1), (7, 48, 6), (48, 44, -1), (48, 53, -1), (7, 55, 7), (55, 54, -1), (55, 56, -1), (7, 59, 8), (59, 57, -1), (59, 62, -1), (7, 67, 9), (67, 63, -1), (67, 71, -1), (7, 76, 10), (76, 72, -1), (76, 80, -1), (7, 85, 11), (85, 81, -1), (85, 89, -1), (7, 92, 12), (92, 90, -1), (92, 99, -1)]

def test_get_edges():
    assert get_edges(ss) == expected_edges

def test_num_edges():
    assert len(get_edges(ss)) == 4*7-2 + 2*7-2

def test_fold_tree_from_dssp_string_valid():
    myft=fold_tree_from_dssp_string(ss)
    assert myft.check_fold_tree() == True

expected_ft_string = "FOLD_TREE  EDGE 7 1 -1  EDGE 7 10 -1  EDGE 7 12 1  EDGE 12 11 -1  EDGE 12 14 -1  EDGE 7 18 2  EDGE 18 15 -1  EDGE 18 21 -1  EDGE 7 26 3  EDGE 26 22 -1  EDGE 26 30 -1  EDGE 7 35 4  EDGE 35 31 -1  EDGE 35 39 -1  EDGE 7 41 5  EDGE 41 40 -1  EDGE 41 43 -1  EDGE 7 48 6  EDGE 48 44 -1  EDGE 48 53 -1  EDGE 7 55 7  EDGE 55 54 -1  EDGE 55 56 -1  EDGE 7 59 8  EDGE 59 57 -1  EDGE 59 62 -1  EDGE 7 67 9  EDGE 67 63 -1  EDGE 67 71 -1  EDGE 7 76 10  EDGE 76 72 -1  EDGE 76 80 -1  EDGE 7 85 11  EDGE 85 81 -1  EDGE 85 89 -1  EDGE 7 92 12  EDGE 92 90 -1  EDGE 92 99 -1 "
def test_fold_tree_from_dssp_string_edges():
    myft = fold_tree_from_dssp_string(ss)
    assert myft.to_string() == expected_ft_string

