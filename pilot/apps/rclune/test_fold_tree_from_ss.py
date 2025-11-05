import pytest
from pyrosetta import *

from bootcamp_app import identify_secondary_structure_spans, get_edges

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
expected_edges = [(7, 1, -1), (7, 10, -1), (7, 12, 1), (12, 11, -1), (12, 14, -1), (7, 18, 1), (18, 15, -1), (18, 21, -1), (7, 26, 1), (26, 22, -1), (26, 30, -1), (7, 35, 1), (35, 31, -1), (35, 39, -1), (7, 41, 1), (41, 40, -1), (41, 43, -1), (7, 48, 1), (48, 44, -1), (48, 53, -1), (7, 55, 1), (55, 54, -1), (55, 56, -1), (7, 59, 1), (59, 57, -1), (59, 62, -1), (7, 67, 1), (67, 63, -1), (67, 71, -1), (7, 76, 1), (76, 72, -1), (76, 80, -1), (7, 85, 1), (85, 81, -1), (85, 89, -1), (7, 92, 1), (92, 90, -1), (92, 99, -1)]

def test_get_edges():
    assert get_edges(ss) == expected_edges

