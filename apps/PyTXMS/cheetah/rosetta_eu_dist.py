from rosetta import *
from pyrosetta import *
from rosetta.protocols.rigid import *
from math import sqrt, pow, exp

def rosetta_eu_dist(pose, a, b):
    dist = sqrt(pow(( pose.residue( b ).atom( "CB" ).xyz().x - pose.residue( a ).atom( "CB" ).xyz().x ) , 2) \
        + pow(( pose.residue( b ).atom( "CB" ).xyz().y - pose.residue( a ).atom( "CB" ).xyz().y ) , 2) \
	+ pow(( pose.residue( b ).atom( "CB" ).xyz().z - pose.residue( a ).atom( "CB" ).xyz().z ) , 2))
    return dist
