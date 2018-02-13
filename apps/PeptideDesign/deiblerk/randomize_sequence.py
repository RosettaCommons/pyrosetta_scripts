import argparse
import pyrosetta
from pyrosetta import *
import random
import sys 

"""Do I need import sys....yes that is for arg statements"""

"""Attempting to write a code to randomize a sequence of a linear peptide to generate and library random lirary best on a set structure with help from Jacob O'Connor."""
def mutate_peptide (native_pdb):
    peptide = pose_from_pdb(native_pdb)
    
    list = [ "ALA", "GLY", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "VAL", "SER", "THR", "ASN", "GLN", "CYS", "PRO", "ARG", "HIS", "LYS", "ASP", "GLU", "DALA", "DILE", "DLEU", "DMET", "DPHE", "DTRP", "DTYR", "DVAL", "DSER", "DTHR", "DASN", "DGLN", "DCYS", "DPRO", "DARG", "DHIS", "DLYS", "DASP", "DGLU"]
    
    for residue in range (1, len(peptide.sequence())+1):
        randomresidue_number = random.randint(0,len(list)-1)
        #randomresidue_object = pyrosetta.rosetta.numeric.random.RandomGenerator()
        #randomresidue_number = randomresidue_object.random_range(0, len(list)-1)
        mutation = list[randomresidue_number]
        mutateresidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue(residue, mutation)
        mutateresidue.apply(peptide)
    return peptide

def main(argv):

    parser = argparse.ArgumentParser(description="Reprocess boinc output", fromfile_prefix_chars='@')
    parser.add_argument('-native_pdb', action='store', type=str, required=True, help="Name the original pdb")
    parser.add_argument('-pose_outputs', action='store', type=int, required=True, help="How many output poses desired?")
    args = parser.parse_args()

    
    init()
    sequence = []
    for outputs in range (1, args.pose_outputs+1):
        random_pose = mutate_peptide(args.native_pdb)
        unique = 0
        for thing in sequence:
            if thing == random_pose.sequence:
                unique = 1
        if unique == 1:
            continue
        sequence.append(random_pose.sequence())
        random_pose.dump_pdb ("random_%s.pdb" %(outputs))
        #silentfileoptions = pyrosetta.rosetta.core.io.silent.SilentFileOptions()
        #silentfilestruc = pyrosetta.rosetta.core.io.silent.SilentStruct(silentfileoptions)
        #silentfilestruc.fill_struct(random_pose, "random_%s" %(outputs))
        #silentfiledata = pyrosetta.rosetta.core.io.silent.SilentFileData("random.out", 0, 0, "binary", silentfileoptions)
        #silentfiledata.write_silent_struct(silentfilestruc)

if __name__ == '__main__':
    main(sys.argv)

