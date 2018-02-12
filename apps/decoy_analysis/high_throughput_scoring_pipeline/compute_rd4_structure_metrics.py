#!/home/haddox/.conda/envs/1711_term_project/bin/python

"""
Compute scores for structures from the Rocklin et al. study

This script computes scores for all the structures from the fourth round of
designs from the Rocklin et al. study ('./data/Rocklin_rd4_PDB_files/*.pdb').
It uses `jug` to parallelize this computation.
"""

import logging
import glob
from os import path
import os

import pyrosetta
import pyrosetta.distributed.io as io

import jug
import pandas

from scripts import combined_scoring

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s.%(msecs).03d %(name)s %(message)s",
    datefmt='%Y-%m-%dT%H:%M:%S'
)

jug.set_jugdir("compute_rd4_structure_metrics.jugdata")

score_pdb_file = jug.TaskGenerator(combined_scoring.score_pdb_file)

@jug.TaskGenerator
def write_results(results_file, sub_results):
    # Merge the scores of all PDBs into a single dataframe
    scores_all_pdbs = pandas.concat(sub_results)

    # Write the results to an outfile
    if not os.path.isdir(os.path.dirname(results_file)):
        os.makedirs(os.path.dirname(results_file))
    scores_all_pdbs.to_csv(results_file)

# Define input PDB file names and the output scoring file
input_pdb_files = glob.glob('./data/Rocklin_rd4_PDB_files/*.pdb')
output_csv = './results/compute_rd4_structure_metrics/scores.csv'

# Score the PDB files
write_results(
    output_csv,
    [ score_pdb_file(f) for f in input_pdb_files ]
)
