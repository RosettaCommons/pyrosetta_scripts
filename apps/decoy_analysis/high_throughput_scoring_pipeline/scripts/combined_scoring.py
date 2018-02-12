"""
Python module for scoring PDBs with both Rosetta scripts and external Python scripts

This script uses `jug` to parallelize the scoring of multiple PDB files at once.

See the function `score_pdb_file` for more details
"""

import logging
import glob
import os
from os import path
import subprocess

import pyrosetta
import pyrosetta.distributed.io as io

import pandas

from scripts import rs_scoring
from scripts import external_scoring
from scripts import xmls
from xmls import analysis_scripts


def score_pdb_file(input_pdb_file):
    """
    Score an input PDB file with RosettaScripts and with external Python scripts

    This function scores a single input PDB file using two modules: one which
    scores the structure using RosettaScripts (`rs_scoring`) and one which
    scores the structure using a pipeline of custom Python scripts
    (`external_scoring`).

    Args:
        input_pdb_file - the path of an input pdb file (string)

    Returns:
        A dataframe with the scores for the input PDB
    """

    logging.info("generate_pdb_score_summary: %s", input_pdb_file)

    # Conver the input PDB to a pose in a dictionary
    input_pose = {input_pdb_file : io.pose_from_file(input_pdb_file)}

    # Compute RosettaScript scores
    # imported analysis scripts above
    rs_scores = rs_scoring.generate_score_summary_table(input_pose, analysis_scripts)

    # Compute external scores
    scriptsdir = path.dirname(__file__)
    resultsdir  = './temp_{0}/'.format(
        os.path.splitext(input_pdb_file)[0].replace('/', '_')
    )
    output_score_file_prefix = 'score'
    external_scores_dict = external_scoring.generate_enhanced_score_summary_table(
        input_pose,
        scriptsdir,
        resultsdir,
        output_score_file_prefix
    )
    external_scores = pandas.DataFrame.from_dict(external_scores_dict)
    external_scores.set_index('Unnamed: 0', inplace=True)

    # Remove the temporary results directory
    subprocess.check_call(['rm', '-r', resultsdir])

    # Merge scores
    assert all(rs_scores.index.values == external_scores.index.values)
    rs_metrics = list(rs_scores.columns.values)
    external_metrics = list(external_scores.columns.values)
    unique_external_metrics = [m for m in external_metrics if m not in rs_metrics]
    merged_scores = pandas.concat([rs_scores, external_scores[unique_external_metrics]], axis=1)

    return merged_scores
