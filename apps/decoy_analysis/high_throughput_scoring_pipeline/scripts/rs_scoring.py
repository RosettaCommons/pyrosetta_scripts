#!/usr/bin/env python3

"""
This script uses `unittest` to test that input RosettaScripts XML files
from `./scripts/xmls/` can be parsed and used to score a set of test designs
from `./data/Rocklin_rd4_PDB_files/` without giving rise to any errors.

This script also contains varions functions used to score poses both with
XML RosettaScripts and with custom Python scripts.

This script is intended to be run in the base directory for this analysis,
i.e., one directory closer to the root than the directory this script is in,
using the command:
    python test_scoring.py

The top line of this script makes it into an executable that runs in the
the environment:
    `/home/haddox/.conda/envs/1711_term_project/bin/python`
which is the same environment as the kernel of the Jupyter notebook
`test_score_protocols.ipynb`, so that I can use that environment's version
of `pyrosetta`.
"""

# Imports
import unittest
import logging
import os
from os import path
import sys
import glob
import random
import subprocess
import pandas
import pyrosetta
import pyrosetta.distributed.io as io
import pyrosetta.distributed.tasks.score as score
import pyrosetta.distributed.tasks.rosetta_scripts as rosetta_scripts

logger = logging.getLogger(__name__)

#TODO refactor to add informative top-level module name for `scripts` directory
script_dir = path.dirname(__file__)
sys.path.append(script_dir)
from xmls import analysis_scripts
from test_data.rocklin_rd4_examples import test_pdb_files
assert sys.path.pop() == script_dir

# Functions for scoring designs and tabulating the results
def generate_pose_score_summary(input_pose, input_xmls):
    """Create score summary dict for given pose for specific profiles.
    Args:
        input_pose - Packed pose.
        input_xmls - {name : script_content} dict of protocol scripts, which report a
            final score under 'name'.

    Returns:
        {name : metric} dict of metric values.
    """

    score_result = {}
    score_result.update(score.ScorePoseTask()(input_pose).scores)

    # cycle through the scoring xmls, each of which specifies
    # a single filter. Add the results of each filter to the
    # scores for each design.
    for feature_name, script_content in input_xmls.items():
        try:
            xml_score = rosetta_scripts.SingleoutputRosettaScriptsTask(
                    script_content)(input_pose).scores
            score_result[feature_name] = xml_score[feature_name]
        except RuntimeError:
            logger.exception("Error processing feature: %s" % feature_name)
            score_result[feature_name] = None

    return score_result

def generate_score_summary_table(input_poses, input_xmls):
    score_rows = {}
    for name, input_pose in input_poses.items():
        score_result = generate_pose_score_summary(input_pose, input_xmls)

        assert name not in score_rows
        score_rows[name] = score_result

    return pandas.DataFrame.from_dict(score_rows, orient="index")

# Define a class for testing
class TestRosettaScriptsScoring(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.input_poses = {
            f : io.pose_from_file(f)
            for f in test_pdb_files[:1]
        }

        cls.input_xmls = analysis_scripts

    def test_xml_parsing(self):
        """This function tests that Rosetta is able to parse each of
        the input XML files"""

        self.assertTrue( len(self.input_xmls) > 1)

        for test_xml in self.input_xmls:

            # log the XML being tested
            logging.info("Testing the XML: {0}".format(test_xml))

            # setup a rosetta_scripts task to trigger xml->tag load
            # and init of rosetta scripts parser
            dummy_task = rosetta_scripts.SingleoutputRosettaScriptsTask(self.input_xmls[test_xml])
            dummy_task.setup()

            # pass tag into rosetta scripts parser to validate
            # that xml was successfully parsed
            protocol = dummy_task.parser.parse_protocol_tag(
                dummy_task.tag, dummy_task.default_options)

    def test_scoring(self):
        """This function tests whether we can use the XML files to
        score each test design"""

        self.assertTrue( len(self.input_poses) >= 1)

        # use Rosetta to score each of the poses using standard
        # metrics
        # dataframe with metrics as columns and pdbs as rows
        score_summary = \
            generate_score_summary_table(self.input_poses, self.input_xmls)

if __name__ == "__main__":
    logging.basicConfig(filename='test_scoring.log',level=logging.DEBUG)
    unittest.main()
