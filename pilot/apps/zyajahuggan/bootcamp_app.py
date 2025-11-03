import sys
import argparse
import pyrosetta
pyrosetta.init()
from pyrosetta import *

init(extra_options="-ignore_unrecognized_res")

parser = argparse.ArgumentParser()

args = parser.parse_args()
