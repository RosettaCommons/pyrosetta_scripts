import sys
import argparse
from pyrosetta import*

init(extra_options="-ignore_unrecognized_res")

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

print(args.filename)