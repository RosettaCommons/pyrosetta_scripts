import os
from os import path

import glob

__all__ = ["test_pdb_files"]

test_pdb_files = glob.glob(path.join(path.dirname(__file__), "*.pdb"))
