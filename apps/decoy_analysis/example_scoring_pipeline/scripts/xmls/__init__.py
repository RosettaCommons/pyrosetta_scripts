import os
from os import path

import glob

__all__ = ["analysis_scripts"]

script_files = glob.glob(path.join(path.dirname(__file__), "*.xml"))

def _read(f):
    with open(f) as inf:
        return inf.read()

analysis_scripts = {
    path.splitext(path.basename(x))[0] : _read(x)
    for x in script_files 
}
