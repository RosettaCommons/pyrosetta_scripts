#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   test_rna_denovo.py
## @brief  test file rna_denovo.py app
## @author Andy Watkins

from __future__ import print_function

import rna_denovo

import pytest
from pyrosetta import *
from pyrosetta.rosetta import *

import os

def test_app_function():
    init("-constant_seed -rna:denovo:cycles 5")

    r = os.path.dirname(os.path.abspath(__file__)) + '/'

    rna_denovo.denovo([r+"gcaa_tetraloop_HELIX1.pdb"],
                      r+"gcaa_tetraloop_NATIVE_1zih_RNA.pdb",
                      r+"gcaa_tetraloop_NATIVE_1zih_RNA.pdb",
                      r+"gcaa_tetraloop.fasta",
                      1,
                      r+"test_denovo.out")
