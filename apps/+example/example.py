#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   example.py
## @brief  Example PyRosetta app
## @author Sergey Lyskov

from __future__ import print_function

import sys, argparse

import pyrosetta

def app_function(s):
    ''' dummy app function
    '''
    return '__q' + s + 'w__'


def main(args):
    ''' Example PyRosetta app script '''

    parser = argparse.ArgumentParser()
    parser.add_argument('-j', '--jobs', default=1, const=0, nargs="?", type=int, help="Number of processors to use on when building, use '-j' with no arguments to launch job-per-core. (default: 1) ")
    parser.add_argument("--multi-choice", default='Release', choices=['Release', 'Debug', 'MinSizeRel', 'RelWithDebInfo'], help="Specify one of the allowed values")
    parser.add_argument('-s', '--skip-generation-phase', action="store_true", help="Example of true/false option")

    global Options
    Options = parser.parse_args()

    print( 'PyRosetta example app, running with options: {Options}...'.format(**vars()) )
    print('PyRosetta example app... Done!')


if __name__ == "__main__":
    main(sys.argv)
