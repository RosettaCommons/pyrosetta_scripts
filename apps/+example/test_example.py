#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   test_example.py
## @brief  test file example.py app
## @author Sergey Lyskov

from __future__ import print_function

import example

import pytest


def test_app_function():
    ''' defining test function as 'test_' + test-name, such function will be automatically executed by test framework
    '''
    assert example.app_function('A') == '__qAw__'


# def test_app_function2():
#     ''' This test will always fail
#     '''
#     assert False
