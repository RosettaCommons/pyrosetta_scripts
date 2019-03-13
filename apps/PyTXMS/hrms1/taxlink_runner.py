#!/usr/bin/env python

"""
Created on Wed Feb  8 15:26:07 2017

@author: hamedkhakzad
"""

import subprocess

def taxlink_runner(tsv_file_name, XL_set):
    subprocess.check_call('java -jar ~/Applications/taxlink/Taxlink-0.8.1.jar ' + tsv_file_name\
    + ' ' + XL_set, shell=True)
