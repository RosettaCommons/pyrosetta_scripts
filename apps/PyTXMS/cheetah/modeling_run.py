#! /usr/bin/env python
import numpy as np
import sqlite3
import glob
from modeling_core import modeling_core

from applicake.app import BasicApp
from applicake.coreutils.arguments import Argument
from applicake.coreutils.keys import Keys, KeyHelp

class modeling_run(BasicApp):
    """
    ## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	## modeling_run.py
	## PyRosetta modeling by using rosettadock and XL file
	## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    """

    def add_args(self):
        return [
            Argument('merged_pdb', 'merged_pdb'),
            Argument('partners', 'partners'),
            Argument('top_XL_file', 'top_XL_file'),
            Argument('cut_off', 'cut_off'),
            Argument('num_of_models', 'num_of_models'),
            Argument('num_of_top_filters', 'num_of_top_filters'),
            Argument('DATASET_DIR', 'DATASET_DIR'),
            Argument('dir_to_dropbox', 'dir_to_dropbox'),
            Argument(Keys.WORKDIR, KeyHelp.WORKDIR),
        ]

    def run(self, log, info):
        wd = info['WORKDIR']

        if 'dir_to_dropbox' not in info:
            info['dir_to_dropbox'] = []
	##
        inp_pdb_file = info['merged_pdb']
	##

        list_of_XL_files =  modeling_core(inp_pdb_file, info['partners'], info['top_XL_file'], int(info['cut_off']), int(info['num_of_models']), int(info['num_of_top_filters']))

	info['xlink_f'] = list_of_XL_files
        info['dir_to_dropbox'].append('LR_models')
        return info


if __name__ == "__main__":
    modeling_run.main()

