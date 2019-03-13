#! /usr/bin/env python
import numpy as np
import sqlite3
import os

from TaxLink import taxlink

from applicake.app import BasicApp
from applicake.coreutils.arguments import Argument
from applicake.coreutils.keys import Keys, KeyHelp

class ms2(BasicApp):
    """
    ## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	## ms2 analysis
	## Search computational XLs in MSMS data (mgf format) to find and analyze real patterns.

	## <<Inputs>>
	## file containing all xls in kojak format
    ## mgf file
    ## delta value which is 0.05 or 0.01
    ## intensity cut-off 0 to any number

	## <<Output>>
	## A SQlite file which the png figures of top-spectra
	## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    """

    def add_args(self):
        return [
            Argument('xlink_f', 'xlink_f'),
            Argument('mgf_f', 'mgf_f'),
            Argument('delta', 'delta'),
            Argument('intensity', 'intensity'),
            Argument('to_dropbox', 'to_dropbox'),
            Argument('dir_to_dropbox', 'dir_to_dropbox'),
            Argument(Keys.WORKDIR, KeyHelp.WORKDIR),
        ]

    def run(self, log, info):
        wd = info['WORKDIR']

        if 'dir_to_dropbox' not in info:
            info['dir_to_dropbox'] = []
        if 'to_dropbox' not in info:
            info['to_dropbox'] = []

        list_of_XL_files = info['xlink_f']
        score = 0
        score_list = []

        for xl_file in list_of_XL_files:
            score = taxlink(xl_file,info['mgf_f'],float(info['delta']),float(info['intensity']))
            score_list.append(score)
	    print "MS2 analysis is finished for ", xl_file, " with score ", score

        top_model_xl = list_of_XL_files[score_list.index(max(score_list))]
        top_model = top_model_xl[:-7]+"0.pdb"
        print "Final selected model is: ", top_model
        info['top_model'] = top_model

        command1 = "cp %s %s/best_model.pdb" %(top_model, wd)
        os.system(command1)

        command2 = "cp %s %s/best_model_xl.txt" %(top_model_xl, wd)
        os.system(command2)

        info['to_dropbox'].append('best_model.pdb')
        info['to_dropbox'].append('best_model_xl.txt')
        info['to_dropbox'].append('MS2_results.db')
        info['dir_to_dropbox'].append('Top_MSMS_spectra')
        info['MS2_sql_file'] = "MS2_results.db"

        return info


if __name__ == "__main__":
    ms2.main()
