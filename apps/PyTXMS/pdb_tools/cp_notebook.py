#!/usr/bin/env python

import glob
import os
#import ntpath

from applicake.app import BasicApp
from applicake.coreutils.arguments import Argument
from applicake.coreutils.keys import Keys, KeyHelp


class cp_notebook(BasicApp):
    """
    ##
    ##
    """


    def add_args(self):
        return [
            Argument('notebook_dataset', 'notebook_dataset'),
            Argument('html_dataset', 'html_dataset'),
            Argument('to_dropbox', 'to_dropbox'),
            Argument('DATASET_DIR', 'DATASET_DIR'),
            Argument(Keys.WORKDIR, KeyHelp.WORKDIR),
        ]

    def run(self, log, info):
        wd = info['WORKDIR']


        if 'to_dropbox' not in info:
            info['to_dropbox'] = []


        nb_files = glob.glob("%s/%s/*.ipynb"%(info['DATASET_DIR'],info['notebook_dataset']))
        print nb_files
        notebook = nb_files[0]


        html_files = glob.glob("%s/%s/*.html"%(info['DATASET_DIR'],info['html_dataset']))
        print html_files
        html_page = html_files[0]


        command1 = "cp %s %s/raw.000_report.ipynb" %(notebook, wd)
        os.system(command1)

        command2 = "cp %s %s/index.html" %(html_page, wd)
        os.system(command2)

        info['raw_notebooks'] = "raw.000_report.ipynb"
        info['to_dropbox'].append('index.html')

        return info


if __name__ == "__main__":
    cp_notebook.main()

