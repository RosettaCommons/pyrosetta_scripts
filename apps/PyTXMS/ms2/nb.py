#!/usr/bin/env python3.5
import os

from applicake.app import BasicApp
from applicake.coreutils.arguments import Argument
from applicake.coreutils.keys import Keys, KeyHelp

import pyddb.nbunit
import pandas
import glob
import jinja2


class swath_nb(BasicApp):

    def add_args(self):
        return [
            Argument('sqlitedb', 'sqlitedb'),
            Argument('DATASET_DIR', 'DATASET_DIR'),
            Argument('CACHE_DATASET_DIR', 'CACHE_DATASET_DIR'),
            Argument('to_dropbox', 'to_dropbox'),
        ]

    def run(self, log, info):
        """ initialize reports """
        nb = pyddb.nbunit.nbunit(info['sqlitedb'])
        if 'to_dropbox' not in info:
            info['to_dropbox'] = []
        nb.add_template('000_report','/vdb1/lib/applicake/appliapps/swath/cell.python.annotations.j2',info)
        info['to_dropbox'].append("report_grouped_by_transition_group_id.pdf")

        return info

if __name__ == "__main__":
    swath_nb.main()
