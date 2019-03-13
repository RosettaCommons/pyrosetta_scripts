#!/usr/bin/env python3
import os
import re
import glob
import importlib
import inspect
import sys

from applicake.app import WrappedApp
from applicake.apputils import validation
from applicake.coreutils.arguments import Argument
from applicake.coreutils.keys import Keys, KeyHelp
from applicake.coreutils.info import get_handler


class metanode_cheetah_wf(WrappedApp):
    """ metanode """

    def add_args(self):
        return [
            Argument('ms1_datasets', 'ms1_datasets'), ## ms1 mzml file
            Argument('ms2_datasets', 'ms2_datasets'), ## ms2 mzml file
            Argument('java', 'java'),
            Argument('dinosaur_jar', 'dinosaur_jar'),
            Argument('taxlink_jar', 'taxlink_jar'),
            Argument('train_datasets', 'train_datasets'),
            Argument('ensemble', 'ensemble'),
            Argument('pdb_datasets', 'pdb_datasets'),
            Argument('partners', 'partners'),
            Argument('cut_off', 'cut_off'),
            Argument('num_of_models', 'num_of_models'),
            Argument('num_of_top_filters', 'num_of_top_filters'),
            Argument('kfold', 'kfold'),
            Argument('delta', 'delta'),
            Argument('intensity', 'intensity'),
            Argument('msconvert_executable', 'msconvert_executable'),
            #Argument('to_dropbox', 'to_dropbox'),
            #Argument('dir_to_dropbox', 'dir_to_dropbox'),
            Argument(Keys.WORKDIR, KeyHelp.WORKDIR),
        ]

    def run(self, log, info):
        wd = os.getcwd()

        ## running seq2xl to make a list of potential XLs from input seqs
        self._dr('openbis.dss_simple',info['INPUT'],'vars.ini',wd,[])
        self._dr('pdb_tools.single_chain','vars.ini','vars.ini',wd,[])
        self._dr('pdb_tools.cp_notebook','vars.ini','vars.ini',wd,[])
        self._dr('seq2xl.seq2xl_app','vars.ini','vars.ini',wd,[])

        ## running ML based hrms1 analysis to find potential binding interfaces
        self._dr('dinosaur.dinosaur','vars.ini','vars.ini',wd,[])
        self._dr('taxlink.taxlink','vars.ini','vars.ini',wd,[])
        self._dr('hrms1.hrms1','vars.ini','vars.ini',wd,[])

        ## running cheetah modeling on the input structure with potential BIs
        self._dr('cheetah.modeling_run','vars.ini','vars.ini',wd,[])

        ## validating the top model's XLs with ms2 analysis
        self._dr('ms2.msconvert','vars.ini','vars.ini',wd,[])
        self._dr('ms2.ms2','vars.ini','vars.ini',wd,[])

        ## making a pymol session with all the information
        self._dr('cheetah.pymol_run','vars.ini','vars.ini',wd,[])

        ## providing the information to the dropbox and making a report
        # self._dr('cheetah.make_report','vars.ini','vars.ini',wd,[])

        return info

    def validate_run(self, log, info, exit_code, out):
        validation.check_exitcode(log, exit_code)
        return info

    def _dr(self, appliapp, inp, outp, workdir,extra_argv):
        """ dynrunner """

        cls = None
        sys.argv = [sys.argv[0],sys.argv[1],'--INPUT',inp,'--OUTPUT',outp,'--WORKDIR',workdir]
        sys.argv.extend(extra_argv)
        try:
            #load app module, from http://stackoverflow.com/a/10773699
            module = importlib.import_module(appliapp)
            #find main class http://stackoverflow.com/q/1796180/
            for name, obj in inspect.getmembers(module):
                if inspect.isclass(obj) and appliapp in obj.__module__:
                    cls = obj
        except Exception as e:
            raise Exception('Could not find/load app [%s]: %s' % (appliapp, str(e)))
        try:
            cls.main()
        except Exception as e:
            raise Exception('Could not run app [%s]: %s' % (appliapp, str(e)))

if __name__ == "__main__":
    metanode_cheetah_wf.main()

