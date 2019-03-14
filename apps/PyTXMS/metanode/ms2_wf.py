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


class metanode_ms2_wf(WrappedApp):
    """ metanode """

    def add_args(self):
        return [
            Argument('sequence1', 'sequence1'),
            Argument('sequence2', 'sequence2'),
            Argument('datasets', 'datasets'), ## ms2 mzml file
            Argument('delta', 'delta'),
            Argument('intensity', 'intensity'),
            Argument('msconvert_executable', 'msconvert_executable'),
            Argument(Keys.WORKDIR, KeyHelp.WORKDIR),
        ]

    def run(self, log, info):
        wd = os.getcwd()

        ## (1) ## running the seq2xl_app to make a full XL set in kojak format

        self._dr('openbis.dss_simple',info['INPUT'],'vars.ini',wd,[])
        self._dr('seq2xl.seq2xl_app','vars.ini','vars.ini',wd,[])

        #self._dr('seq2xl.seq2xl_app',info['INPUT'],'vars.ini',wd,[])

        ## (2) ## running ms2
        self._dr('ms2.msconvert','vars.ini','vars.ini',wd,[])
        self._dr('ms2.ms2','vars.ini','vars.ini',wd,[])

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
    metanode_ms2_wf.main()
