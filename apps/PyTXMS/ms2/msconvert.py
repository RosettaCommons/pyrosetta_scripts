#!/usr/bin/env python
import os
import glob
from applicake.app import WrappedApp
from applicake.apputils.validation import check_exitcode, check_xml, check_stdout
from applicake.coreutils.arguments import Argument
from applicake.coreutils.keys import Keys, KeyHelp

class msconvert(WrappedApp):
    """
    Wrapper for msconvert tool
    """

    def add_args(self):
        return [
            Argument(Keys.WORKDIR, KeyHelp.WORKDIR),
            Argument('ms2_datasets', 'ms2_datasets'),
            Argument('DATASET_DIR', 'DATASET_DIR'),
            Argument('datasets', 'datasets'),
            Argument('msconvert_executable', 'msconvert_executable'),
        ]

    def prepare_run(self, log, info):
        wd = info[Keys.WORKDIR]

        ##
        mzml_files = glob.glob("%s/%s"%(info['DATASET_DIR'],info['ms2_datasets']))
        msconvert_file = mzml_files[0]
        ##

        command = info['msconvert_executable'] + " " +msconvert_file + " --mgf -o "+ wd
        print command
        return info, command

    def validate_run(self, log, info, exit_code, stdout):
        check_stdout(log, stdout)
        check_exitcode(log, exit_code)
        mgf_files = glob.glob("*.mgf");
        print("Found: %s"%(",".join(mgf_files)))
        info['mgf_f'] = mgf_files[0]
        return info


if __name__ == "__main__":
    msconvert.main()
