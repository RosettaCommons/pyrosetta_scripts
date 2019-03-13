#!/usr/bin/env python
import os
import glob
from applicake.app import WrappedApp
from applicake.apputils.validation import check_exitcode, check_xml, check_stdout
from applicake.coreutils.arguments import Argument
from applicake.coreutils.keys import Keys, KeyHelp

class dinosaur(WrappedApp):
    """
    Wrapper for dinosaur tool
    """

    def add_args(self):
        return [
            Argument(Keys.WORKDIR, KeyHelp.WORKDIR),
            Argument('ms1_datasets', 'ms1_datasets'),
            Argument('java', 'java'),
            Argument('dinosaur_jar', 'dinosaur_jar'),
            Argument('DATASET_DIR', 'DATASET_DIR'),
            Argument('datasets', 'datasets'),
        ]

    def prepare_run(self, log, info):
        wd = info[Keys.WORKDIR]

        ##
        # if 'ms1_datasets' in info:
	#     for mzml in info['ms1_datasets']:
        #         mzmlfiles = glob.glob("%s/%s/*.mzML*"%(info['DATASET_DIR'],mzml))
        #         print("Found %s in %s"%(mzmlfiles,mzml))
        #     dinosaur_file = mzmlfiles[0]
        mzml_files = glob.glob("%s/%s/*.mzML*"%(info['DATASET_DIR'],info['ms1_datasets']))
        dinosaur_file = mzml_files[0]
        ##
        
        command = info["java"] + " -jar " + info['dinosaur_jar'] + " --verbose --profiling --concurrency=4 --outDir="+ wd+ " "+ dinosaur_file
        return info, command

    def validate_run(self, log, info, exit_code, stdout):
        check_stdout(log, stdout)
        check_exitcode(log, exit_code)
        feature_files = glob.glob("%s/*.features.tsv"%(info[Keys.WORKDIR]));
        os.system("mv qc %s"%(os.path.basename(info['ms1_datasets']).split('.')[0]))
        info['dinosaur_feature_file'] = feature_files[0]
        return info


if __name__ == "__main__":
    dinosaur.main()
