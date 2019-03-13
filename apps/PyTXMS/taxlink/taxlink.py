#!/usr/bin/env python
import os
import glob
from applicake.app import BasicApp
from applicake.apputils.validation import check_exitcode, check_xml, check_stdout
from applicake.coreutils.arguments import Argument
from applicake.coreutils.keys import Keys, KeyHelp

class taxlink(BasicApp):
    """
    Wrapper for taxlink tool
    """

    def add_args(self):
        return [
            Argument(Keys.WORKDIR, KeyHelp.WORKDIR),
            Argument('java', 'java'),
            Argument('taxlink_jar', 'taxlink_jar'),
            Argument('kojak_xl_file', 'kojak_xl_file'),
            Argument('dinosaur_feature_file', 'dinosaur_feature_file'),
        ]

    def run(self, log, info):
        wd = info[Keys.WORKDIR]

#        if not isinstance(info['kojak_xl_file'], list):
#            info['kojak_xl_file'] = [info['kojak_xl_file']]
#        if not isinstance(info['dinosaur_feature_file'], list):
#            info['dinosaur_feature_file'] = [info['dinosaur_feature_file']]
        
        command = info["java"] + " -jar " + info['taxlink_jar'] + " " + info["dinosaur_feature_file"] + " " + info["kojak_xl_file"]
        print("Command is: %s"%(command))
        os.system(command)
        
        taxlink_files = glob.glob("*.isopairs.csv");
        print("Found: %s"%(",".join(taxlink_files)))
        info['taxlink_data'] = taxlink_files[0]
        return info



if __name__ == "__main__":
    taxlink.main()
