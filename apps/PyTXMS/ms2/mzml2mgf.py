#!/usr/bin/env python
import os
import re
import sys
from applicake.app import BasicApp
#import jinja2
import glob
from applicake.apputils.validation import check_exitcode, check_xml, check_stdout
from applicake.coreutils.arguments import Argument
from applicake.coreutils.keys import Keys, KeyHelp

class mzml2mgf(BasicApp):
    """
    Wrapper for pwiz msconvert
    """

    def add_args(self):
        return [
            Argument(Keys.WORKDIR, KeyHelp.WORKDIR),
            Argument('DATASET_DIR', 'DATASET_DIR'),
            Argument('CACHE_DATASET_DIR', 'CACHE_DATASET_DIR'),
            Argument('datasets', 'datasets'),
            Argument('mzml_datasets', 'mzml_datasets'),
            Argument('cache_in_dss', 'cache_in_dss'),
            Argument('msconvert_executable', 'msconvert_executable'),
        ]

    def run(self, log, info):
        wd = info[Keys.WORKDIR]
        mgfs = []
        mzmls = []
        dss = []
        if 'datasets' in info:
            dss = info['datasets']
        if 'mzml_datasets' in info: # overwrite datasets and this is wanted
            dss = info['mzml_datasets']
        if not isinstance(dss,list):
            dss = [dss]
        print("Number of datasets found: %s"%(len(dss)))
        for ds in dss:
            print("Attempt to convert: %s"%(ds))
            if ds == "20161101134412891-150291":
                continue
            mzxml_files = glob.glob("%s/%s/*.mzML.gz"%(info['DATASET_DIR'],ds));
            #mgf = "%s/%s/%s"%(info['CACHE_DATASET_DIR'],ds,os.path.basename(mzxml_files[0]).replace('mzML.gz','pwiz.mgf'))
            mgf = "%s/%s/%s"%(info['CACHE_DATASET_DIR'],ds,os.path.basename(mzxml_files[0]).replace('.mzML.gz','.mgf'))
            mzmls.append(mzxml_files[0])
            if (os.path.isfile(mgf)):
                print("have mgf %s"%(mgf))
                find_scan = re.compile("TITLE=[^\.]+\.(\d+)\.")
                newfile = "%s/%s.mgf"%(wd,os.path.basename(mgf).split('.')[0])
                out = open(newfile,'w')
                mgfs.append(newfile)
                with open(mgf) as f:
                    for line in f:
                        out.write(line)
                        match = find_scan.match(line)
                        if match:
                            out.write("SCANS=%s\n"%(match.group(1)))
            else:
                #raise IOError("Not now...")
                output_dir = wd
                if 'cache_in_dss' in info and info['cache_in_dss']:
                    output_dir = "%s/%s"%(info['CACHE_DATASET_DIR'],ds)
                expected_mgf_file = os.path.join(output_dir,mzxml_files[0].split("/")[-1].replace(".mzML.gz",".mgf"))
                print("Expected: %s"%(expected_mgf_file))
                print("BEFORE")
                if os.path.isfile(expected_mgf_file):
                    print("HAVE %s, not converting"%(expected_mgf_file))
                    os.system("ln -s %s"%(expected_mgf_file))
                else:
                    os.system("%s %s -o %s --%s --filter 'msLevel 2'" % (info['msconvert_executable'], mzxml_files[0],output_dir,'mgf'))
                print("AFTER")
                mgfs.append(mgf)
        info['MZML'] = mzmls
        info['mgf'] = mgfs
        return info

    def validate_run(self, log, info, exit_code, stdout):
        check_stdout(log, stdout)
        check_exitcode(log, exit_code)
        return info


if __name__ == "__main__":
    mzml2mgf.main()
