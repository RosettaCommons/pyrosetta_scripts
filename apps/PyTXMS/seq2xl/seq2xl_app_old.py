#! /usr/bin/env python
import numpy as np
import sqlite3
from kojak_generator import kojak_generator
import glob

from applicake.app import BasicApp
from applicake.coreutils.arguments import Argument
from applicake.coreutils.keys import Keys, KeyHelp

class seq2xl(BasicApp):
    """
    ## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	## seq2xl is a python app to generate all inte XLs
	## between two input aa sequences.

	## <<Input>>
	## Two aa sequences from command line

	## <<Output>>
	## A txt file contains all inter XLs
	## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    """

    def add_args(self):
        return [
            Argument('sequence1', 'sequence1'),
            Argument('sequence2', 'sequence2'),
            Argument('seq_datasets', 'seq_datasets'),
            #Argument('mzml_datasets', 'mzml_datasets'),
            Argument('DATASET_DIR', 'DATASET_DIR'),
            Argument('to_dropbox', 'to_dropbox'),
            Argument(Keys.WORKDIR, KeyHelp.WORKDIR),
        ]

    def run(self, log, info):
        wd = info['WORKDIR']
        seqfiles = []
        if not 'to_dropbox' in info:
            info['to_dropbox'] = []
        if not 'java' in info:
            info['java'] = 'java'
        if not 'dinosaur_jar' in info:
            info['dinosaur_jar'] = '/dinosaur/target/Dinosaur-1.1.3.free.jar'
        #mzml_files = glob.glob("%s/%s/*.mzML.gz"%(info['DATASET_DIR'],info['mzml_datasets']))
        #info['MZML'] = mzml_files[0]
        if 'seq_datasets' in info:
            for sds in info['seq_datasets']:
                files = glob.glob("%s/%s/*.txt"%(info['DATASET_DIR'],sds))
                print("Found %s in %s"%(files,sds))
                seqfiles.append(files[0])
            seq1_file = seqfiles[0]
            seq2_file = seqfiles[1]
            print(seq1_file,seq2_file)
        else:
            seq1_file = info['sequence1']
            seq2_file = info['sequence2']

        PEP_LEN = 4
	seqlist1 = []
        seqlist2 = []
	kojak_list = []

        ## Importing sequence files file
        with open(seq1_file,'r') as sequence1:
            seq1_list = sequence1.read()
        with open(seq2_file,'r') as sequence2:
            seq2_list = sequence2.read()    

        ALLXL_file = open('ALL_XL.txt', 'w')
        ALLXL_file.write("sequence\n")
        
        ## finding and storing all inter_XLs
        for i in range(len(seq1_list)):
            
            if seq1_list[i] == ('K'):
                pep1,K_pos1 = kojak_generator(seq1_list,i)

                if len(pep1) >= PEP_LEN: 
                
                    for j in range(len(seq2_list)):
                        
                        if seq2_list[j] == ('K'):
                            pep2,K_pos2 = kojak_generator(seq2_list,j)

                            if len(pep2) >= PEP_LEN:

                                if pep1 != pep2:
                                
                                    XL_kojak_format = "-."+pep1+"("+str(K_pos1)+")--"+pep2+"("+str(K_pos2)+").-"
                                    XL_kojak_format_rev = "-."+pep2+"("+str(K_pos2)+")--"+pep1+"("+str(K_pos1)+").-"

                                    if (XL_kojak_format not in kojak_list) and (XL_kojak_format_rev not in kojak_list):
                                        kojak_list.append(XL_kojak_format)
                                        ALLXL_file.write("%s\n" % XL_kojak_format)

        info['kojak_xl_file'] = 'ALL_XL.txt'
        info['xlink_f'] = 'ALL_XL.txt'
        ## send to dropbox
        info['to_dropbox'].append(seq1_file)
        info['to_dropbox'].append(seq2_file)
        return info

if __name__ == "__main__":
    seq2xl.main()
