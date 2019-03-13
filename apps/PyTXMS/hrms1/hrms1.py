#!/usr/bin/env python

#from taxlink_runner import taxlink_runner
from TX_Learning_Train import TX_Learning_Train
from TX_Learning_Test import TX_Learning_Test

from applicake.app import BasicApp
from applicake.coreutils.arguments import Argument
from applicake.coreutils.keys import Keys, KeyHelp

import glob

class hrms1(BasicApp):
    """
    ## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ## hrms1 Machine Learning based MS1 analysis to find XL 
    ## patterns in mass spec MS1 data.

    ## <<Input>>
    ## Traning set data, All_XL file, and mzML data in isopairs format (output of taxlink).

    ## <<Output>>
    ## list of top XLs in a txt and sqlite3 format.
    ## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    """

    def add_args(self):
        return [
            Argument('train_datasets', 'train_datasets'),
            Argument('taxlink_data', 'taxlink_data'),
            Argument('ensemble', 'ensemble'),
            Argument('kfold', 'kfold'),
            Argument('DATASET_DIR', 'DATASET_DIR'),
            Argument('to_dropbox', 'to_dropbox'),
            Argument(Keys.WORKDIR, KeyHelp.WORKDIR),
        ]

    def run(self, log, info):
        wd = info['WORKDIR']

        if 'to_dropbox' not in info:
            info['to_dropbox'] = []

        print "Downloading the training file ..."
	##
        training_files = glob.glob("%s/%s/*.csv"%(info['DATASET_DIR'],info['train_datasets']))
	print training_files
        training_set = training_files[0]
        ##

        print "training phase:"
        model = TX_Learning_Train(training_set, int(info['ensemble']), int(info['kfold']))
        
        print "training is done, wait for the results ..."
        TX_Learning_Test(model, info['taxlink_data'], 1)

        info['hrms1_top_xl_file'] = 'top_XL'+info['taxlink_data']+'.txt'
        info['top_XL_file'] = 'top_XL'+info['taxlink_data']+'.txt'
        info['to_dropbox'].append('hrMS1_result1.db')
        return info


if __name__ == "__main__":
    hrms1.main()

