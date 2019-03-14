#! /usr/bin/env python
import numpy as np
import os
import glob

from rosetta import *
from pyrosetta import *
from rosetta.protocols.rigid import *

from applicake.app import BasicApp
from applicake.coreutils.arguments import Argument
from applicake.coreutils.keys import Keys, KeyHelp

class single_chain(BasicApp):
    """
    ## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	## single_chain is a applicake node to convert pdb files to a single chain pdb.
    ## only works with two pdbs!

	## <<Input>>
	## pdb_files datasets

	## <<Output>>
	## single chain pdb files, sequences
	## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    """

    def add_args(self):
        return [
            Argument('pdb_datasets', 'pdb_datasets'),
            Argument('to_dropbox', 'to_dropbox'),
            Argument('DATASET_DIR', 'DATASET_DIR'),
            Argument(Keys.WORKDIR, KeyHelp.WORKDIR),
        ]

    def run(self, log, info):
        wd = info['WORKDIR']
        info["single_chain_pdbs"] = []
        info["single_chain_seqs"] = []

        if not 'to_dropbox' in info:
            info['to_dropbox'] = [] 

        init()

        pdb_files = []
        if 'pdb_datasets' in info:
            for pds in info['pdb_datasets']:
                files = glob.glob("%s"%(pds)) ## without dss_simple
                print("Found %s in %s"%(files,pds))
                pdb_files.append(files[0])

            pdb1_file = pdb_files[0]
            pdb2_file = pdb_files[1]

	    print ("PDB files: ", pdb_files)
        else:
            print ("ERROR! There is no PDB_datasets file in the input.")

        out_chains = ['A', 'B']
        for num, pdbs in enumerate(pdb_files):
            ## make a single chain, cleaned pdb
            pose = Pose()
            pose_from_file(pose, pdbs)
            nres = pose.total_residue()
            sequence = pose.sequence()
            pdb_info = pose.pdb_info()
            chains = [pdb_info.chain(i) for i in range(1, nres + 1)]
            unique_chains = []
            for c in chains:
                if c not in unique_chains:
                    unique_chains.append(str(c))
	    
            chains_all = ','.join(unique_chains)

            ## making the sequence file
            seq_file_name = "sequence_"+out_chains[num]+".txt"
            seq_file = open(seq_file_name, "w")
            seq_file.write("%s" %sequence)
            seq_file.close()
            info["single_chain_seqs"].append(seq_file_name)
	    info['to_dropbox'].append(seq_file_name)

            ## making the pdb file
            out_pdb_file_name = "pdb_"+out_chains[num]+".pdb"
            command = "pdb_reres -1 %s | pdb_selchain -%s | pdb_chain -%s | pdb_delhetatm > %s" \
                %(pdbs, chains_all, out_chains[num], out_pdb_file_name)            
            print(command)
            os.system(command)
            info["single_chain_pdbs"].append(out_pdb_file_name)
	    info['to_dropbox'].append(out_pdb_file_name)
        
        command2 = "pdb_merge pdb_A.pdb pdb_B.pdb | pdb_reres -1 > pdb_C.pdb"
        os.system(command2)

	info['merged_pdb'] = "pdb_C.pdb"
	info['partners'] = "A_B"
	info['to_dropbox'].append("pdb_C.pdb")
        return info

if __name__ == "__main__":
    single_chain.main()

