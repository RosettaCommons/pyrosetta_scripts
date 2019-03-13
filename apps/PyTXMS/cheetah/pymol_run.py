#!/usr/bin/env python

import sqlite3
import numpy as np
import __main__
__main__.pymol_argv = [ 'pymol', '-qc' ]

import pymol
from pymol import cmd

from rosetta import *
from pyrosetta import *
from rosetta.protocols.rigid import *

from string import digits
from math import sqrt, pow, exp
from rosetta_eu_dist import rosetta_eu_dist

from applicake.app import BasicApp
from applicake.coreutils.arguments import Argument
from applicake.coreutils.keys import Keys, KeyHelp


class pymol_run(BasicApp):
    """
    ##
    making a pymol session of the top model and top XLs mapped on the structure.
    the output file contains those XLs that are within the predefined cut-off threshold.
    ##
    """
    # def rosetta_eu_dist(pose, a, b):
    #     dist = sqrt(pow(( pose.residue( b ).atom( "CB" ).xyz().x - pose.residue( a ).atom( "CB" ).xyz().x ) , 2) \
    #         + pow(( pose.residue( b ).atom( "CB" ).xyz().y - pose.residue( a ).atom( "CB" ).xyz().y ) , 2) \
    #         + pow(( pose.residue( b ).atom( "CB" ).xyz().z - pose.residue( a ).atom( "CB" ).xyz().z ) , 2))
    #     return dist


    def add_args(self):
        return [
            Argument('top_model', 'top_model'),
            Argument('cut_off', 'cut_off'),
            Argument('MS2_sql_file', 'MS2_sql_file'),
            Argument('to_dropbox', 'to_dropbox'),
            Argument(Keys.WORKDIR, KeyHelp.WORKDIR),
        ]

    def run(self, log, info):
        init()
        wd = info['WORKDIR']

        pdb_filename = info['top_model']
        top_XL_file = pdb_filename[:-5]+"XLs.txt"
        found_xl_file = open('mapped_xl.txt', 'w')
        cut_off = info['cut_off']

        ## 0. Parsing MS2 SQLite file to select only top XLs from MS/MS analysis
        ms2_table = info['MS2_sql_file']
        conn0 = sqlite3.connect(ms2_table)
        c0 = conn0.cursor()
        c0.execute("select XL from MS2Data where mgf_file is '"+top_XL_file+"'")
        table0 = c0.fetchall()

        ms2_xl_list = []
        for row in table0:
            ms2_xl_list.append(str(row[0]))
        print ("XL selected from MS/MS analysis as the inout of pymol session\n")
        print (ms2_xl_list)


        ## 1. Reading and storing the pdb in a pose
        pose = Pose()
        pose_from_file(pose, pdb_filename)
        sequence = pose.sequence()
        pdb_info = pose.pdb_info()

        ## pymol initialization
        pymol.finish_launching() 
        p_name = pdb_filename[:-4]
        #cmd.do("set internal_gui_width, 250")
        cmd.bg_color("white")
        cmd.load(pdb_filename)
        cmd.set_title(p_name, 1, '')
        cmd.show_as("cartoon",p_name)
        cmd.color("gray",p_name+" and name CA")    

        print "\nThe protein sequence is: \n", sequence, "\n"


        ## 2. Reading the xl_file
        # with open(top_XL_file,'r') as xlfile:
        #     top_XL = xlfile.read().splitlines() # each rows as one element of the list

        output_xl_number = 0
        normal_dist_score = 0.0
        found_XL_list = []

        for num_xl,xl in enumerate(ms2_xl_list):
            print num_xl+1, xl

            K_pos_P1 = 0
            K_pos_P2 = 0
            eulidean_dist = 10000.0

            xl_without_digit = xl.translate(None, digits)
            peptide1 = xl_without_digit.split('--')[0].replace("-","").replace(".","")[:-2]
            peptide2 = xl_without_digit.split('--')[1].replace("-","")[:-3]

            K_pos_P1 = peptide1.find('K') + 1
            K_pos_P2 = peptide2.find('K') + 1

            multiple_occ_list1 = [x for x in xrange(len(sequence)) if sequence.find(peptide1, x) == x]
            multiple_occ_list2 = [y for y in xrange(len(sequence)) if sequence.find(peptide2, y) == y]

            ## finding minimum distance if multiple occurance happened
            seq_pos_p1_k = 0
            seq_pos_p2_k = 0
            tmp_dist = 10000.0
            for pos1 in multiple_occ_list1:
                for pos2 in multiple_occ_list2:
                    if tmp_dist > rosetta_eu_dist(pose, pos1+K_pos_P1, pos2+K_pos_P2):
                        tmp_dist = rosetta_eu_dist(pose, pos1+K_pos_P1, pos2+K_pos_P2)
                        seq_pos_p1_k = pos1+K_pos_P1
                        seq_pos_p2_k = pos2+K_pos_P2
                     
            print "aa positions on the sequence: ", seq_pos_p1_k, seq_pos_p2_k

            if ((peptide1 in sequence) and (peptide2 in sequence)):
                eulidean_dist = rosetta_eu_dist(pose, seq_pos_p1_k, seq_pos_p2_k)

                print "Euclidean distance is:  ", eulidean_dist, "\n"
            
            else:
                print "The XL is not found on the protein sequence. Check each peptide to be valid!\n"
            
            if eulidean_dist <= cut_off:
                output_xl_number += 1
                found_XL_list.append(xl)
                cmd.distance( "dist_"+str(num_xl+1), str(seq_pos_p1_k)+"/CA", str(seq_pos_p2_k)+"/CA")
                    

        ## writing XLs below threshold in a file
        for item in found_XL_list:
            found_xl_file.write("%s\n" % item)

        found_xl_file.close()

        cmd.set('dash_color', 'red')
        cmd.set('dash_width', 4)
        cmd.set('label_size', 22)
        pymol.cmd.save("pymol_result.pse")

        # Get out!
        pymol.cmd.quit()

        info['mapped_XL_on_TopModel'] = 'mapped_xl.txt'
        info['pymol_file'] = "pymol_result.pse"
        info['to_dropbox'].append("pymol_result.pse")
        info['to_dropbox'].append('mapped_xl.txt')

        return info


if __name__ == "__main__":
    pymol_run.main()

