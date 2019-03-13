#!/usr/bin/env python
"""
Created on Fri Dec  9 11:28:23 2016

Modified on Tue Sep 5 17:11:00 2017

@author: hamedkhakzad
"""

## TaxLink
## Search computational XLs in MSMS data (mgf format) to find and analyze real patterns.

from __future__ import division
from pyteomics import mgf, mass
from fragment_generator import fragment_generator
from fig_maker import fig_maker
from DDA_clean import DDA_clean
import multiprocessing as mp
import argparse
import sqlite3


def taxlink(top_XL_file, mgf_file, delta, intensity_filter):
    ## input arguments
    ## COMMAND: ./MSMS_analysis.py -x top_xl.txt -m file.mgf -d 0.05
    # parser = argparse.ArgumentParser()
    # parser.add_argument("-x","--xlink", help="input top_xl file")
    # parser.add_argument("-m","--mgf", help="input mgf file")
    # parser.add_argument("-d","--delta", help="delta window for fragment detection; default=0.01", type=float, default=0.01)
    # parser.add_argument("-p","--predelta", help="delta for pre-cursor m/z filtering; default=0.01", type=float, default=0.01)
    # parser.add_argument("-f","--figthreshold", help="generate png for pectra with this number of detected fragments; default=20", type=int, default=20)
    # parser.add_argument("-i","--intensity", help="intensity threshold", type=float, default=0.0)
    # parser.add_argument("-l","--linkertype", help="xlinker type; DSS, DSG, or EGS", default='DSS')
    # args = parser.parse_args()
    
    # top_XL_file = args.xlink
    # mgf_file = args.mgf
    # delta = args.delta
    # intensity_filter = args.intensity
    # xlinker_type = args.linkertype
    # PRE_CURSOR_DELTA = args.predelta
    # fig_threshodl = args.figthreshold
    
    DSS_mass = 138.06808 # mass.calculate_mass(formula='C16H20N2O8') #disuccinimidyl suberate
    DSG_mass = 96.02059
    EGS_mass = 226.04774

    # if xlinker_type == 'EGS':
    #     xlinker_mass = EGS_mass
    # elif xlinker_type == 'DSG':
    #     xlinker_mass = DSG_mass
    # else:
    #     xlinker_mass = DSS_mass

    xlinker_mass = DSS_mass
    fig_threshodl = 10
    final_score = 0
    PRE_CURSOR_DELTA = delta

    ## PRE-PROCESSING STAGE
    ## Considering the list of input XLs, we remove all extra spectra
    ## from input MGF file according to the pre-cursor m/z and make
    ## a new 'cleaned version' MGF file to analyze.
    ## This step makes the whole process a lot faster!
    print "cleaning stage! ..."
    print "pre-cursor delta is: ", PRE_CURSOR_DELTA
    cleaned_mgf_file = DDA_clean(top_XL_file, mgf_file, PRE_CURSOR_DELTA, xlinker_mass)
    
    print "cleaning is done!\n"
    

    ## SQLITE3 database
    con = sqlite3.connect('MS2_results.db')
    c = con.cursor()
    # Create table
    c.execute('''CREATE TABLE IF NOT EXISTS MS2Data
                 (XL text, mgf_file text, spectrum_id text, spectrum_num integer, \
                 delta real, pre_charge integer, H_L integer, fragSc integer, \
                 coverage real, covered_Frags text, covered_Mz text, covered_int text,\
                  main_Mz text, main_int text, count integer)''')
    c.execute('PRAGMA journal_mode = WAL')
    c.execute('PRAGMA synchronous = NORMAL')
            
    ## Read the MGF (MS/MS) file
    list_of_mgfDic = list(mgf.read(cleaned_mgf_file))
      
    ## Importing top_XLs file
    with open(top_XL_file,'r') as xlfile:
        top_XL = xlfile.read().splitlines() # each rows as one element of the list
   
    output_file1_name = "detected_spectra.txt"
    output_file1 = open(output_file1_name, 'w')

    for num_xl,xl in enumerate(top_XL):        

        print "====================================================================================="
        print num_xl, xl    

        p1 = ""
        p2 = ""
        fragment_all = []
        mz_light_all = []
        mz_heavy_all = []
        precursor_dic = {}

        precursor_dic, fragment_all, mz_light_all, mz_heavy_all, p1, p2 = fragment_generator(xl, xlinker_mass)
    
        mgf_spec_score = [0] * len(list_of_mgfDic)

        ## Considering all fragments and find the match peaks in MS/MS data
        mgf_spec_list = [None] * len(mz_light_all)
        output_file1.write(xl+"\n") 
        
        for specnum,anySpectra in enumerate(list_of_mgfDic):
            
            SpecPEPMASS = anySpectra['params']['pepmass'][0]
            SpecCharge = anySpectra['params']['charge'][0]

            if SpecCharge in range (3,9):

                mz_difference_Light = 1.000
                mz_difference_Heavy = 1.000
                mz_difference_Light = min(abs(SpecPEPMASS - pre_mz) for pre_mz in precursor_dic[SpecCharge][0:4])
                mz_difference_Heavy = min(abs(SpecPEPMASS - pre_mz) for pre_mz in precursor_dic[SpecCharge][5:9])

                if (mz_difference_Light <= PRE_CURSOR_DELTA): # we fixed delta for precursor to 0.01

                    output_file1.write(anySpectra['params']['title']+"\n")
                    for num_mz,anymZ in enumerate(mz_light_all):

                        for y in range(len(anySpectra['m/z array'])):
                            min_MZ = abs(anymZ - anySpectra['m/z array'][y])
                            if min_MZ <= delta:
                                maxintensity = max(anySpectra['intensity array'][0:])
                                mgf_spec_list[num_mz] = anySpectra
                                
                                ## intensity-based scoring of the detected peaks
                                if anySpectra['intensity array'][y]>= maxintensity*3/4:
                                    mgf_spec_score[specnum] += 12
                                elif anySpectra['intensity array'][y]>= maxintensity/2:
                                    mgf_spec_score[specnum] += 8
                                elif anySpectra['intensity array'][y]>= maxintensity/4:
                                    mgf_spec_score[specnum] += 4
                                elif anySpectra['intensity array'][y]>=intensity_filter:
                                        mgf_spec_score[specnum] += 2

                elif (mz_difference_Heavy <= PRE_CURSOR_DELTA): # we fixed delta for precursor to 0.01

                    output_file1.write(anySpectra['params']['title']+"\n")
                    for num_mz,anymZ in enumerate(mz_heavy_all):

                        for y in range(len(anySpectra['m/z array'])):
                            min_MZ = abs(anymZ - anySpectra['m/z array'][y])
                            if min_MZ <= delta:
                                maxintensity = max(anySpectra['intensity array'][0:])
                                mgf_spec_list[num_mz] = anySpectra
                                
                                ## intensity-based scoring of the detected peaks
                                if anySpectra['intensity array'][y]>= maxintensity*3/4:
                                    mgf_spec_score[specnum] += 12
                                elif anySpectra['intensity array'][y]>= maxintensity/2:
                                    mgf_spec_score[specnum] += 8
                                elif anySpectra['intensity array'][y]>= maxintensity/4:
                                    mgf_spec_score[specnum] += 4
                                elif anySpectra['intensity array'][y]>=intensity_filter:
                                        mgf_spec_score[specnum] += 2
        
        spectrum_id = "NA"
        spectrum_num = 0
        peaks_num = 0
        coverage = 0.0

        for num_mgf,item_mgf in enumerate(list_of_mgfDic):
            if max(mgf_spec_score) != 0:
                if mgf_spec_score[num_mgf] == max(mgf_spec_score):
                    heavy_light_flag = ""
                    SpecPEPMASS = item_mgf['params']['pepmass'][0]
                    SpecCharge = item_mgf['params']['charge'][0]
                    print "SpecPEPMASS, SpecCharge", SpecPEPMASS, SpecCharge
                    mz_difference_Light = 1.000
                    mz_difference_Heavy = 1.000
                    mz_difference_Light = min(abs(SpecPEPMASS - pre_mz) for pre_mz in precursor_dic[SpecCharge][0:4])
                    mz_difference_Heavy = min(abs(SpecPEPMASS - pre_mz) for pre_mz in precursor_dic[SpecCharge][5:9])

                    spectrum_id = item_mgf['params']['title']
                    spectrum_num = num_mgf
                    peaks_num = len(item_mgf['m/z array'])

                    ## to send to fig_maker
                    main_spectra = item_mgf['m/z array']
                    main_intensity = item_mgf['intensity array']
                    
                    covered_Frags = []
                    covered_mz = []
                    covered_intensity = []

                    if (mz_difference_Light <= PRE_CURSOR_DELTA): # we fixed delta for precursor to 0.01

                        heavy_light_flag = "Light"
                        for numz,mz in enumerate(mz_light_all):
                            min_MZ2 = min(abs(mz - y) for y in item_mgf['m/z array'][0:])
                            if min_MZ2 <= delta:
                                covered_Frags.append(fragment_all[numz])
                                covered_mz.append(mz_light_all[numz])

                    elif (mz_difference_Heavy <= PRE_CURSOR_DELTA): # we fixed delta for precursor to 0.01

                        heavy_light_flag = "Heavy"
                        for numz,mz in enumerate(mz_heavy_all):
                            min_MZ2 = min(abs(mz - y) for y in item_mgf['m/z array'][0:])
                            if min_MZ2 <= delta:
                                covered_Frags.append(fragment_all[numz])
                                covered_mz.append(mz_heavy_all[numz])

                    for y_num,ymz in enumerate(covered_mz):    
                        for mz_num,anymz in enumerate(main_spectra):
                            min_MZ = abs(anymz - ymz)
                            if min_MZ <= delta:
                                if main_intensity[mz_num] != 0.0:
                                    covered_intensity.append(main_intensity[mz_num])
                                    break
                        
                    if peaks_num != 0:
                        coverage = (len(covered_Frags) * 100)/peaks_num
                    
                    if (spectrum_id != "NA" and len(covered_mz) >= 5):
                        c.execute("INSERT INTO MS2Data(XL, mgf_file,spectrum_id, spectrum_num, \
                        delta, pre_charge, H_L, fragSc, coverage, covered_Frags, covered_Mz, covered_int, \
                        main_Mz, main_int, count) VALUES \
                        (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)" \
                        ,(xl, top_XL_file, spectrum_id, spectrum_num, delta, SpecCharge, heavy_light_flag, \
                        max(mgf_spec_score), coverage, ",".join(str(x) for x in covered_Frags), \
                        ",".join(str(w) for w in covered_mz), ",".join(str(x) for x in covered_intensity), \
                        ",".join(str(x) for x in main_spectra), ",".join(str(x) for x in main_intensity), \
                        len(covered_Frags) ) )
                        
                        con.commit()

                    if len(covered_mz) >= fig_threshodl: ## here we can set a threshold for number of figures.
        
                        fig_maker(main_spectra, main_intensity, covered_Frags, covered_mz, p1, p2, xl, num_mgf, delta)
                        final_score += 1

    output_file1.close()
    con.close()
    return final_score

if __name__ == '__main__':
    
    main ()


