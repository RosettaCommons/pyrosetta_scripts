## DDA_clean
## Filter out the input mgf file according to the list of XLs and produce a shorter version.

def mgf_writer(mgffile, spectra):

    mgffile.write("BEGIN IONS\n")
    mgffile.write("TITLE="+str(spectra['params']['title'])+"\n")
    mgffile.write("PEPMASS="+str(spectra['params']['pepmass'][0])+"\n")
    mgffile.write("RTINSECONDS="+str(spectra['params']['rtinseconds'])+"\n")
    mgffile.write("CHARGE="+(str(spectra['params']['charge'][0])).replace("+", "")+"\n")
    for m_z in range(len(spectra['m/z array'])):
        mgffile.write(str(spectra['m/z array'][m_z])+" "+str(spectra['intensity array'][m_z])+"\n")
    mgffile.write("END IONS\n")
    mgffile.write("\n")
    return

def DDA_clean(top_XL_file, mgf_file, PRE_CURSOR_DELTA, xlinker_mass):
                
    #from __future__ import division
    from fragment_generator import fragment_generator
    from pyteomics import mgf, mass

    print "Input information:", mgf_file, PRE_CURSOR_DELTA

    ## Read the MGF (MS/MS) file
    list_of_mgfDic = list(mgf.read(mgf_file))
      
    ## Importing top_XLs file
    with open(top_XL_file,'r') as xlfile:
        top_XL = xlfile.read().splitlines() # each rows as one element of the list
   
    output_file1_name = (str(mgf_file))[:-4]+"_cleaned.mgf"
    output_file1 = open(output_file1_name, 'w')

    for num_xl,xl in enumerate(top_XL):        

        print num_xl, xl    

        p1 = ""
        p2 = ""
        fragment_all = []
        mz_light_all = []
        mz_heavy_all = []
        precursor_dic = {}

        precursor_dic, fragment_all, mz_light_all, mz_heavy_all, p1, p2 = fragment_generator(xl, xlinker_mass)
        
        for specnum,anySpectra in enumerate(list_of_mgfDic):
            
            SpecPEPMASS = anySpectra['params']['pepmass'][0]
            SpecCharge = anySpectra['params']['charge'][0]

            if SpecCharge in range (3,9):

                mz_difference_Light = 1.000
                mz_difference_Heavy = 1.000
                mz_difference_Light = min(abs(SpecPEPMASS - pre_mz) for pre_mz in precursor_dic[SpecCharge][0:4])
                mz_difference_Heavy = min(abs(SpecPEPMASS - pre_mz) for pre_mz in precursor_dic[SpecCharge][5:9])

                if (mz_difference_Light <= PRE_CURSOR_DELTA):

                    mgf_writer(output_file1, anySpectra)

                elif (mz_difference_Heavy <= PRE_CURSOR_DELTA):

                    mgf_writer(output_file1, anySpectra)

    return output_file1_name

