def fragment_generator(xl, DSS_mass):

    from pyteomics import mgf, mass
    from fragments import fragments
    from split_xl import split_xl

    """
    input should be a XL in kojak format.
    example: -.PEPKTIDER(4)--PEPKTIDER(4).-
    """

    ## constants
    AA_string = "ARNDCQEGHILKMFPSTWYV-"
    DSS_string = "X"
    #DSS_mass = 138.06808 # mass.calculate_mass(formula='C16H20N2O8') #disuccinimidyl suberate
    #DSS_mass_H = 138.06808 + 12*1.006276746
    Mass_of_H = 1.008

    fragment_list1 = []
    mZ_list1 = []
    fragment_list2 = []
    mZ_list2 = []
    fragment_list = []
    mZ_list = []

    comb_frag1 = []
    comb_frag2 = []

    point = 0
    K_pos_P1 = 0
    K_pos_P2 = 0

    temp = []
    temp1 = []
    temp2 = []

    for item in xl:
        if item in AA_string:
            temp.append(item)
            point = split_xl(temp)

    temp1 = temp[1:point]
    temp2 = temp[point+2:len(temp)-1]
    peptide1 = ''.join(temp1)
    peptide2 = ''.join(temp2)

    K_pos_P1 = peptide1.find('K')
    K_pos_P2 = peptide2.find('K')

    ## for filtering spectra based on comparison of mass of XL and the precursor mass
    ## precursor charge range is 3 to 8
    ## ==================================================
    precursor_dic = {}

    precursor_frag_list = []
    add_to_frag_list_L = []
    add_to_frag_list_H = []

    for precursor_charge in range (3,9):

        precursor_mZ_list_H = []
        precursor_mZ_list_L = []

        precursor_charge_letter = ''.join(['+' for s in xrange(precursor_charge)])

        precursor_frag_list.append(peptide1 + peptide2 + precursor_charge_letter)

        XLMASS = mass.calculate_mass(sequence=peptide1) +\
            DSS_mass+\
            mass.calculate_mass(sequence=peptide2)

        # adding Carbamidomethylation for Cys residues
        XLMASS = XLMASS + (57.021464*(peptide1.count("C")+peptide2.count("C")))
        XLMASS = (XLMASS + (precursor_charge * Mass_of_H))/precursor_charge
        #print "Precursor monoisotopic M/Z", XLMASS

        ## calculating different isotopic mass/z for the precursor
        precursor_mZ_list_L.append(XLMASS)
        precursor_mZ_list_L.append(XLMASS + float(1/precursor_charge))
        precursor_mZ_list_L.append(XLMASS + float(2/precursor_charge))
        precursor_mZ_list_L.append(XLMASS + float(3/precursor_charge))
        precursor_mZ_list_L.append(XLMASS + float(4/precursor_charge))

        ## calculating M heavy (DSS_D12)
        XLMASS_heavy = XLMASS + ((12*1.006276746)/precursor_charge)
        precursor_mZ_list_H.append(XLMASS_heavy)

        precursor_mZ_list_H.append(XLMASS_heavy + float(1/precursor_charge))
        precursor_mZ_list_H.append(XLMASS_heavy + float(2/precursor_charge))
        precursor_mZ_list_H.append(XLMASS_heavy + float(3/precursor_charge))
        precursor_mZ_list_H.append(XLMASS_heavy + float(4/precursor_charge))

        ## adding the list's to the dictionary
        ## 0,1,2,3,4 = Light indices
        ## 5,6,7,8,9 = Heavy indices

        precursor_dic[precursor_charge] = precursor_mZ_list_L + precursor_mZ_list_H

        ## here we store light and heavy m/z (each list 6 numbers for different charge values)
        ## to add them below to the fragment list
        add_to_frag_list_L.append(XLMASS)
        add_to_frag_list_H.append(XLMASS_heavy)

    #print "precursor_frag_list", precursor_frag_list
    #print "add_to_frag_list_L", add_to_frag_list_L
    #print "add_to_frag_list_H", add_to_frag_list_H

    ## ==================================================
    ## lists to store fragment's and m/z for all charge values
    fragment_list = []
    mZ_list_H = []
    mZ_list_L = []
    mZ_list = []

    for charge in range (1,5):

        #print "======================================================================================"

        #print "CHARGE VALUE IS:", charge
        fragment_list1,mZ_list1 = fragments(peptide1, charge)
        fragment_list2,mZ_list2 = fragments(peptide2, charge)

        charge_letter = ''.join(['+' for s in xrange(charge)])
        #print charge_letter
        ## charge 1 should just consider fragments without XL
        #if charge == 1:

        # print fragment_list1
        # print mZ_list1
        # print fragment_list2
        # print mZ_list2

        ## considering the combination of fragments with full peptide
        ## FRAG1 + XL + PEPTIDE2
        fragment_list_temp1 = list(fragment_list1)
        for num1,frag1 in enumerate(fragment_list_temp1):
            if (peptide1.find(frag1) <= K_pos_P1) and ((peptide1.find(frag1) + len(frag1) - 1) >= K_pos_P1):
                fragment_list.append(frag1 + peptide2 + charge_letter)
                comb_frag1.append(frag1)
                comb_mZ = mass.calculate_mass(sequence=frag1) - 18.010565 + DSS_mass + mass.calculate_mass(sequence=peptide2)
                comb_mZ = comb_mZ + (57.021464*(frag1.count("C")+peptide2.count("C")))
                comb_mZ = (comb_mZ + (charge * Mass_of_H))/charge
                mZ_list_L.append(comb_mZ)
                mZ_list_H.append(comb_mZ + float((12*1.006276746)/charge))
            else:
                fragment_list.append(frag1 + charge_letter)
                mZ_list_L.append(mZ_list1[num1])
                mZ_list_H.append(mZ_list1[num1])


        ## FRAG1 + XL + PEPTIDE2
        fragment_list_temp2 = list(fragment_list2)
        for num2,frag2 in enumerate(fragment_list_temp2):
            if (peptide2.find(frag2) <= K_pos_P2) and (peptide2.find(frag2) + len(frag2) - 1) >= K_pos_P2:
                fragment_list.append(frag2 + peptide1 + charge_letter)
                comb_frag2.append(frag2)
                comb_mZ = mass.calculate_mass(sequence=frag2) - 18.010565 + DSS_mass + mass.calculate_mass(sequence=peptide1)
                comb_mZ = comb_mZ + (57.021464*(frag2.count("C")+peptide1.count("C")))
                comb_mZ = (comb_mZ + (charge * Mass_of_H))/charge
                mZ_list_L.append(comb_mZ)
                mZ_list_H.append(comb_mZ + float((12*1.006276746)/charge))
            else:
                fragment_list.append(frag2 + charge_letter)
                mZ_list_L.append(mZ_list2[num2])
                mZ_list_H.append(mZ_list2[num2])


        ## FRAG1 + XL + FRAG2
        ## here we make all combination of fragments that contain xlinker arm
        # for item1 in comb_frag1:
        #     for item2 in comb_frag2:
        #         if item1+item2 not in fragment_list:
        #             print item1, '+' , item2
        #             fragment_list.append(item1+item2+charge_letter)
        #             comb_mZ = mass.calculate_mass(sequence=item1) - 18.010565 + DSS_mass + mass.calculate_mass(sequence=item2) - 18.010565
        #             comb_mZ = comb_mZ + (57.021464*(item1.count("C")+item2.count("C")))
        #             comb_mZ = (comb_mZ + (charge * Mass_of_H))/charge
        #             mZ_list.append(comb_mZ)

    fragment_list.extend(precursor_frag_list)
    mZ_list_L.extend(add_to_frag_list_L)
    mZ_list_H.extend(add_to_frag_list_H)

    #print "fragment_list_full:\n", fragment_list
    #print "mZ_list_Light:\n", mZ_list_L
    #print "mZ_list_Heavy:\n", mZ_list_H


    return precursor_dic, fragment_list, mZ_list_L, mZ_list_H, peptide1, peptide2
