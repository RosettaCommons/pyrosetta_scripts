def fragments(peptide, charge):

    from pyteomics import mgf, mass

    """
    The function generates all possible m/z for fragments of types 
    `types` and of charges from 1 to `maxharge`.
    """
    pep_mass = 0.0
    fragment_list = []
    mZ_list = []
    types=('b', 'y')

    for i in xrange(1, len(peptide)):

        for ion_type in types:

            if ion_type[0] in 'abc':

                if "C" in peptide: # adding Carbamidomethylation for Cys residues
                    pep_mass = mass.fast_mass(peptide[:i], ion_type=ion_type, charge=charge) \
                    + (57.021464*peptide[:i].count("C"))
                else:
                    pep_mass = mass.fast_mass(peptide[:i], ion_type=ion_type, charge=charge)
                fragment_list.append(peptide[:i])
                mZ_list.append(pep_mass)

            else:

                if "C" in peptide: # adding Carbamidomethylation for Cys residues
                    pep_mass = mass.fast_mass(peptide[i:], ion_type=ion_type, charge=charge) \
                    + (57.021464*peptide[i:].count("C"))
                else:
                    pep_mass = mass.fast_mass(peptide[i:], ion_type=ion_type, charge=charge)

                fragment_list.append(peptide[i:])
                mZ_list.append(pep_mass)

    return fragment_list,mZ_list
