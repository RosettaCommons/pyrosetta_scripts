The script, once running will ask you for all the information it needs. It is designed such tthat someone less knowledgable about code need not have to deal with command line issues. However, to make it this way, a lot of assumptions had to be made, such as a particular directory configuration. It works by running GenerateFavorableMutationsRun.sh.

All involved scripts are:
GenerateFavorableMutationsRun.sh
generate_favorable_mutations.py
contact_num.xml
ProcessMutationsBatch.sh

The PSSM folder containing all associated scripts for running J. Klesmith's method for calculate PSSM is also required, and it must be located in the same directory as the above scripts. 

Rosetta and PyRosetta are also required to be installed in the directory above these scripts and the PSSM folder. 