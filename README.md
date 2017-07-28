#PyRosettaScripts   
#Authors: All members of the RosettaCommons  
#Date: 28/Apr/2017    
#Description: Transient repository for python scripts that use pyrosetta  
  
This is an early experiment, everything has to be figured out yet.  
  
Your python scripts MUST go into:  
  
pyRosettaScripts/designCategory/$USERNAME  
  
  
Tree-structure:  
 .  
 ├── README.md  
 ├── benchmarks  
 ├── doc  
 │   ├── examples  
 │   └── tutorials  
 ├── lib  
 ├── pyRosettaScripts  
 │   ├── motifGrafting  
 │   ├── peptideDesign  
 │   │   └── dbaker  
 │   │       ├── #make_peptides.py#  
 │   │       ├── hash_subclass.py  
 │   │       ├── hash_subclass.pyc  
 │   │       ├── make_peptides.py  
 │   │       ├── make_peptides.py~  
 │   │       ├── make_peptides_vary_repeat_length.py  
 │   │       ├── make_peptides_vary_repeat_length.py~  
 │   │       ├── pdb_utils_noclass.py  
 │   │       ├── pdb_utils_noclass.pyc  
 │   │       ├── repeat+xtl.list  
 │   │       ├── repeat_protein.list  
 │   │       ├── repeat_protein_xtls.list  
 │   │       ├── repeat_utils.py  
 │   │       ├── repeat_utils.pyc  
 │   │       └── repeat_utils.py~  
 │   ├── proteinIdealization  
 │   └── toolbox  
 ├── tests  
 └── tools  
  
NOTE 1 : It is not OK to just touch commas, spaces and other cosmetic details or   
         small bug-fixes and then add your name to the list of authors of some code.  
         Please do it only if you have a substantial contribution to such code development.  
  
