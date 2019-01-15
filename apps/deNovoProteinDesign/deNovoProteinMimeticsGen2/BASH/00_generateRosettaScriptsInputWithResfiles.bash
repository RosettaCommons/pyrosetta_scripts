#!/bin/bash

#De Novo Protein Mimic Designer, by D.A.S.
#Citation: 
#  Daniel-Adriano Silva*, Shawn Yu*, Umut Y. Ulge*, Jamie B. Spangler, et. al., De novo design of potent and selective mimics of IL-2 and IL-15, Nature, 2019. https://doi.org/10.1038/s41586-018-0830-7

# OPTIONS START #
rosetta_bin="rosetta_scripts "
n_design_repeats=3
xml_file="${PWD}/XMLrosettaScripts/deNovoPreprofiledMimeticsDesign_gen2.xml"
rosetta_options="\
-beta \
-ex1 \
-ex2aro \
-overwrite \
-out:file:renumber_pdb false \
-ignore_zero_occupancy false \
-out:output -parser:protocol ${xml_file} \
-nstruct 3"
# OPTIONS END #

#DO NOT MODIFY BEYOND HERE
echo "#############################################################################################################################################"
echo "Script to generate a RosettaScripts setup for sequence design the PDB results from STEP1 (design of fully profiled de novo mimetic backbones)"
echo "Cite: Daniel-Adriano Silva*, Shawn Yu*, Umut Y. Ulge*, Jamie B. Spangler, et. al., De novo design of potent and selective mimics of IL-2 and IL-15, Nature, 2019." 
echo "      https://doi.org/10.1038/s41586-018-0830-7"
echo "#############################################################################################################################################"

if [ $# -ne 2 ]; 
then 
  echo "!!! Illegal number of parameters. It should be:"
  echo "     Arg1: Input directory with the PDB results from STEP1 (design of fully profiled de novo mimetic backbones) "
  echo "     Arg2: Output directory for the setup for rosettaScripts design "
  exit 1
fi

tasklist_file="${PWD}/${2}/designs.tasklist"
echo "#!/bin/bash" > ${tasklist_file}
for i in `ls ${PWD}/${1}/*.pdb`; 
do 
  echo "Generating RosettaScripts setup for file: ${i}"
  basename_pdb=`basename $i`; 
  this_design_dir="${PWD}/${2}/${basename_pdb}/";
  resfile_name="${basename_pdb}.resfile";
  mkdir -p ${this_design_dir}; 
  ln -s ${i} ${this_design_dir}; 
  grep "^REMARK PDBinfo-LABEL" ${i} | awk 'BEGIN{print "ALLAA\nSTART" }{outline=""; for(i=4;i<=NF;i++){if(substr($i,1,6)=="PIKAA_"){outline=outline substr($i,7,1)}}if(length(outline)>0){print $3 " A PIKAA " outline}}' > ${this_design_dir}/${resfile_name};
  for j in `seq 1 ${n_design_repeats}`
  do
    rad_id=`cat /dev/urandom | env LC_CTYPE=C tr -cd 'a-f0-9' | head -c 8`
    echo "cd ${this_design_dir} && ${rosetta_bin} ${rosetta_options} -s ${this_design_dir}/${basename_pdb} -out:path:pdb ${this_design_dir} -out:path:score ${this_design_dir} -parser:script_vars resfile_name=${this_design_dir}/${resfile_name} -suffix _${rad_id} >> design.log" >> ${tasklist_file}
  done
done
echo "All Done"

echo "Generation of the setup for rosettascripts is done. Run your designs!"
