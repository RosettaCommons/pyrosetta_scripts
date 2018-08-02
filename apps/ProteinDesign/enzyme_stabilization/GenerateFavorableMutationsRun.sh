#!/bin/bash

read -p "Enter 'b' to process pdb files in a batch. If you'd like to process one pdb at a time, press enter: " batch
if [[ ( $batch == "b" ) || ( $batch == "B" )]]
then
	read -p "Enter the project name: " project_name
	mkdir $project_name
	cd $project_name
	marker=1
	while [ $((marker)) == 1 ]
	do 
		read -p "Enter whether distances to active site and contact numbers have been generated (True or False): " distandcn
		if [[ ( $distandcn == "True" ) || ( $distandcn == "T" ) || ( $distandcn == "t" ) || ( $distandcn == "true" ) || ( $distandcn == "False" ) || ( $distandcn == "false" ) || ( $distandcn == "F" ) || ( $distandcn == "f" ) ]]
		then
			marker=2
		else
			echo "Invalid input given. Please check your spelling."
			marker=1
		fi
	done
	if [[ ( $distandcn == "False" ) || ( $distandcn == "F") || ( $distandcn == "f" ) || ( $distandcn == "false" ) ]]
	then
		bash ../ProcessMutationsBatch.sh
	fi
	if (( $? != 0 )) ; then
		echo "ERROR IN CALCULATING THE DISTANCE TO ACTIVE SITE AND CONTACT NUMBER OF THE BATCH PROVIDED. Please check your spelling for provided filenames. Exiting..."
		exit
	fi
	cd ..
	read -p "Enter PSSM score threshold or press enter to accept the default (0): " pssm_threshold
	if [ "$pssm_threshold" == "" ]
	then
		pssm_threshold=0	
	fi
	read -p "Enter the distance to active site threshod or press enter to accept the default (15.0): " dist_threshold
	if [ "$dist_threshold" == "" ]
	then
		dist_threshold=15.0	
	fi
	read -p "Enter the contact number threshold or press enter to accept the default (16): " contact_threshold
	if [ "$contact_threshold" == "" ]
		then
		contact_threshold=16
	fi
	distance_threshold=$(echo $dist_threshold | bc)
	marker=0
	while [ $((marker)) == 0 ]
	do
		read -p "Enter whether PSSM scores have been generated (True or False): " scores
		if [[ ( $scores == "True" ) || ( $scores == "T" ) || ( $scores == "t" ) || ( $scores == "true" ) || ( $scores == "False" ) || ( $scores == "false" ) || ( $scores == "F" ) || ( $scores == "f" ) ]]
		then
			marker=1
		else
			echo "Invalid input given. Please check your spelling."
			marker=0
		fi
	done
	if [[ ( $scores == "False" ) || ( $scores == "F") || ( $scores == "f" ) || ( $scores == "false" ) ]] 
	then
		echo NOW GENERATING PSSM SCORES
		cd PSSM
		bash PSSMPipelineRun.sh 
		if (( $? != 0 )) ; then
			echo "ERROR DETECTED IN CALCULATING PSSM SCORES. Terminating..."
			exit
		fi
		pssm_filename=PSSM_0001.csv
		cd ..
	fi
	mv PSSM_0001.csv $project_name 
	cd $project_name
	python ../generate_favorable_mutations.py MutationsinBatch -s $((pssm_threshold)) -w $dist_threshold -l $((contact_threshold))
	cd ..
	exit
fi
read -p "Enter 'x' to generate the fasta sequence from a pdb file to use in a BlastP search. If you already have the fasta sequence in the PSSM folder, press enter to continue: " fasta
if [[ ( $fasta == "x" ) || ( $fasta == "X" ) ]]
then
	read -p "Enter the pdb name in the form 'filename.pdb': " pdbname
	python ../Rosetta/main/source/scripts/python/public/pdb2fasta.py $pdbname > $pdbname.fasta
	head -n2 $pdbname.fasta
	echo " "
	rm $pdbname.fasta
	exit
fi
marker=0
while [ $((marker)) == 0 ]
do
	read -p "Enter whether PSSM scores have been generated (True or False): " scores
	if [[ ( $scores == "True" ) || ( $scores == "T" ) || ( $scores == "t" ) || ( $scores == "true" ) || ( $scores == "False" ) || ( $scores == "false" ) || ( $scores == "F" ) || ( $scores == "f" ) ]]
	then
		marker=1
	else
		echo "Invalid input given. Please check your spelling."
		marker=0
	fi
done 
while [ $((marker)) == 1 ]
do 
	read -p "Enter whether distances to active site and contact numbers have been generated (True or False): " distandcn
	if [[ ( $distandcn == "True" ) || ( $distandcn == "T" ) || ( $distandcn == "t" ) || ( $distandcn == "true" ) || ( $distandcn == "False" ) || ( $distandcn == "false" ) || ( $distandcn == "F" ) || ( $distandcn == "f" ) ]]
	then
		marker=2
	else
		echo "Invalid input given. Please check your spelling."
		marker=1
	fi
done
read -p "Enter PSSM score threshold or press enter to accept the default (0): " pssm_threshold
if [ "$pssm_threshold" == "" ]
then
	pssm_threshold=0	
fi
read -p "Enter the distance to active site threshod or press enter to accept the default (15.0): " dist_threshold
if [ "$dist_threshold" == "" ]
then
	dist_threshold=15.0	
fi
read -p "Enter the contact number threshold or press enter to accept the default (16): " contact_threshold
if [ "$contact_threshold" == "" ]
	then
	contact_threshold=16
fi
distance_threshold=$(echo $dist_threshold | bc)

if [[ ( $scores == "True" ) || ( $scores == "T") || ( $scores == "t" ) || ( $scores == "true" ) ]]
then
	if [[ ( $distandcn == "True" ) || ( $distandcn == "T") || ( $distandcn == "t" ) || ( $distandcn == "true" ) ]] 
	then
		read -p "Enter the project name: " project_name
		mkdir $project_name
		read -p "Enter the name of the file containing the PSSM scores: " pssm_filename
		read -p "Enter the name of the file containing the distance to the active site and contact number of each residue: " dist_and_contact_filename
		read -p "Enter the desired name for the favorable mutations output file or press enter to accept the default (favorable_mutations.csv): " outfilename
		if [ "$outfilename" == '' ]
		then
			outfilename="favorable_mutations.csv"
		fi
		mv $dist_and_contact_filename $project_name
		cd $project_name
		python ../generate_favorable_mutations.py GenerateMutations -d $distance_threshold -p $((pssm_threshold)) -c $((contact_threshold)) -o $outfilename -m ../$pssm_filename -n $dist_and_contact_filename
		cd ..
		exit
	fi
fi
if [[ ( $scores == "True" ) || ( $scores == "T") || ( $scores == "t" ) || ( $scores == "true" ) ]] 
then
	if [[ ( $distandcn == "False" ) || ( $distandcn == "F") || ( $distandcn == "f" ) || ( $distandcn == "false" ) ]]
	then
		read -p "Enter the project name: " project_name
		mkdir $project_name
		read -p "Enter the name of the file containing the PSSM scores: " pssm_filename
		read -p "Enter the desired name for the distance and contact number output file or press enter to accept the default (results.csv): " distcn_outfilename
		read -p "Enter the desired name for the favorable mutations output file or press enter to accept the default (favorable_mutations.csv): " outfilename
		if [ "$outfilename" == '' ]
		then
			outfilename="favorable_mutations.csv"
		fi
		if [ "$distcn_outfilename" == "" ]
		then 
			distcn_outfilename="results.csv"
		fi
		echo NOW CALCULATING THE DISTANCE TO ACTIVE SITE AND CONTACT NUMBER
		cd $project_name
		python ../generate_favorable_mutations.py CalculateDistAndCN -o $distcn_outfilename
		cd ..
		if (( $? != 0 )) ; then
			echo "ERROR DETECTED IN CALCULATING DISTANCE TO ACTIVE SITE AND CONTACT NUMBER. Terminating..."
			exit
		fi
		echo NOW GENERATING FAVORABLE MUTATIONS
		cd $project_name
		python ../generate_favorable_mutations.py GenerateMutations -p $((pssm_threshold)) -d $distance_threshold -c $((contact_threshold)) -o $outfilename -m ../$pssm_filename -n $distcn_outfilename
		cd ..
		exit
	fi
fi
if [[ ( $scores == "False" ) || ( $scores == "F") || ( $scores == "f" ) || ( $scores == "false" ) ]]
then
	if [[ ( $distandcn == "True" ) || ( $distandcn == "T") || ( $distandcn == "t" ) || ( $distandcn == "true" ) ]]
	then
		read -p "Enter the project name: " project_name
		mkdir $project_name
		read -p "Enter the name of the file containing the distance to active site and contact number for each residue: " dist_and_contact_filename
		read -p "Enter the desired name for the favorable mutation output file or press enter to accept the default (favorable_mutations.csv): " outfilename
		if [ "$outfilename" == "" ]
		then 
			outfilename="favorable_mutations.csv"
		fi
		echo NOW GENERATING PSSM SCORES
		cd PSSM
		bash PSSMPipelineRun.sh 
		if (( $? != 0 )) ; then
			echo "ERROR DETECTED IN CALCULATING PSSM SCORES. Terminating..."
			exit
		fi
		echo NOW GENERATING FAVORABLE MUTATIONS
		cd ..
		mv $dist_and_contact_filename $project_name
		cd $project_name
		python ../generate_favorable_mutations.py GenerateMutations -p $((pssm_threshold)) -d $distance_threshold -c $((contact_threshold)) -o $outfilename -m ../PSSM_0001.csv -n $dist_and_contact_filename
		cd ..
		exit
	fi
fi
if [[ ( $scores == "False" ) || ( $scores == "F") || ( $scores == "f" ) || ( $scores == "false" ) ]]
then
	if [[ ( $distandcn == "False" ) || ( $distandcn == "F") || ( $distandcn == "f" ) || ( $distandcn == "false" ) ]]
	then
		read -p "Enter the project name: " project_name
		mkdir $project_name
		read -p "Enter the desired name for the distance to active site and contact number output file or press enter to accept the default (results.csv): " dist_and_contact_filename
		if [ "$dist_and_contact_filename" == "" ]
		then
			dist_and_contact_filename="results.csv"
		fi
		read -p "Enter the desired name for the favorable mutations output file or press enter to accept the default (favorable_mutations.csv): " outfilename
		if [ "$outfilename" == "" ]
		then
			outfilename="favorable_mutations.csv"
		fi
		echo NOW CALCULATING THE DISTANCE TO ACTIVE SITE AND CONTACT NUMBER
		cd $project_name
		python ../generate_favorable_mutations.py CalculateDistAndCN -o $dist_and_contact_filename
		cd ..
		if (( $? != 0 )) ; then
			echo "ERROR DETECTED IN CALCULATING DISTANCE TO ACTIVE SITE AND CONTACT NUMBER. Terminating..."
			exit
		fi
		echo NOW GENERATING PSSM SCORES
		cd PSSM
		bash PSSMPipelineRun.sh
		if (( $? != 0 )) ; then
			echo "ERROR DETECTED IN CALCULATING PSSM SCORES. Terminating..."
			exit
		fi
		echo NOW GENERATING FAVORABLE MUTATIONS
		cd ..
		cd $project_name
		python ../generate_favorable_mutations.py GenerateMutations -p $((pssm_threshold)) -d $distance_threshold -c $((contact_threshold)) -o $outfilename -m ../PSSM_0001.csv -n $dist_and_contact_filename
		cd ..
		exit
	fi
fi
if [[ ( ( $scores != "True" ) || ( $scores != "T") || ( $scores != "t" ) || ( $scores != "true" ) ) && ( ( $distandcn != "False" ) || ( $distandcn != "F") || ( $distandcn != "f" ) || ( $distandcn != "false" ) ) ]]
then
	if [[ ( ( $scores != "False" ) || ( $scores != "F") || ( $scores != "f" ) || ( $scores != "false" ) ) && ( ( $distandcn != "True" ) || ( $distandcn != "T") || ( $distandcn != "t" ) || ( $distandcn != "true" ) ) ]]
	then
		echo "One or more invalid answers given. Acceptable answers: 'True', 'T', 't', 'true', 'False', 'F', 'f', 'false' "
	fi
fi
	