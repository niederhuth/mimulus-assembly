#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --job-name count-HOGs
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
orthofinder_dir= #Directory to orthogroup results, if left blank will look in current directory

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/phylo/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/phylo/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically

#Find orthofinder_dir and cd to it
if [ -z ${orthofinder_dir} ]
then
	orthofinder_dir=$(ls orthofinder/Results_*)
	echo "Previous orthofinder results found: ${orthofinder_dir}"
	cd ${orthogroup_dir}
else
	cd ${orthofinder_dir}
	echo "Previous orthofinder results found: ${orthofinder_dir}"
fi

#Set path to the SpeciesID.txt file
if [ -z ${speciesIDs} ]
then
	speciesIDs=$(ls WorkingDirectory/SpeciesIDs.txt)
	col_num=$(cat ${speciesIDs} | wc -l)
fi
#Set path to N0.tsv orthofinder output
if [ -z ${orthogroups} ]
then
	orthogroups=$(ls Phylogenetic_Hierarchical_Orthogroups/N0.tsv)
fi

#Loop over each HOG and get the relevant data
while read line 
do
	og=$(echo ${line} | cut -d ' ' -f1)
	column=4
	#Loop over each species
	until [[ ${column} -gt ${col_num} ]]
	do
		#Count the genes for that species
		count=$(cut -f ${column} tmp | sed 's/\ //g' | sed '/^$/d' | wc -l)
		og="${og} ${count}" 
		#Increase the column number by 1
		column=$(expr ${column} + 1)
	done
	echo ${og} | tr ' ' '\t' >> HOG_counts.tsv
done	

echo "Done"