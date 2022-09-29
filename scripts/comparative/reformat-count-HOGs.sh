#!/bin/bash --login
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --job-name reformat-count-HOGs
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
orthofinder_dir= #Directory to orthogroup results, if left blank will look in current directory
trim_p=TRUE #TRUE/FALSE, trim the trailing ".p" from the protein name so matches cds name

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/phylo/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/phylo/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically

#Find orthofinder_dir and cd to it
if [ -z ${orthofinder_dir} ]
then
	orthofinder_dir="orthofinder"
	echo "Previous orthofinder results found: ${orthofinder_dir}"
	cd ${orthofinder_dir}
else
	cd ${orthofinder_dir}
	echo "Previous orthofinder results found: ${orthofinder_dir}"
fi

#Set path to the SpeciesID.txt file
speciesIDs=$(ls */WorkingDirectory/SpeciesIDs.txt)
#Count number of species
col_num=$(cat ${speciesIDs} | wc -l)
#Set path to Phylogenetic_Hierarchical_Orthogroups
orthogroups=$(ls */Phylogenetic_Hierarchical_Orthogroups/N0.tsv)
#Create the output files with header line
head -1 ${orthogroups} | cut -f1,4- > orthogroup_counts.tsv
echo "Transcript Orthogroup" | tr ' ' '\t' > transcript_orthogroup.tsv

#Reformat the N0.tsv
echo "Reformatting HOG N0.tsv"
cut -f1,4- ${orthogroups} | sed s/\,\ /\;/g | tr '\t' ',' | tr ';' ' ' > orthogroups.csv

#Loop over each HOG and get the relevant data
sed '1d' ${orthogroups} | while read line 
do
	og=$(echo ${line} | cut -d ' ' -f1)
	og_counts=${og}
	grep ${og} ${orthogroups} > tmp
	column=4
	#Loop over each species
	until [[ ${column} -gt $(expr ${col_num} + 4) ]]
	do
		#Get the genes for that species
		genes=$(cut -f ${column} tmp | sed 's/\ //g' | sed s/\,$//)
		#Trim off the ".p" found at end of some protein sequences
		if [ ${trim_p} = TRUE ]
		then
			genes=$(echo ${genes} | tr ',' '\n' | sed '/^$/d' | sed s/\.p$// | tr '\n' ' ')
		fi
		#Output each gene and its orthogroup
		for i in $(echo ${genes} | tr ',' '\n' | sed '/^$/d')
		do
			echo "${i} ${og}" | tr ' ' '\t' >> transcript_orthogroup.tsv
		done
		#Count the genes for that species
		count=$(echo ${genes} | tr ',' '\n' | sed '/^$/d' | wc -l)
		og_counts="${og_counts} ${count}" 
		#Increase the column number by 1
		column=$(expr ${column} + 1)
	done
	#Output gene counts in each species for each orthogroup
	echo ${og_counts} | tr ' ' '\t' >> orthogroup_counts.tsv
	rm tmp
done	

echo "Done"