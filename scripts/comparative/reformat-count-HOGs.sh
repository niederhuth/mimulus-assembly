#!/bin/bash --login
#SBATCH --time=128:00:00
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
orthogroups=$(ls Results*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv)

#Reformat the N0.tsv
echo "Reformatting HOG N0.tsv"
if [ ${trim_p} = TRUE ]
then
	cut -f1,4- ${orthogroups} | sed s/\,\ /\;/g | tr '\t' ',' | sed s/\.p\,/\,/g | sed s/\.p\;/\;/g | tr ';' ' ' > orthogroups.csv
else
	cut -f1,4- ${orthogroups} | sed s/\,\ /\;/g | tr '\t' ',' | tr ';' ' ' > orthogroups.csv
fi

#Create the output count file with header line
head -1 ${orthogroups} | cut -f1,4- > orthogroup_counts.tsv

#Loop over each HOG and count
sed '1d' ${orthogroups} | while read line 
do
	og=$(echo ${line} | cut -d ' ' -f1)
	grep ${og} ${orthogroups} > tmp
	column=4
	#Loop over each species
	until [[ ${column} -gt $(expr ${col_num} + 4) ]]
	do
		#Count the genes for that species
		count=$(cut -f ${column} tmp | tr ',' '\n' | sed '/^$/d' | wc -l)
		og="${og} ${count}" 
		#Increase the column number by 1
		column=$(expr ${column} + 1)
	done
	#Output gene counts in each species for each orthogroup
	echo ${og} | tr ' ' '\t' >> orthogroup_counts.tsv
	rm tmp
done	

echo "Done"
