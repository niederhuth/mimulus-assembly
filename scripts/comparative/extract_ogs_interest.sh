
#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=400GB
#SBATCH --job-name anchorwave
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
target_species="Mguttatus_L1 Mguttatus_S1"
goi= #
speciesIDs= #
orthogroups= #
seqs= #
datatype="proteins-primary" #Set this to datatype used for orthofinder, normally this will be proteins-primary

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/orthofinder/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/orthofinder/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
path2=$(pwd)
path3=goi

#
if [ -z ${goi} ]
then
	goi=${path1}/goi.csv
fi
#
if [ -z ${speciesIDs} ]
then
	speciesIDs="${path2/}orthofinder/*/WorkingDirectory/SpeciesIDs.txt"
fi
#
if [ -z ${orthogroups} ]
then
	orthogroups="${path2}/orthofinder/*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"
fi
#
if [ -z ${seqs} ]
then
	seqs="${path2}/${datatype}/${a}.fa"

#Check for and make/cd working directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#
while read line 
do
	name=$(echo ${line} | cut -d ' ' -f1)
	ref_genes=$(echo ${line} | cut -d ' ' -f2)
	mkdir ${name}
	#Get the orthogroups
	for x in ${ref_genes}
	do
		grep ${x} ${orthogroups} >> tmp
	done
	#First species starts on column 4
	column=4
	#Loop over each species
	cat ${speciesIDs} | while read line
	do
		#Get species name
		a=$(echo ${line} | cut -d ' ' -f2 | sed s/\.fa//)
		#Get the genes for that species
		genes=$(cut -f ${column} tmp | sed 's/\ //g')
		#Count the genes for that species
		count=$(echo ${genes} | tr ',' '\n' | sed '/^$/d' | wc -l)
		#Output to table
		echo "${a} ${count} ${genes}" | tr ' ' '\t' >> ${i}/${i}_table.tsv
		#Extract the sequences
		samtools faidx $(echo ${genes} | tr ',' '\n' | sed 's/^ *//') >> ${i}/${i}.fa
		#Increase the column number by 1
		column=$(expr ${column} + 1)
	done
	rm tmp
done

echo "Done"


