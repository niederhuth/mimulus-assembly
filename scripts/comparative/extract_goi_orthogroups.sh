#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --job-name extract_goi
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
goi= #path to a csv of genes of interest, if blank, will look for goi.csv in misc dir
speciesIDs= #orthogroup speciesID (SpeciesID.txt), if blank, will look for orthofinder results in current dir
orthogroups= #orthogroups (N0.tsv), if blank, will look for orthofinder results in current dir
datatype="proteins-primary" #Set this to datatype used for orthofinder, default proteins-primary
seqs= #Path to input sequences, if blank, will look for dir named ${datatype} in the current dir
CDS=TRUE #True/FALSE, extract CDS sequences, assumes datatype is proteins-primary or proteins
trim_p=TRUE #TRUE/FALSE, trim the trailing ".p" from the protein name so matches cds name

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/phylo/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/phylo/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
path2=$(pwd)
path3=goi

#Set path to gene of interest list
if [ -z ${goi} ]
then
	goi=${path1}/goi.csv
fi
#set path to the SpeciesID.txt file
if [ -z ${speciesIDs} ]
then
	speciesIDs=$(ls ${path2}/orthofinder/*/WorkingDirectory/SpeciesIDs.txt)
fi
#Set path to N0.tsv orthofinder output
if [ -z ${orthogroups} ]
then
	orthogroups=$(ls ${path2}/orthofinder/*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv)
fi
#Set path for input sequences
if [ -z ${seqs} ]
then
	#seqs="${path2}/${datatype}/${a}.fa"
	path4=$(pwd | sed s/data.*/data/)
fi

#Check for and make/cd working directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#Loop over each goi and get the relevant data
sed '1d' ${goi} | grep -v \# | while read line 
do
	name=$(echo ${line} | cut -d ',' -f1)
	echo "Extracting sequences for ${name}"
	ref_genes=$(echo ${line} | cut -d ',' -f2)
	mkdir ${name} ${name}/seqs
	outseqs="${name}/seqs/${name}"
	#Get the orthogroups
	for x in ${ref_genes}
	do
		grep ${x} ${orthogroups} >> tmp
	done
	#Remove redundant lines
	sort tmp | uniq > tmp2
	mv tmp2 tmp
	#Get the orthogroup(s)
	ogs=$(cut -f1 tmp | sort | uniq | tr '\n' ',' | sed s/\,$//)
	#First species starts on column 4
	column=4
	#Loop over each species
	cat ${speciesIDs} | while read line2
	do
		#Get species name
		a=$(echo ${line2} | cut -d ' ' -f2 | sed s/\.fa//)
		#Get the genes for that species
		genes=$(cut -f ${column} tmp | sed 's/\ //g' | sed s/\,$//)
		proteins=${genes}
		#Trim off the ".p" found at end of some protein sequences
		if [ ${trim_p} = TRUE ]
		then
			genes=$(echo ${genes} | tr ',' '\n' | sed '/^$/d' | sed s/\.p$// | tr '\n' ',' | sed s/\,$//)
		fi
		#Count the genes for that species
		count=$(echo ${genes} | tr ',' '\n' | sed '/^$/d' | wc -l)
		#Output to table
		echo "${a} ${ogs} ${count} $(echo ${genes} | tr ' ' ',')" | tr ' ' '\t' >> ${name}/${name}_table.tsv
		#Extract the sequences
		seqs=$(ls ${path4}/${a/_*/}/${a/*_/}/ref/annotations/${a/*_/}*-${datatype}.fa)
		samtools faidx ${seqs} ${proteins} >> ${outseqs}-${datatype}.fa
		#Extract CDS sequences?
		if [ ${CDS} = TRUE ]
		then
			cds=$(ls ${path4}/${a/_*/}/${a/*_/}/ref/annotations/${a/*_/}*-${datatype/proteins/cds}.fa)
			samtools faidx ${cds} $(echo ${genes} | tr ' ' '\n') >> ${outseqs}-${datatype/proteins/cds}.fa
		fi
		#Increase the column number by 1
		column=$(expr ${column} + 1)
	done
	rm tmp
done

echo "Done"


