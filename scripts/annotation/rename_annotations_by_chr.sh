#!/bin/bash

#Set variables from command line
while getopts p:j:z:


prefix=$1
justify=$2
zeros_at_end=$3
input_gff=
input_proteins=
input_transcripts=

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Export paths to conda
export PATH="${conda}/envs/maker/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/maker/lib:${LD_LIBRARY_PATH}"

#Get chromosome list
cut -f1 ${input_gff} | grep -v \# | sort | uniq > chr_list

#Loop over chr_list and rename genes for each chromosome
mkdir tmp
cd tmp
cat ../chr_list | while read line
do
	echo "Working on ${line}"
	awk -v a=${line} '$1==a' ../${input_gff} > ${line}.gff
	if [[ ${line:0:3} = "chr" ]]
	then
		if [ $(echo ${line} | wc -c) -gt 5 ]
		then
			chr=$(echo ${line} | sed s/chr//)
		else
			chr=$(echo ${line} | sed s/chr/0/)
		fi
	else
		chr="UN"
	fi
	maker_map_ids \
		--prefix ${prefix}${chr}g \
		--justify ${justify} \
		--iterate 1 \
		${line}.gff | awk -v a=${zeros_at_end} '{if ($2 ~ /-R/) print $0; else print $0a}' | sed s/-R/${zeros_at_end}./ > ${line}_renamed.map
	chr=
done

#Combine files
cd ..
cat tmp/${line}_renamed.map > ${species}_${genotype}_${sample}_renamed_genes.map
#rm -R tmp

#Rename gff & fasta files
map_gff_ids ${species}_${genotype}_${sample}_renamed_genes.map ${input_gff} > ${output}_renamed.gff
map_fasta_ids ${species}_${genotype}_${sample}_renamed_genes.map ${input_proteins} > ${output}-transcripts_renamed.fa
map_fasta_ids ${species}_${genotype}_${sample}_renamed_genes.map ${input_proteins} > ${output}-proteins_renamed.fa

echo "Done"
