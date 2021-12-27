#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=54
#SBATCH --mem=80GB
#SBATCH --job-name align_radseq2
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=50
threads2=4
fasta="flye/medaka_racon_3/tigmint_long_pilon_2/pilon_2_tigmint.fa"

#In general dont change this, unless using a similar datatype
#This should match the dataype in the misc/samples.csv file
datatype="radseq"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/scaffolding/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/scaffolding/lib:${LD_LIBRARY_PATH}"
#Path to picard
picard="${conda}/envs/variant-calling/share/picard-*/picard.jar"
#Path to trimmomatic fastas 
adapter_path="${conda}/envs/variant-calling/share/trimmomatic/adapters"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
assembly=$(pwd | sed s/^.*\\///)
path2=$(pwd | sed s/${genotype}\\/${sample}\\/.*/${genotype}\\/${sample}/)
path3=${path1}/genetic_map/${datatype}
ref=$(echo ${fasta} | sed s/.*\\///)

#Make and cd to analysis directory
if [ -d ${datatype}_2 ]
then
	cd ${datatype}_2
else
	mkdir ${datatype}_2
	cd ${datatype}_2
fi

#Check for and copy fasta file
if [ ! -f ${ref} ]
then
	cp ${path2}/${fasta} ${ref}
fi

#Make index
if [ -s ${ref}.sa ]
then
	echo "BWA index found"
else
	echo "Indexing fasta"
	bwa index ${ref}
fi

#Loop over each dataset and align
for i in ${path3}/*
do
	name=$(echo ${i} | sed s/^.*\\/// | sed s/_SRR.txt//)
	echo "Working on dataset ${i}"
	if [ -d ${name} ]
	then
		cd ${name}
	else
		mkdir ${name}
		cd ${name}
	fi
	cat ${i} | while read line
	do
		if [ ! -s ${line}.bam.bai ]
		then
			echo "Aligning ${line}"
			r1=${path2}/fastq/${datatype}/${name}/${line}.fastq.gz
			#Define Read Group
			if [ $(echo ${line} | cut -c 1,2,3) = SRR ]
			then
				ID=$(zcat ${r1} | head -1 | cut -d ' ' -f1 | sed s/\@//)
				PU=$(zcat ${r1} | head -1 | cut -d ' ' -f1 | sed s/\@//)
				SM=$(pwd | sed s/^.*\\///)
				PL="ILLUMINA"
				LB="lib1"
			else
				ID=$(zcat ${r1} | head -1 | cut -d ':' -f 3,4 | tr ':' '.')
				PU=$(zcat ${r1} | head -1 | cut -d ':' -f 3,4,10 | tr ':' '.')
				SM=$(pwd | sed s/^.*\\///)
				PL="ILLUMINA"
				LB="lib1"
			fi
			#Align data
			bwa mem -t ${threads} \
				-R "@RG\tID:${ID}\tLB:${LB}\tPL:${PL}\tSM:${SM}\tPU:${PU}" \
				-M ../${ref} ${r1} | \
				samtools view -@ ${threads2} -bSh | \
				samtools sort -@ ${threads2} > ${line}.bam
			#Index data
			samtools index ${line}.bam
		else
			echo "${line} already aligned, skipping to next sample"
		fi
	done
	cd ../
	echo "Dataset ${i} complete"
done

echo "Done"
