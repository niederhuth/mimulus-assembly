#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=40GB
#SBATCH --job-name format-ont-for-tombo-mods
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=40
fast5_dir="$(pwd)/fastq/ont/fast5_pass/"
fastq="$(pwd)/fastq/ont/combined.fastq.gz"
sequencing_summary="$(pwd)/fastq/ont/sequencing_summary_PAG35136_e70c34ec.txt"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/deepsignal/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/deepsignal/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/data\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/data\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
datatype="ont"
path2="ont_mods"

#Check for and make/cd working directory
if [ -d ${path2} ]
then
	cd ${path2}
else
	mkdir ${path2}
	cd ${path2}
fi

#Convert multi-fast5 to single-fast5
echo "Converting multi-fast5 to single-fast5"
multi_to_single_fast5 \
	-t ${threads} \
	-i ${fast5_dir} \
	-s fast5

#Unzip the fastq data if it is compressed
if [[ $(pwd)/${fastq} =~ \.gz$ ]]
then
	if [ -f combined.fastq ]
	then
		echo "Uncompressed fastq file found"
	else
		echo "File is gzip file, uncompressing fastq"
		gunzip -c ${fastq} > combined.fastq
	fi
elif [[ ${fastq} =~ \.fastq$ ]]
then
	echo "File is uncompressed, copying over"
	cp ${fastq} combined.fastq
fi

#Remove reads that did not pass filtering from the sequencing summary
if [[ ! -f passed_filter_sequencing_summary.txt ]]
then
	awk -v OFS="\t" '{if ($1 == "filename_fastq") print $0; 
		else if ($10 != "FALSE") print $1,$3".fast5",$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' \
		${sequencing_summary} > passed_filter_sequencing_summary.txt
fi

echo "Done"
