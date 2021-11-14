#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --job-name download_radseq
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#In general dont change this, unless using a similar datatype
#This should match the dataype in the misc/samples.csv file
datatype="radseq"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/scaffolding/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/scaffolding/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
assembly=$(pwd | sed s/^.*\\///)
path2=$(pwd | sed s/${genotype}\\/${sample}\\/.*/${genotype}\\/${sample}/)
path3=${path1}/genetic_map/${datatype}/


#Make data directory
if [ -d fastq/${datatype} ]
then
	cd fastq/${datatype}
else
	mkdir fastq/${datatype}
	cd fastq/${datatype}
fi

#Loop through SRR list and download
for i in ${path3}/*
do
	name=$(echo ${i} | sed s/^.*\\/// | sed s/_SRR.txt//)
	echo "Downloading data for ${name}"
	if [ -d ${name} ]
	then
		cd ${name}
	else
		mkdir ${name}
		cd ${name}
	fi
	cat ${i} | while read line
	do
		if [ ! -f ${line}.gz ]
		then
			if [ ! -f ${line}_2.gz ]
			then
				fastq-dump --split-3 ${line}
				gzip ${line}*
			fi
		fi
	done
	cd ../
done	

echo "Done"

