#!/bin/bash --login
#SBATCH --time=128:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=500GB
#SBATCH --job-name ragtag_correct
#SBATCH --output=../../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=5
datatype="ont"
aligner="unimap"
min_len=1000
merge_dist=100000
break_dist=5000
cov_win_size=8000
min_cov=
max_cov=
ref=$(pwd | sed s/data.*/data/)/Mguttatus/IM62/ref/IM62-v2.fa
reads=clean.fastq.gz #Set blank to run without read validation
#reads=

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/test/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/test/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
assembly=$(pwd | sed s/^.*\\///)
path2=$(pwd | sed s/${genotype}\\/${sample}.*/${genotype}\\/${sample}/)

#Check for coverage limits
if [ -z ${min_cov} ]
then
	if [ -z ${max_cov} ]
	then
		cov_vars="-v ${cov_win_size}"
	else
		cov_vars="-v ${cov_win_size} --max-cov ${max_cov}"
	fi
else
	if [ -z ${max_cov} ]
	then
		cov_vars="-v ${cov_win_size} --min-cov ${min_cov}"
	else
		cov_vars="-v ${cov_win_size} --min-cov ${min_cov} --max-cov ${max_cov}"
	fi
fi

#Look for fasta file, there can only be one!
if [ -z ${input} ]
then
	echo "No input fasta provided, looking for fasta"
	if ls *.fa >/dev/null 2>&1
	then
		input=$(ls *fa | sed s/.*\ //)
		echo "Fasta file ${input} found"
	elif ls *.fasta >/dev/null 2>&1
	then
		input=$(ls *fasta | sed s/.*\ //)
		echo "Fasta file ${input} found"
	elif ls *.fna >/dev/null 2>&1
	then
		input=$(ls *fna | sed s/.*\ //)
		echo "Fasta file ${input} found"
	else
		echo "No fasta file found, please check and restart"
	fi
else
	echo "Input fasta: ${input}"
fi

#Run ragtag correct
#Check if reads are provided
if [ -z ${reads} ]
then
	#If not, then run based on just alignment to the genome
	echo "No reads provided, running without validation"
	echo "Running ragtag correct on ${assembly} against ${ref}"
		ragtag.py correct \
		-t ${threads} \
		--aligner ${aligner} \
		-f ${min_len} \
		--remove-small \
		-d ${merge_dist} \
		-b ${break_dist} \
		-u \
		${ref} ${input}
else
	#If reads provided, run with read validation
	echo "Using ${reads} for validation"
	echo "Running ragtag correct on ${assembly} against ${ref}"
	ragtag.py correct \
		-t ${threads} \
		--aligner ${aligner} \
		-f ${min_len} \
		--remove-small \
		-d ${merge_dist} \
		-b ${break_dist} \
		-u \
		--read-aligner minimap2 \
		-R ${path2}/fastq/${datatype}/${reads} \
		-T ${datatype} \
		${cov_vars} \
		${ref} ${input}
fi

echo "Done"


