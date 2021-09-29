#!/bin/bash --login
#SBATCH --time=128:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100GB
#SBATCH --job-name ragtag_correct
#SBATCH --output=%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
ref=$(pwd | sed s/data.*/data/)/Mguttatus/TOL/ref/TOL-v5.fa
threads=10
aligner="unimap"
min_len=1000
merge_dist=100000
min_gap=100
max_gap=100000
infer_gaps="yes" #yes or no
min_orient_score=0.0
min_location_score=0.0
min_grouping_score=0.2

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

if [ ${infer_gaps} ]
then
	gaps="-r -g ${min_gap} -m ${max_gap} "
else
	gaps=
fi

#Run ragtag scaffold
echo "Running ragtag scaffold"
ragtag.py scaffold \
	-f ${min_len} \
	-d ${merge_dist} \
	-i ${min_grouping_score} \
	-a ${min_location_score} \
	-s ${min_orient_score} ${gaps}\
	${ref} ${input}

echo "Done"


