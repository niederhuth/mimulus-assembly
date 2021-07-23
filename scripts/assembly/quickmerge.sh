#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200GB
#SBATCH --job-name quickmerge
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
fasta_list="" #If left blank, will search for fasta sequences in submitted directory
HCO=5.0
C=1.5
length_cutoff=0
merging_length_cutoff=5000

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/scaffolding/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/scaffolding/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${genotype}\\/// | sed s/\\/.*//)
path1=$(pwd | sed s/${genotype}.*/${genotype}/)

#If fasta_list not provided, find fastas
if [ -z ${fasta_list} ]
then
	fasta_list=$(ls *fa *fasta *fna | sed s/.*\ // | tr '\n' ' ')
fi

#Run quickmerge on all combinations
for a in ${fasta_list}
do
	for b in ${fasta_list}
	do
		if [ ${a} != ${b} ]
		then
			echo "Mapping ${a} against ${b}"
			mkdir ${a}_vs_${b}
			cd ${a}_vs_${b}
			merge_wrapper.py \
				--prefix ${a}_vs_${b} \
				--hco ${HCO} \
				--c ${C} \
<<<<<<< HEAD
				--length_cutoff ${LENGTH_CUTOFF} \
				--merging_length_cutoff ${MERGING_LENGTH_CUTOFF} \
				../${a} ../${b}
=======
				--length_cutoff ${length_cutoff} \
				--merging_length_cutoff ${merging_length_cutoff} \
				${a} ${b}
>>>>>>> 5e49e91c70989af9dcd2189b80c0e4881e16d7ba
			cd ../
		fi
	done
done

