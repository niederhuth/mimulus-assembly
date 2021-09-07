#!/bin/bash --login
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --job-name soft_mask
#SBATCH --output=../job_reports/%x-%j.SLURMout

#This script is for a quick (relatively) soft-masking of genomes for whole-genome alignment
#If you have the time, I'd recommend running the repeatmasker.sh in the scripts/assembly directory instead.


#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
engine="rmblast" #crossmatch, wublast, abblast, ncbi, rmblast, hmmer
threads=5 #Actual cores used by RM are thread number multiplied by the cores for search engine used.
			#These are: RMBlast=4 cores, ABBlast=4 cores, nhmmer=2 cores, crossmatch=1 core
			#So 10 threads with RMBlast actually needs 40 cores!

#For better masking, I recommend using EDTA or another program to generate a set of sequences.
#However, if you are going to do that, why not use scripts/assembly/repeatmasker.sh for a proper job?
repeats="${conda}/envs/EDTA/share/EDTA/database/maizeTE11122019" 

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/EDTA/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/EDTA/lib:${LD_LIBRARY_PATH}"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${genotype}\\/// | sed s/\\/.*//)
path2="repeatmasker"

#Look for fasta file, there can only be one!
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

#mask with repeatmasker
if ls ${path2} >/dev/null 2>&1
then
	echo "Previous RepeatMasker results found."
	echo "Quitting. To repeat this analysis, delete ${path2} and resubmit."
else
	mkdir ${path2}
	RepeatMasker \
		-e ${engine} \
		-pa ${threads} \
		-qq \
		-no_is \
		-norna \
		-xsmall \
		-dir ${path2} \
		-lib ${repeats} \
		${input}
fi

echo "Done"
