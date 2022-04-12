#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200GB
#SBATCH --job-name=lastz
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
fa1=scaffold_507.fa
fa2=contig_1056.fa

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/mummer/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/mummer/lib:${LD_LIBRARY_PATH}"

#Align Sequences
lastz \
	${fa1} \
	${fa2} \
	--strand=both \
	--inner=1000 \
	--format=maf \
	--output=${fa2}_${fa1}_lastz

echo "Done"
