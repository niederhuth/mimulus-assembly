#!/bin/bash --login
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --job-name align-genomes
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20
fa1=../S1/ref/S1-vflye.fa
fa2=../:L1/ref/L1-vflye.fa

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/scaffolding/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/scaffolding/lib:${LD_LIBRARY_PATH}"

#Align Genomes
echo "Running nucmer"
nucmer \
	--threads ${threads} \
	--prefix ${output} \
	${fa1} ${fa2}

echo "Done"
