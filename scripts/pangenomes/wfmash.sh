#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20GB
#SBATCH --job-name mafft
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
pidentity=90 #Minimum percent identity

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/pangenome/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/pangenome/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=

#
echo "Running wfmash"
wfmash \
	-p ${pidentity} \
	-n 6 \
	-t 4 \
	-m \
	${species}.fa.gz > ${species}-mapping.paf

echo "Done"