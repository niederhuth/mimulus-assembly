#!/bin/bash --login
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --job-name bam2pileup
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
datatype="poolseq"
ref_fasta="$(pwd | sed s/Mguttatus.*/Mguttatus/)/L1/ref/L1-v1.fa"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/pangenome/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/pangenome/lib:$LD_LIBRARY_PATH"

#Other variables, these should not have to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*data\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*data\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)

#Set input
input="giraffe/${sample}_${datatype}_sorted.bam"
#Set output
output="giraffe/${sample}_${datatype}.pileup"

#Set arguments
arguments="${input}"
#Add reference fasta
if [[ ! -z ${ref_fasta} ]]
then
	arguments="-f ${ref_fasta} ${arguments}"
fi

# Compute the read support from the gam
echo "Converting to pileup"
samtools mpileup ${arguments} > ${output}

echo "Done"