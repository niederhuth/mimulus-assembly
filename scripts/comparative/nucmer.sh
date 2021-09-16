#!/bin/bash --login
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=200GB
#SBATCH --job-name nucmer
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=10
fa1=S1.fa
fa2=L1.fa

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/mummer/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/mummer/lib:${LD_LIBRARY_PATH}"

#
a=$(echo ${fa1} | sed s/.fa//)
b=$(echo ${fa2} | sed s/.fa//)

#Align Genomes
echo "Running nucmer ${fa2} against ${fa1}"
mkdir ${b}_${a}
nucmer \
	-t ${threads} \
	-p ${b}_${a}/${b}_${a} \
	${fa1} ${fa2}

echo "Running nucmer ${fa2} against ${fa1}"
mkdir ${a}_${b}
nucmer \
        -t ${threads} \
        -p ${a}_${b}/${a}_${b} \
        ${fa2} ${fa1}

echo "Done"
