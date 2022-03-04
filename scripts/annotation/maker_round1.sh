#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100GB
#SBATCH --job-name maker_round1
#SBATCH --output=%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
mpi=4 #number of threads for mpi
fasta= #input fasta, if left blank, will look for it in current directory
blast_threads=1 #Leave 1 for MPI

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/maker/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/maker/lib:${LD_LIBRARY_PATH}"
#Export path to agusutus config files
#export ZOE="${conda}/envs/maker" #Need to check
export AUGUSTUS_CONFIG_PATH="${conda}/envs/maker/config/"
#export REPEATMASKER_LIB_DIR=
#export REPEATMASKER_MATRICES_DIR=

#Set temporary directories for large memory operations
export TMPDIR=$(pwd | sed s/data.*/data/)
export TMP=$(pwd | sed s/data.*/data/)
export TEMP=$(pwd | sed s/data.*/data/)

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*${species}\\/${genotype}\\/// | sed s/\\/.*//)
path2="maker"

#Look for fasta file, there can only be one!
if ls *.fa >/dev/null 2>&1
then
	fasta=$(ls *fa | sed s/.*\ //)
	echo "Fasta file ${fasta} found"
elif ls *.fasta >/dev/null 2>&1
then
	fasta=$(ls *fasta | sed s/.*\ //)
	echo "Fasta file ${fasta} found"
elif ls *.fna >/dev/null 2>&1
then
	fasta=$(ls *fna | sed s/.*\ //)
	echo "Fasta file ${fasta} found"
else
	echo "No fasta file found, please check and restart"
fi

#Make & cd to directory
if [ -d ${path2} ]
then
	cd ${path2}
else
	mkdir ${path2}
	cd ${path2}
fi

#Run maker
mpiexec -n ${mpi} maker \
	-genome ../${fasta} \
	-cpus ${blast_threads} \
	${path1}/annotation/maker_round1/*

echo "Done"

