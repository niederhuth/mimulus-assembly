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
use_mpi=TRUE #TRUE or FALSE, if FALSE, then run without mpi and use single process
mpi_threads=4 #number of threads for mpi
fasta= #input fasta, if left blank, will look for it in current directory
blast_threads=1 #Leave 1 for MPI

#Load modules from cluster used in setting up maker for mpi
#You will need to modify this according to your own setup
if [ ${use_mpi} = TRUE ]
then
	module load GCC/11.1.0 OpenMPI/4.1.1
fi

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
export TMPDIR=$(pwd)
export TMP=$(pwd)
export TEMP=$(pwd)

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

#Set temporary directories for large memory operations
export TMPDIR=$(pwd)
export TMP=$(pwd)
export TEMP=$(pwd)

#Copy over rmlib
if [ ! -s TElib.fa ]
then
	cp ../edta/*.fa.mod.EDTA.TElib.fa TElib.fa
fi

#Run maker
if [ ${use_mpi} = TRUE ]
then
	#Run maker
	mpiexec -n ${mpi_threads} maker \
		-genome ../${fasta} \
		${path1}/annotation/maker_round1/*
else
	#Run maker without mpi
	maker \
		-genome ../${fasta} \
		${path1}/annotation/maker_round1/*
fi

echo "Done"

