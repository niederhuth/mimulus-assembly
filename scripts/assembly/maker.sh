#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=100GB
#SBATCH --job-name maker
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=40

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/maker/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/maker/lib:$LD_LIBRARY_PATH"
#Export path to agusutus config files
export ZOE="${conda}/envs/maker" #Need tp check
export AUGUSTUS_CONFIG_PATH="${conda}/envs/maker/config/"
export REPEATMASKER_LIB_DIR=
export REPEATMASKER_MATRICES_DIR=

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
CTL_FILES="${path1}/maker_ctl/*"
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/^.*\\///)
path2="maker"

#Run Maker
maker \
	-base ${path2} \
	-cpus ${threads} \
	${CTL_FILES}

