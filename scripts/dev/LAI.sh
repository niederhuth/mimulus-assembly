#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=400GB
#SBATCH --job-name edta
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=40

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/EDTA/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/EDTA/lib:${LD_LIBRARY_PATH}"

#Set temporary directories for large memory operations
export TMPDIR=$(pwd | sed s/data.*/data/)
export TMP=$(pwd | sed s/data.*/data/)
export TEMP=$(pwd | sed s/data.*/data/)

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*/${species}\\/${genotype}\\/// | sed s/\\/.*//)
version=$(ls ${sample}-v*.fa | sed s/.*\-v// | sed s/.fa//) 
path2="LAI"

#Make output directory

#Run LTR_retriever using EDTA results
~/Dev/LTR_retriever/LTR_retriever \
	-genome ../${sample}-v${version}.fa \
	-inharvest ../edta/S1-v1.fa.mod.EDTA.raw/LTR/S1-v1.fa.mod.rawLTR.scn \
	-threads ${threads}

#Run LAI
~/Dev/LTR_retriever/LAI \
	-genome ${sample}-v${version}.fa \
	-intact ${sample}-v${version}.fa.pass.list \
	-all ${sample}-v${version}.fa.out

~/Dev/LTR_retriever/LAI -genome S1-v1.fa.mod -intact S1-v1.fa.mod.pass.list -all S1-v1.fa.mod.out