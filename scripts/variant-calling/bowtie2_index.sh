#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50GB
#SBATCH --job-name bowtie2-index
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/variant-calling/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/variant-calling/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*data\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*data\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
version=$(ls ${genotype}*-v*.fa | sed s/.*\-v// | sed s/\.fa//)

#Make directory and cd to it
if [ -d bowtie2 ]
then
	cd bowtie2
else
	mkdir bowtie2
	cd bowtie2
fi

#Index the genome
echo "Building bowtie2 index for ${species} ${genotype} version ${version}"
bowtie2-build ../${genotype}-v${version}.fa ${genotype}-v${version}

echo "Done"
