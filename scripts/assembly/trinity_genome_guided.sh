#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=100GB
#SBATCH --job-name trinity-genome-guided
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=40
max_memory=
intron=

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/transcript-assembly/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/transcript-assembly/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/^.*\\///)
path2="trinity_gg"




#Run Trinity
Trinity \
	--seqType fq \
	--max_memory ${max_memory} \
	--left ${r1} \
	--right ${r2} \
	--SS_lib_typ RF \
	--CPU ${threads} \
	--genome_guided_bam ${bam} \
	--genome_guided_max_intron ${intron} \
	--trimmomatic \
	--output ${path2}