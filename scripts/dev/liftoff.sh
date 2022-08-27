#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500GB
#SBATCH --job-name liftoff
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables


#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/liftoff/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/liftoff/lib:${LD_LIBRARY_PATH}"

#Set temporary directories for large memory operations
export TMPDIR=$(pwd | sed s/data.*/data/)
export TMP=$(pwd | sed s/data.*/data/)
export TEMP=$(pwd | sed s/data.*/data/)

#Run liftoff
echo ""
liftoff \
	-polish \
	-o mapped \
	-u unmapped \
	-g ${gff} \
	-chroms ${chroms} \
	${target} \
	${reference}


gffread -S -y mapped-protein.fa -g S1-v1.fa mapped_polished