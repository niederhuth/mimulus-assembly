#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200GB
#SBATCH --job-name purge_dups
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20
fasta_list=""

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/quickmerge/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/quickmerge/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${genotype}\\/// | sed s/\\/.*//)
path1=$(pwd | sed s/${genotype}.*/${genotype}/)
path2="merge"

#Make output directory
if [ -d ${path2} ]
then
	cd ${path2}
else
	mkdir ${path2}
	cd ${path2}
fi

#Run nucmer
for a in fasta_list
do
	for b in fasta_list
	nucmer \
		-l 100 \
		-p ${out} \
		${fasta1} \
		${fasta2}

#Delta-filter
delta-filter -r -q -l 10000 out.delta > out.rq.delta

#Quickmerge
quickmerge -d out.rq.delta -q hybrid_assembly.fasta -r self_assembly.fasta -hco 5.0 -c 1.5 -l n -ml m -p prefix


echo "Done"