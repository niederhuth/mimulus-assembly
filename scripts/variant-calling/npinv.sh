#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200GB
#SBATCH --job-name medaka
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
region="all" #Specify a region for analysis, e.g. chr9:1-1000 OR chr9 OR all
minAln=500 #minimum size for Alignment & Inv Default[500]
min=500 #minimum size of an inversion Default[500]
max=5000000 #maximum size of an inversion.  Default[10000]
window=2000 #minimun window size (bp) to merge inversion breakpoints. Default[2000]
threshold=3 #minimum number of supporting reads for an inversion. Default[3]

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/polishing/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/polishing/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${genotype}\\/// | sed s/\\/.*//)

#Run npInv
for i in *bam
do
	vcf=$(echo ${i} | sed s/\.bam/-npinv.vcf/)
	if [ -s ${vcf} ]
	then
		echo "${vcf} file found."
		echo "To rerun this analysis delete ${vcf} and resubmit"
	else
		echo "Beginning npInv on ${i}"
		npinv \
			--input ${i} \
			--output ${vcf} \
			--region ${region} \
			--minAln ${minAln} \
			--min ${min} \
			--max ${max} \
			--window ${window} \
			--threshold ${threshold}
	fi
done

echo "Done"


