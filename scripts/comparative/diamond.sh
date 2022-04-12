#!/bin/bash --login
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=50GB
#SBATCH --job-name diamond
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
targets=$()
threads=50
evalue=
max_target_seqs=5

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/EDTA/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/EDTA/lib:${LD_LIBRARY_PATH}"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
path2=$(pwd | sed s/data.*/data/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${genotype}\\/// | sed s/\\/.*//)
path3="diamond"

#Make output directory and change to it
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#Run diamond
for i in ${targets}
do
	echo "Checking for ${i} diamond database"
	if [ -f ${path2}/${}/${}/ref/annotations/${}.dmnd ]
	then
		echo "Found ${i} diamond database"
	else
		echo "Makind ${i} diamond database"
		diamond makedb \
			--in ${path2}/${}/${}/ref/annotations/${}-proteins.fa \
			--db ${path2}/${}/${}/ref/annotations/${}.dmnd
	fi
	
	if [ -f ${x}-${i}.m8 ]
	then
		echo "${x}-${i}.m8 already exists, skipping"
		echo "To rerun diamond blastp on ${x}-${i}.m8, then delete ${x}-${i}.m8 and resubmit"
	else
		echo "Running ${species} ${i} diamond blastp"
		diamond blastp \
			--threads ${threads} \
			--db ${path2}/${}/${}/ref/annotations/${}.dmnd \
			--query ${x}-protein.fa \
			--out ${x}-${i}.m8 \
			--un ${x}-${i}-un.fa \
 			--more-sensitive \
			--evalue ${evalue} \
			--max-target-seqs ${max_target_seqs} \
			--unal 0
	fi
done

