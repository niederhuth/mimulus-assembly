#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20GB
#SBATCH --job-name partition-before-pggb
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20
pidentity=90 #Minimum percent identity

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/pangenome/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/pangenome/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=

#Change sequence names to PanSN-spec: https://github.com/pangenome/PanSN-spec
echo "Renaming Sequences"
for seq in ${genome_list}
do
	echo "${seq} ... Starting"
	sed -i s/\>/\>${genotype}\#${haplotype}\#/ ${seq}
	echo "${seq} ... Done"
done

#Combine the fasta files
echo "Combining input fastas"
cat ${genome_list} | bgzip --threads ${threads} > ${species}.fa.gz
#Index the combined fasta
samtools faidx ${species}.fa.gz

#
echo "Running partition-before-pggb"
partition-before-pggb \
	-i ${species}.fa.gz \
	-o partition \
	-




echo "Done"