#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --job-name effective_genome_size
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
kmers="35 50 75 100 150"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/atac/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/atac/lib:${LD_LIBRARY_PATH}"

#Other variables, these should not have to be changed, but should set automatically
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
version=$(ls ref/${genotype}-v*.fa | sed s/.*\-v// | sed s/.fa//) 
fasta="ref/${genotype}-v${version}.fa"

#Run unique-kmers.py
for i in ${kmers}
do
	unique-kmers.py -q -R tmp -k ${i} ${fasta}
	head -1 tmp | tr ' ' '\t' | cut -f2,1 >> ref/${genotype}-v${version}_effective_genome_sizes.txt
	rm tmp
done

echo "Done"
