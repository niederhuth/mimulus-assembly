#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --job-name genotype_giraffe
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20
index="$(pwd | sed s/Mguttatus.*/Mguttatus/)/pangenome/giraffe/index.giraffe.gbz"
qual_cutoff=5 #-Q X ignore mapping and base qualitiy < X
ignore_bp=5 #-s X ignore first and last X bp from each read
datatype="poolseq"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/pangenome/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/pangenome/lib:$LD_LIBRARY_PATH"

#Other variables, these should not have to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*data\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*data\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)

#Set output
output="giraffe/${sample}_${datatype}"

# Compute the read support from the gam
echo "Computing read support"
vg pack \
	-t ${threads} \
	-x ${index} \
	-g giraffe/${sample}_${datatype}.gam \
	-Q ${qual_cutoff} \
	-s ${ignore_bp} \
	-o ${output}.pack

# Generate a VCF from the support
echo "Calling variants"
vg call \
	-t ${threads} \
	-a \
	${index} \
	-k ${output}.pack > ${output}.vcf

echo "Done"