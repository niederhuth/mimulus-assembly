#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=200GB
#SBATCH --job-name tombo-resquiggle
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=50
fast5_dir="fast5"
ref_fasta= #reference genome, if left blank, will look for in the ref directory for that genotype

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/deepsignal/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/deepsignal/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
path2=$(pwd | sed s/data.*/data/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/data\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/data\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
datatype="ont"
path3="methylC_ont"

#Check for and make/cd working directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#find the reference fasta if not provided
if [ -z ${ref_fasta} ]
then
	#Get version number for genome fasta
	genome_ver=$(ls ${path2}/${species}/${genotype}/ref/${genotype}-v*.fa | head -1 | \
		sed s/.*\-v// | sed s/\.fa//)
	#Set the ref_fasta path
	ref_fasta="${path2}/${species}/${genotype}/ref/${genotype}-v${genome_ver}.fa"
fi 

#Run tombo annotate_raw_with_fastqs
echo "Running tombo resquiggle"
tombo resquiggle \
	${fast5_dir}/ \
	${ref_fasta} \
	--processes ${threads} \
	--corrected-group RawGenomeCorrected_000 \
	--basecall-group Basecall_1D_000 \
	--overwrite

echo "Done"
