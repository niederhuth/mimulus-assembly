#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=800GB
#SBATCH --job-name tombo-annotate-raw-with-fastqs
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20
fast5_dir="fast5_pass/"
fastq="pass.fastq.gz"
sequencing_summary="sequencing_summary_PAG35136_e70c34ec.txt"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/deepsignal/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/deepsignal/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/data\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/data\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
datatype="ont"
path2=

#Change to fastq file
cd fastq/ont

#Run tombo annotate_raw_with_fastqs
echo ""
tombo preprocess annotate_raw_with_fastqs \
	--processes ${threads} \
	--fast5-basedir ${fast5_dir} \
	--fastq-filename ${fastq} \
	--sequencing-summary-filenames ${sequencing_summary} \
	--basecall-group Basecall_1D_000 \
	--basecall-subgroup BaseCalled_template \
	--overwrite

echo "Done"