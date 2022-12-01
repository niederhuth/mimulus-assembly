#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=400GB
#SBATCH --job-name tombo-annotate-raw-with-fastqs
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=40
fast5_dir="fast5"
fastq="combined.fastq"
sequencing_summary="$(pwd)/fastq/ont/sequencing_summary_PAG35136_e70c34ec.txt"
convert_to_single=TRUE #TRUE or FALSE, fast5 is multi-fast5 and should be converted to single-fast5

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
path2="methylC_ont"

#Check for and make/cd working directory
if [ -d ${path2} ]
then
	cd ${path2}
else
	mkdir ${path2}
	cd ${path2}
fi

#Remove reads that did not pass filtering from the sequencing summary
awk '$10 != "FALSE"' ${sequencing_summary} > passed_filter_sequencing_summary.txt

#Run tombo annotate_raw_with_fastqs
echo "Running tombo preprocess annotate_raw_with_fastqs"
tombo preprocess annotate_raw_with_fastqs \
	--processes ${threads} \
	--fast5-basedir ${fast5_dir} \
	--fastq-filename ${fastq} \
	--sequencing-summary-filenames passed_filter_sequencing_summary.txt \
	--basecall-group Basecall_1D_000 \
	--basecall-subgroup BaseCalled_template \
	--overwrite

echo "Done"
