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
fast5_dir="$(pwd)/fastq/ont/fast5_pass/"
fastq="$(pwd)/fastq/ont/pass.fastq"
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

#If data is in multi-fast5 format convert to single-fast5 format
if [ ${convert_to_single} = TRUE ]
then
	echo "Converting multi-fast5 to single-fast5"
	multi_to_single_fast5 \
		-t ${threads} \
		-i ${path2}/${fast5_dir} \
		-s fast5/
	fast5_dir="$(pwd)/fast5"
fi

#Unzip the fastq data if it is compressed
if [[ $(pwd)/${fastq} =~ \.gz$ ]]
then
	echo "File is gzip file, uncompressing fastq"
	gunzip -c ${fastq} > combined.fastq
	fastq="$(pwd)/combined.fastq"
elif [[ ${fastq} =~ \.fastq$ ]]
	echo "File is uncompressed"
fi

#Run tombo annotate_raw_with_fastqs
echo "Running tombo preprocess annotate_raw_with_fastqs"
tombo preprocess annotate_raw_with_fastqs \
	--processes ${threads} \
	--fast5-basedir ${fast5_dir} \
	--fastq-filename ${fastq} \
	--sequencing-summary-filenames ${sequencing_summary} \
	--basecall-group Basecall_1D_000 \
	--basecall-subgroup BaseCalled_template \
	--overwrite

echo "Done"
