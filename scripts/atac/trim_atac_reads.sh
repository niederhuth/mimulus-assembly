#!/bin/bash --login
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20GB
#SBATCH --job-name trim_atac_reads
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set variables
threads=10
min_length=30 #If left blank, will use defaul length of 15 nucleotides
adapter_fasta= #"${conda}/envs/atac/share/trimmomatic/adapters" #Path to adapter fasta if set, will use instead of automatic detection

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"
#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/atac/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/atac/lib:${LD_LIBRARY_PATH}"

#Other variables, these should not have to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*data\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*data\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)

#Change directories to data
cd fastq/atac

#Names of output fastq files, these should not have to be changed
r1="combined.1.fastq.gz"
r2="combined.2.fastq.gz"
t1="fastp_trimmed.1.fastq.gz"
t2="fastp_trimmed.2.fastq.gz"
t3="trimmed.1.single.fastq.gz"
t4="trimmed.2.single.fastq.gz"

#Check if data is single-end or paired-end and combine reads from multiple files
if ls *_R1_001.fastq.gz >/dev/null 2>&1
then
	if ls *_R2_001.fastq.gz >/dev/null 2>&1
	then
		echo "Data is Paired-end"
		#Set PE to TRUE to run trimmomatic in paired-end mode
		PE="TRUE"
		#Combine the reads for trimming
		cat *_R1_001.fastq.gz > $r1
		cat *_R2_001.fastq.gz > $r2
	else
		echo "Data is Single-end"
		#Set PE to FALSE to run trimmomatic in single-end mode
		PE="FALSE"
		#Combine the reads for trimming
		cat *_R1_001.fastq.gz > $r1
	fi
elif ls *_1.fastq.gz >/dev/null 2>&1
then
	if ls *_2.fastq.gz >/dev/null 2>&1	
	then
		echo "Data is Paired-end"
		#Set PE to TRUE to run trimmomatic in paired-end mode
		PE="TRUE"
		#Combine the reads for trimming
		cat *_1.fastq.gz > $r1
		cat *_2.fastq.gz > $r2
	else
		echo "Data is Single-end"
		#Set PE to FALSE to run trimmomatic in single-end mode
		PE="FALSE"
		#Combine the reads for trimming
		cat *_1.fastq.gz > $r1
	fi
else
	echo "Data Missing"
fi

#Set parameters
parameters="--thread ${threads}"
if [[ ! -z ${min_length} ]]
then
	parameters="${parameters} --length_required ${min_length}"
fi
if [[ ! -z ${adapter_path} ]]
then
	parameters="${parameters} --adapter_fasta ${adapter_fasta}"
fi

#Trim & QC reads
if [ ${PE} = "TRUE" ]
then
	echo "Running fastp"
	#Run in paired-end mode
	fastp \
		${parameters} \
		--html ${sample}.html \
		--json ${sample}.json \
		--in1 ${r1} \
		--in2 ${r2} \
		--out1 ${t1} \
		--out2 ${t2} \
		--unpaired1 ${t3} \
		--unpaired2 ${t4}
	echo "Running fastqc"
	mkdir trimmed_fastqc
	fastqc -t ${threads} -o fastp_fastqc/ ${t1} ${t2} ${t3} ${t4}
elif [ ${PE} = "FALSE" ]
then
	echo "Running fastp"
	#Run in single-end mode
		fastp \
		${parameters} \
		--html ${sample}.html \
		--json ${sample}.json \
		--in1 ${r1} \
		--out1 ${t1} 
	#Run fastqc on trimmed data
	echo "Running fastqc"
	mkdir trimmed_fastqc
	fastqc -t ${threads} -o fastp_fastqc/ ${t1}
fi

echo "Done"
