#!/bin/bash --login
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=44
#SBATCH --mem=80GB
#SBATCH --job-name macs2
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20

#In general dont change this, unless using a similar datatype, e.g. DNase-seq
#This should match the dataype in the misc/samples.csv file
datatype="atac"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/atac/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/atac/lib:${LD_LIBRARY_PATH}"
#Path to picard
picard="${conda}/envs/atac/share/picard-2.23.9-0/picard.jar"
#Path to trimmomatic fastas 
adapter_path="${conda}/envs/atac/share/trimmomatic/adapters"

#Other variables, these should not have to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/^.*\\///)
path2=${datatype}

#Set filetype
#For now we assume data will be in BAM format
path3="fastq/${datatype}"
if ls ${path3}/*_R1_001.fastq.gz >/dev/null 2>&1
then
	if ls ${path3}/*_R2_001.fastq.gz >/dev/null 2>&1
	then
		echo "Data is Paired-end"
		filetype="BAMPE"
	else
		echo "Data is Single-end"
		echo "Single-end ${datatype}? If this is wrong, double check and restart"
		filetype="BAM"
	fi
elif ls ${path3}/*_1.fastq.gz >/dev/null 2>&1
then
	if ls ${path3}/*_2.fastq.gz >/dev/null 2>&1	
	then
		echo "Data is Paired-end"
		filetype="BAMPE"
	else
		echo "Data is Single-end"
		echo "Single-end ${datatype}? If this is wrong, double check and restart"
		filetype="BAM"
	fi
else
	echo "Data Missing"
fi

#Get genomes
genome=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $4 == d) print $7}' \
	${path1}/samples.csv)

#Get mappable genome size
genomeSize=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $4 == d) print $9}' \
	${path1}/samples.csv)

#MACS2
a=1
for i in ${genomes}
do
	if <something>
	then
		echo "Something"
		gs=$(cut -d ' ' -f ${a} ${genomeSize})
		echo "Effective genome size is ${gs}"
		macs2 callpeak \
			-T ${input} \
			-C ${control} \
			-N ${sample} \
			--outdir <something> \
			-F ${filetype} \
			-G ${gs} \
			-Q 0.05 
	fi
	a=$(expr $a + 1)
done

echo "Done"

