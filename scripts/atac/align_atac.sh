#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=60GB
#SBATCH --job-name align_atac
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set variables
threads=10 #threads for bowtie2
sort_threads=4 #threads for sorting bam file
max_insert=1000
secondary_alignments=10 #Sets -k parameter to keep secondary alignments, necessary for genrich

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
datatype="atac"

#Set various options for bowtie2
bowtie2_options="-p ${threads} --very-sensitive"
if [[ ! -z ${max_insert} ]]
then
	bowtie2_options="${bowtie2_options} -X 1000"
fi
if [[ ! -z ${secondary_alignments} ]]
then
	bowtie2_options="${bowtie2_options} -k ${secondary_alignments}"
fi

#Set fastq files
t1="$(pwd)/fastq/atac/trimmed.1.fastq.gz"
t2="$(pwd)/fastq/atac/trimmed.2.fastq.gz"
#Check if paired-end or single-end
if [ -f ${t2} ]
then
	echo "Data is paired-end"
	PE="TRUE"
else
	echo "Data is single-end"
	PE="FALSE"
fi

#Creat output directory & cd to it
if [ -d ${datatype} ]
then 
	cd ${datatype}
else
	mkdir ${datatype}
	cd ${datatype}
fi

#Get list of genomes
genomes=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $5 == d) print $7}' \
	${path1}/samples.csv)

#Loop over each reference genome and align to it
for i in ${genomes}
do
	#Set path to bowtie2 index
	path2=$(pwd | sed s/data\\/${species}\\/.*/data\\/${species}\\/${i}/)
	version=$(ls ${path2}/ref/${i}-v*.fa | sed s/.*\-v// | sed s/.fa//)
	bowtie2_index="${path2}/ref/bowtie2/${i}-v${version}"
	#Set bam name
	bam="${sample}_ref_${i}-v${version}.bam"

	#Check if the bam index already exists, if it does, then indicates data was already successfully aligned
	if [ -f ${bam}.bai ]
	then
		echo "Existing bam file found, skipping" to mark duplicates""
		echo "To rerun this step, delete ${bam} & ${bam}.bai and resubmit"
	elif [[ ! -f ${bam}.bai && ${PE} = "TRUE" ]]
	then
		#Align paired-end data
		echo "Aligning to ${i}-v${version}"
		bowtie2 \
			${bowtie2_options} \
			-x ${bowtie2_index} \
			-1 ${t1} \
			-2 ${t2} | samtools view -@ ${sort_threads} -bSh | samtools sort -@ ${sort_threads} > ${bam}
		#Index the bam file
		echo "Indexing ${bam}"	
		samtools index ${bam}
	elif [[ ! -f ${bam}.bai && ${PE} = "FALSE" ]]
	then
		#Align single-end data
		echo "Aligning to ${i}-v${version}"
		bowtie2 \
			${bowtie2_options} \
			-x ${bowtie2_index} \
			-U ${t1} | samtools view -@ ${sort_threads} -bSh | samtools sort -@ ${sort_threads} > ${bam}
		#Index the bam file
		echo "Indexing ${bam}"	
		samtools index ${bam}
	fi
	echo "Alignment to ${i}-v${version} complete"
done

echo "Done"

