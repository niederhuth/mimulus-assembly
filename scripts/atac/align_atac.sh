#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=60GB
#SBATCH --job-name align_atac
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set variables
threads=20 #threads for bowtie2
sort_threads=4 #threads for sorting bam file
compression_threads=4 #threads for compressing bam file
max_insert=1000
secondary_alignments= #Sets -k parameter to keep secondary alignments, necessary for genrich
sort_by_read_names=TRUE
dovetail=FALSE
local_alignment=FALSE

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
bowtie2_options="-p ${threads}"
#
if [ ${local_alignment} = "TRUE" ]
then
	bowtie2_options="${bowtie2_options} --very-sensitive-local"
else
	bowtie2_options="${bowtie2_options} --very-sensitive"
fi
#Set max insert size
if [[ ! -z ${max_insert} ]]
then
	bowtie2_options="${bowtie2_options} -X 1000"
fi
#Set number of secondary alignments for genrich
if [[ ! -z ${secondary_alignments} ]]
then
	bowtie2_options="${bowtie2_options} -k ${secondary_alignments}"
fi
#Allow dovetail aignments
if [ ${dovetail} = "TRUE" ]
then
	bowtie2_options="${bowtie2_options} --dovetail"
fi

#Set samtools sort options
sort_options="-@ ${sort_threads}"
if [ ${sort_by_read_names} = "TRUE" ]
then
	sort_options="${sort_options} -n"
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
	path2="$(pwd | sed s/data\\/.*/data/)/${i/_*}/${i/*_/}"
	version=$(ls ${path2}/ref/${i/*_/}-v*.fa | sed s/.*\-v// | sed s/.fa//)
	bowtie2_index="${path2}/ref/bowtie2/${i/*_/}-v${version}"
	#Set bam name
	bam="${sample}_ref_${i}-v${version}.bam"

	#Check if the bam index already exists, if it does, then indicates data was already successfully aligned
	if [ -f ${bam} ]
	then
		echo "Existing bam file found, skipping" to mark duplicates""
		echo "To rerun this step, delete ${bam} and resubmit"
	elif [[ ! -f ${bam} && ${PE} = "TRUE" ]]
	then
		#Align paired-end data
		echo "Aligning to ${i}-v${version}"
		bowtie2 \
			${bowtie2_options} \
			-x ${bowtie2_index} \
			-1 ${t1} \
			-2 ${t2} | samtools view -@ ${compression_threads} -bSh | samtools sort ${sort_options} > ${bam}
		if [ ${sort_by_read_names} = "TRUE" ]
		then
			#Index the bam file
			echo "Indexing ${bam}"	
			samtools index ${bam}
		fi
	elif [[ ! -f ${bam} && ${PE} = "FALSE" ]]
	then
		#Align single-end data
		echo "Aligning to ${i}-v${version}"
		bowtie2 \
			${bowtie2_options} \
			-x ${bowtie2_index} \
			-U ${t1} | samtools view -@ ${compression_threads} -bSh | samtools sort ${sort_options} > ${bam}
		if [ ${sort_by_read_names} = "TRUE" ]
		then
			#Index the bam file
			echo "Indexing ${bam}"	
			samtools index ${bam}
		fi
	fi
	echo "Alignment to ${i}-v${version} complete"
done

echo "Done"

