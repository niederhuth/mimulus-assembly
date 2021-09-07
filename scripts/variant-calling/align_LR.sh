#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=500GB
#SBATCH --job-name racon
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=40
rounds=3
datatype="ont"
input="" #Can set to empty and script will find fasta in directory submitted

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/lr_variants/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/lr_variants/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${genotype}\\/// | sed s/\\/.*//)
assembly=$(pwd | sed s/^.*\\///)
reads="fastq/${datatype}/clean.fastq.gz"
path2="LRvariants"

#Check for reads file
if [ -f ${reads} ]
then
	echo "${reads} found"
else
	echo "Could not find ${reads}"
	echo "Check if reads are cleaned reads are present"
	exit 1
fi

#Change preset based on datatype
if [ ${datatype} = "ont" ]
then
	preset="map-ont"
elif [ ${datatype} = "ont-cor" ]
then
	preset="map-ont"
elif [ ${datatype} = "pac" ]
then
	preset="map-pb"
elif [ ${datatype} = "pac-cor" ]
then
	preset="map-pb"
elif [ ${datatype} = "hifi" ]
then
	preset="map-pb"
else
	echo "Do not recognize ${datatype}"
	echo "Please check and resubmit"
fi

#Get list of genomes
genomes=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $5 == d) print $7}' \
	${path1}/samples.csv)
 
#Create output directory
if [ -d ${path2} ]
do
	cd ${path2}
else
	mkdir ${path2}
	cd ${path2}
fi

#Align data & Mark Duplicates
for i in ${genomes}
do
	#Set variables
	path4=$(pwd | sed s/${species}\\/.*/${species}\\/${i}/)
	version=$(ls ${path4}/ref/${i}-v*.fa | sed s/.*\-v// | sed s/.fa//)
	ref="${path4}/ref/${i}-v${version}.fa"
	sam="${sample}-${i}.sam"
	bam="${sample}-${i}.bam"

	#Align reads to assembly
	if [ -s ${sam} ]
	then
		echo "Aligned reads found, proceeding to bam conversion & sorting."
		echo "To repeat this step, delete ${sam} and resubmit."
	else
		echo "Aligning reads to assembly"
		minimap2 \
			-a \
			-t ${threads} \
			-x ${preset} \
			${ref} \
			${path1}/${reads} > ${sam}
	fi

	if [ -s ${bam} ]
	then
		echo "Bam file found, proceeding to read-depth histogram."
		echo "To repeat this step, delete ${bam} and resubmit."
	else
		echo "Sorting and converting to bam"
		samtools view -@ ${threads2} -bSh aligned.sam | samtools sort -@ ${threads2} > ${bam}
		echo "Indexing aligned.bam"
		samtools index ${bam}
	fi
done

echo "Done"
