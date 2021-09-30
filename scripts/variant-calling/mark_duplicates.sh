#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --job-name mark_duplicates
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
java_options="-Xmx16G"

#In general dont change this, unless using a similar datatype
#This should match the dataype in the misc/samples.csv file
datatype="wgs"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/variant-calling/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/variant-calling/lib:${LD_LIBRARY_PATH}"
#Path to picard
picard="${conda}/envs/variant-calling/share/picard-2.23.9-0/picard.jar"

#Other variables, these should not have to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/^.*\\///)
path2=${datatype}

#Get list of genomes
genomes=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $5 == d) print $7}' \
	${path1}/samples.csv)

#Mark Duplicates
for i in ${genomes}
do
	bam="${path2}/${sample}-${i}.bam"
	md="${path2}/${sample}-${i}-md.bam"
	metrics="${path2}/${sample}-${i}-md-metrics.txt"
	#Mark Duplicates
	if [ -s ${md} ]
	then
		echo "Existing marked duplicate bam found"
		echo "To rerun this step, delete ${md} and resubmit"
	else	
		echo "Marking duplicates for ${sample}-${i}"
		java ${java_options} -jar ${picard} MarkDuplicates \
			-I ${bam} \
			-O ${md} \
			-M ${metrics}
		echo "Indexing ${sample}-${i} marked duplicate bam"
		samtools index ${md}
		#Alignment Stats
		echo "Getting ${sample}-${i} marked duplicate alignment stats"
		samtools flagstat ${md} > ${md}.flagstats
		samtools stats ${md} > ${md}.stats
	fi
done

echo "Done"

