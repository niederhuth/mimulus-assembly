#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20GB
#SBATCH --job-name gatk
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=4
java_options="-Xmx16G"

#In general dont change this, unless using a similar datatype
#This should match the dataype in the misc/samples.csv file
datatype="wgs"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/gatk/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/gatk/lib:${LD_LIBRARY_PATH}"
#Path to GATK Jar file
export GATK_LOCAL_JAR="${conda}/envs/gatk/share/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar"

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

#Call variants using haplotype caller
for i in ${genomes}
do
	echo "Calling variants for ${sample} & ${i}"
	#Set variables
	path3=$(pwd | sed s/${species}\\/.*/${species}\\/${i}/)
	version=$(ls ${path3}/ref/${i}-v*.fa | sed s/.*\-v// | sed s/.fa//) 
	ref="${path3}/ref/${i}-v${version}.fa"
	known_sites="${path3}/ref/known_sites/*vcf.gz"
	bam="${path2}/${sample}-${i}-md.bam"
	vcf="${path2}/${sample}-${i}.vcf.gz"
	
	if [ -s ${known_sites} ]
	then
		echo "Known sites found"
		bqsr="${path2}/${sample}-${i}_bqsr.table"
		bam2="${path2}/${sample}-${i}_bqsr.bam"
		if [ -f ${bam2} ]
		then
			echo "Recalibrated bam file found, skipping base recalibration"
			echo "To rerun this analysis, delete ${bam2} and resubmit"
			echo "Proceeding to variant calling"
		else
			echo "Recalibrating base quality scores"
			echo "Creating recalibrated base table"
			gatk --java-options ${java_options} BaseRecalibrator \
   				-R ${ref} \
				-I ${bam} \
   				--known-sites ${known_sites} \
   				-O ${bqsr}
			echo "Creating recalibrated bam file"
			gatk --java-options ${java_options} ApplyBQSR \
				-R ${ref} \
   				-I ${bam} \
   				--bqsr-recal-file ${bqsr} \
   				-O ${bam2}
			echo "Base recalibration done"
		fi
		bam="${bam2}"
	else
		echo "Known sites not available, skipping base recalibration"
		echo "Proceeding to HaplotypeCaller"
	fi
	if [ -s ${vcf} ]
	then
		echo "${vcf} already exists, skipping variant calling"
		echo "To rerun this analysis, delete ${vcf} and rerun"
	else
		echo "Running gatk HaplotypeCaller"
		gatk --java-options ${java_options} HaplotypeCaller \
			--native-pair-hmm-threads ${threads} \
			-R ${ref} \
			-I ${bam} \
			-O ${vcf}
		echo "Variant calling for ${sample} & ${i} complete"
	fi
done

echo "Done"

