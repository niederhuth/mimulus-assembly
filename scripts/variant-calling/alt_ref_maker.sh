#!/bin/bash --login
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --job-name alt_ref_maker
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

for i in ${genomes}
do
	echo "Working on ${sample} & ${i}"
	#Set variables
	path3=$(pwd | sed s/${species}\\/.*/${species}\\/${i}/)
	version=$(ls ${path3}/ref/${i}-v*.fa | sed s/.*\-v// | sed s/.fa//) 
	ref="${path3}/ref/${i}-v${version}.fa"
	vcf="${path2}/${sample}-${i}.vcf.gz"
	tmp_vcf="${path2}/tmp.vcf.gz"
	filtered_vcf="${path2}/${sample}-${i}-filtered.vcf.gz"
	alt_ref="${path2}/${sample}-${i}-v${version}.fa"
	
	#Variant Filtration. I Recommend this to remove heterozygous calls and cut-down on low-quality
	#variants that might affect mapping and downstream analysis
	echo "Checking for filtered vcf file"
	if [ -s ${filtered_vcf} >/dev/null 2>&1
	then
		echo "Filtered vcf found, proceeding to next step"
		echo "To rerun this step, delete ${filtered_vcf} and resubmit"
	else
		echo "Filtering variants with VariantFiltration"
		#Currently set to remove heterozygous and low-quality genotype calls
		gatk --java-options ${java_options} VariantFiltration \
			--variant ${vcf} \
			--output ${tmp_vcf} \
			--genotype-filter-expression "isHet == 1" --genotype-filter-name "isHet" \
			--genotype-filter-expression "GQ < 30" --genotype-filter-name "lowGQ"
		echo "Removing filtered variants and setting heterozygous to no-call"
		gatk SelectVariants \
			--variant ${tmp_vcf} \
			--exclude-filtered \
			--set-filtered-gt-to-nocall \
			-output ${filtered_vcf}
		rm ${tmp_vcf} ${tmp_vcf}.tbi
	fi
	#Create alternative reference
	if [ -s ${alt_ref} ]
	then
		echo "Variant-substituted reference for ${sample} & ${i} found"
		echo "To rerun this step, delete ${alt_ref} and resubmit"
	else
		echo "Running gatk FastaAlternativeReferenceMaker"
		gatk --java-options ${java_options} FastaAlternateReferenceMaker \
			--output ${alt_ref} \
			--reference ${ref} \
			--variant ${filtered_vcf}
		#Need to fix fasta names, because keeping original names would be too logical, wouldn't it GATK?
		rm ${alt_ref}.fai $(echo ${alt_ref} | sed s/.fa$/.dict/)
		sed s/\>.*\ /\>/ ${alt_ref} | sed s/\:.*// > tmp
		mv tmp ${alt_ref}
	fi
done

