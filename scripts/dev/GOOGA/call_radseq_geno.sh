#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=50GB
#SBATCH --job-name call_radseq_geno
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20 #threads for mpileup
threads2=5 #threads for call
max_alleles=2 #Max allele count for filtering, recommend 2 for now
maf=0.1 #Minor allele frequency, assuming this data is from F2s or F3s, so keeping high
max_non_ref_af=0.9 #Maximum non-ref allele frequencyi, use to eliminate allels without ref allele
max_missing=0.2 #Max number of samples can be missing calls

#In general dont change this, unless using a similar datatype
#This should match the dataype in the misc/samples.csv file
datatype="radseq"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/scaffolding/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/scaffolding/lib:${LD_LIBRARY_PATH}"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
assembly=$(pwd | sed s/^.*\\///)
path2=$(pwd | sed s/${genotype}\\/${sample}\\/.*/${genotype}\\/${sample}/)
path3=${path1}/genetic_map/${datatype}

#Look for fasta file, there can only be one!
if [ -z ${fasta} ]
then
	echo "No input fasta provided, looking for fasta"
	if ls *.fa >/dev/null 2>&1
	then
		fasta=$(ls *fa | sed s/.*\ //)
		echo "Fasta file ${input} found"
	elif ls *.fasta >/dev/null 2>&1
	then
		fasta=$(ls *fasta | sed s/.*\ //)
		echo "Fasta file ${input} found"
	elif ls *.fna >/dev/null 2>&1
	then
		fasta=$(ls *fna | sed s/.*\ //)
		echo "Fasta file ${input} found"
	else
		echo "No fasta file found, please check and restart"
	fi
else
	echo "Input fasta: ${input}"
fi

for i in ${path3}/*
do
	name=$(echo ${i} | sed s/^.*\\/// | sed s/_SRR.txt//)
	echo "Working on dataset ${i}"
	cd ${name}
	if [ ! -d bam_files ]
	then
		mkdir bam_files
		mv *bam bam_files/
		mv *bam.bai bam_files/
	fi
	#Call genotypes
	if [ ! -s ${name}.vcf ]
	then
		echo "Calling variants with bcftools"
		bcftools mpileup --threads ${threads} --ignore-RG --fasta-ref ../${fasta} bam_files/*.bam |\
		bcftools call --threads ${threads2} -m -v -o ${name}.vcf
	fi
	#Filter sites
	echo "Filtering variants"
	vcftools \
		--vcf ${name}.vcf \
		--out ${name} \
		--recode \
		--remove-indels \
		--max-alleles ${max_alleles} \
		--maf ${maf} \
		--max-non-ref-af ${max_non_ref_af} \
		--max-missing ${max_missing}
	#Modify genotypes
	echo "Reformatting genotypes"
	cat ${name}.recode.vcf | \
	sed s/\\.\\/\\.\\:/NN\:/ | \
	sed s/0\\/0\\:/AA\:/g | \
	sed s/0\\/1\\:/AB\:/g | \
	sed s/1\\/1\\:/BB\:/g > ${name}.modified.recode.vcf
	if [ -s ${name}.modified.recode.vcf ]
	then
		vcf=${name}.modified.recode.vcf
		#Make g_files directory if not present
		if [ ! -d g_files ]
		then
			mkdir g_files
		fi
		#Get sample number from vcf
		ncol=$(grep "#CHROM" ${name}.modified.recode.vcf | awk '{print NF; exit}')
		#we start at column 10, where first sample is
		a=10
		echo "Outputing individual sample genotypes"
		#Loop over each sample column and output a g_gile
		until [ ${a} -gt ${ncol} ]
		do
			#Get sample for that column
			sample=$(grep "#CHROM" ${vcf} | cut -f${a} | sed s/bam_files\\/// | sed s/.bam//)
			#Cut the propter columsn and output g_file
			cut -f1,2,${a} ${vcf} | grep -v \# | sed s/\:.*// |\
			awk -v OFS="\t" -v x=${sample} '{print x,$0}' > g_files/${sample}_genotype.txt
			#Add 1 to the column
			a=$(expr ${a} + 1)
		done
	fi
	echo "Done processing ${name}"
	cd ../
done

echo "Done"
