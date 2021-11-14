#!/bin/bash --login
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=44
#SBATCH --mem=80GB
#SBATCH --job-name align_wgs
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=40
threads2=4
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
picard="${conda}/envs/variant-calling/share/picard-*/picard.jar"
#Path to trimmomatic fastas 
adapter_path="${conda}/envs/variant-calling/share/trimmomatic/adapters"

#Other variables, these should not have to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/^.*\\///)
path2=${datatype}

#Adapter fasta, set automatically from misc/samples.csv
adapters=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $5 == d) print $8}' \
	${path1}/samples.csv)

#Fastq files, these should not have to be changed, but should set automatically
path3="fastq/${datatype}"
r1="${path3}/combined.1.fastq.gz"
r2="${path3}/combined.2.fastq.gz"
t1="${path3}/trimmed.1.fastq.gz"
t2="${path3}/trimmed.2.fastq.gz"
t3="${path3}/trimmed.1.single.fastq.gz"
t4="${path3}/trimmed.2.single.fastq.gz"
if ls ${path3}/*_R1_001.fastq.gz >/dev/null 2>&1
then
	if ls ${path3}/*_R2_001.fastq.gz >/dev/null 2>&1
	then
		echo "Data is Paired-end"
		PE="TRUE"
		if [ -f ${t1} ]
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path3}/*_R1_001.fastq.gz > $r1
			cat ${path3}/*_R2_001.fastq.gz > $r2
		fi
	else
		echo "Data is Single-end"
		PE="FALSE"
		if [ -f ${t1} ]
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path3}/*_R1_001.fastq.gz > $r1
		fi
	fi
elif ls ${path3}/*_1.fastq.gz >/dev/null 2>&1
then
	if ls ${path3}/*_2.fastq.gz >/dev/null 2>&1	
	then
		echo "Data is Paired-end"
		PE="TRUE"
		if [ -f ${t1} ]
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path3}/*_1.fastq.gz > $r1
			cat ${path3}/*_2.fastq.gz > $r2
		fi
	else
		echo "Data is Single-end"
		PE="FALSE"
		if [ -f ${t1} ]
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path3}/*_1.fastq.gz > $r1
		fi
	fi
else
	echo "Data Missing"
fi

#Trim & QC reads
if [ -f ${t1} ]
then
	if [ ${PE} = "TRUE" ]
	then
		echo "To rerun this step, please delete ${t1} & ${t2} and resubmit"
	else
		echo "To rerun this step, please delete ${t1} and resubmit"
	fi
else
	if [ ${PE} = "TRUE" ]
	then
		echo "Running trimmomatic PE"
		trimmomatic PE \
			-threads ${threads} \
			-phred33 \
			-trimlog ${path3}/trim_log.txt \
			-summary ${path3}/trim_summary.txt \
			${r1} ${r2} ${t1} ${t3} ${t2} ${t4} \
			ILLUMINACLIP:${adapter_path}/${adapters}.fa:2:30:10:4:TRUE \
			LEADING:3 \
			TRAILING:3 \
			SLIDINGWINDOW:4:15 \
			MINLEN:30
		echo "Running fastqc"
		mkdir ${path3}/fastqc
		fastqc -t ${threads} -o ${path3}/fastqc/ ${t1} ${t2} ${r1} ${r2}
	elif [ ${PE} = "FALSE" ]
	then
		echo "Running trimmomatic SE"
		trimmomatic SE \
			-threads ${threads} \
			-phred33 \
			-trimlog ${path3}/trim_log.txt \
			-summary ${path3}/trim_summary.txt \
			${r1} ${t1} \
			ILLUMINACLIP:${adapter_path}/${adapters}.fa:2:30:10:4:TRUE \
			LEADING:3 \
			TRAILING:3 \
			SLIDINGWINDOW:4:15 \
			MINLEN:30 
		echo "Running fastqc"
		mkdir ${path3}/fastqc
		fastqc -t ${threads} -o ${path3}/fastqc/ ${t1} ${r1}
	fi
fi
rm ${r1} ${r2}

#Define Read Group
if ls ${path3}/SRR*.fastq.gz >/dev/null 2>&1
then
	ID=$(zcat ${t1} | head -1 | cut -d ' ' -f1 | sed s/\@//)
	PU=$(zcat ${t1} | head -1 | cut -d ' ' -f1 | sed s/\@//)
	SM=$(pwd | sed s/^.*\\///)
	PL="ILLUMINA"
	LB="lib1"
else
	ID=$(zcat ${t1} | head -1 | cut -d ':' -f 3,4 | tr ':' '.')
	PU=$(zcat ${t1} | head -1 | cut -d ':' -f 3,4,10 | tr ':' '.')
	SM=$(pwd | sed s/^.*\\///)
	PL="ILLUMINA"
	LB="lib1"
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
mkdir ${path2}

#Align data & Mark Duplicates
for i in ${genomes}
do
	#Set variables
	path4=$(pwd | sed s/${species}\\/.*/${species}\\/${i}/)
	version=$(ls ${path4}/ref/${i}-v*.fa | sed s/.*\-v// | sed s/.fa//)
	ref="${path4}/ref/bwa/${i}-v${version}.fa"
	bam="${path2}/${sample}-${i}.bam"
	md="${path2}/${sample}-${i}-md.bam"
	metrics="${path2}/${sample}-${i}-md-metrics.txt"
	
	#Align Data
	if [ -s ${bam} ]
	then
		echo "Existing bam file found, skipping to mark duplicates"
		echo "To rerun this step, delete ${bam} and resubmit"
	else
		echo "Running BWA for ${sample} to ${i} genome"
		if [ ${PE} = "TRUE" ]
		then
			bwa mem -t ${threads} \
				-R "@RG\tID:${ID}\tLB:${LB}\tPL:${PL}\tSM:${SM}\tPU:${PU}" \
				-M ${ref} ${t1} ${t2} | \
				samtools view -@ ${threads2} -bSh | \
				samtools sort -@ ${threads2} > ${bam}
		elif [ ${PE} = "False" ]
		then	
			bwa mem -t ${threads} \
				-R "@RG\tID:${ID}\tLB:${LB}\tPL:${PL}\tSM:${SM}\tPU:${PU}" \
				-M ${ref} ${t1} | \
				samtools view -@ ${threads2} -bSh | \
				samtools sort -@ ${threads2} > ${bam}
		fi
		echo "Indexing ${sample}-${i} bam file"	
		samtools index ${bam}
		echo "Getting ${sample}-${i} alignment stats"
		samtools flagstat ${bam} > ${bam}.flagstats
		samtools stat ${bam} > ${bam}.stats
	fi
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

