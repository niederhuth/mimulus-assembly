#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=510GB
#SBATCH --job-name align-rna-SR-annotation
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=40
fasta=$(ls -l *.fa | sed s/.*\ //)
genomeSAindexNbases=13 #14 default
IntronMin=10
IntronMax=5000
BAMsortRAM=7035318712
datatype="rna"


#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/transcript-assembly/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/transcript-assembly/lib:${LD_LIBRARY_PATH}"
#Path to trimmomatic fastas 
adapter_path="${conda}/envs/transcript-assembly/share/trimmomatic/adapters"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="annotation"
assembly=$(pwd | sed s/^.*\\///)
path2="$(pwd | sed s/annotation.*//)/fastq/${datatype}"
path3="SRrna"

#Make & cd to directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#Make STAR index
if [ ! -d STAR ]
then
	echo "Making STAR index"
	STAR \
		--runThreadN ${threads} \
		--runMode genomeGenerate \
		--genomeDir STAR/ \
		--genomeFastaFiles ../${fasta} \
		--genomeSAindexNbases ${genomeSAindexNbases}
else
	echo "STAR index found"
fi

#Adapter fasta, set automatically from misc/samples.csv
adapters=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${condition} \
	-v e=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $4 == d && $5 == e) print $8}' \
	${path1}/samples.csv)

#Fastq files, these should not have to be changed, but should set automatically
r1="${path2}/combined.1.fastq.gz"
r2="${path2}/combined.2.fastq.gz"
t1="${path2}/trimmed.1.fastq.gz"
t2="${path2}/trimmed.2.fastq.gz"
t3="${path2}/trimmed.1.single.fastq.gz"
t4="${path2}/trimmed.2.single.fastq.gz"
if ls ${path2}/*_R1_001.fastq.gz >/dev/null 2>&1
then
	if ls ${path2}/*_R2_001.fastq.gz >/dev/null 2>&1
	then
		echo "Data is Paired-end"
		PE="TRUE"
		if ls ${t1} >/dev/null 2>&1
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path2}/*_R1_001.fastq.gz > $r1
			cat ${path2}/*_R2_001.fastq.gz > $r2
		fi
	else
		echo "Data is Single-end"
		echo "Single-end ${datatype}? If this is wrong, double check and restart"
		PE="FALSE"
		if ls ${t1} >/dev/null 2>&1
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path2}/*_R1_001.fastq.gz > $r1
		fi
	fi
elif ls ${path2}/*_1.fastq.gz >/dev/null 2>&1
then
	if ls ${path2}/*_2.fastq.gz >/dev/null 2>&1	
	then
		echo "Data is Paired-end"
		PE="TRUE"
		if ls ${t1} >/dev/null 2>&1
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path2}/*_1.fastq.gz > $r1
			cat ${path2}/*_2.fastq.gz > $r2
		fi
	else
		echo "Data is Single-end"
		echo "Single-end ${datatype}? If this is wrong, double check and restart"
		PE="FALSE"
		if ls ${t1} >/dev/null 2>&1
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path2}/*_1.fastq.gz > $r1
		fi
	fi
elif ls ${path2}/SRR*.fastq.gz >/dev/null 2>&1
then
	echo "Data is Single-end"
	echo "Single-end ${datatype}? If this is wrong, double check and restart"
	PE="FALSE"
	if ls ${t1} >/dev/null 2>&1
	then
		echo "Trimmed reads found, skipping trimming"
	else
		cat ${path2}/SRR*.fastq.gz > $r1
	fi
else
	echo "Data Missing"
fi

#Trim & QC reads
if ls ${t1} >/dev/null 2>&1
then
	if [ ${PE} = "TRUE" ]
	then
		echo "To rerun this step, please delete ${t1} & ${t2} and resubmit"
		fastq="${t1} ${t2}"
	else
		echo "To rerun this step, please delete ${t1} and resubmit"
		fastq=${t1}
	fi
else
	if [ ${PE} = "TRUE" ]
	then
		echo "Running trimmomatic PE"
		trimmomatic PE \
			-threads ${threads} \
			-phred33 \
			-trimlog ${path2}/trim_log.txt \
			-summary ${path2}/trim_summary.txt \
			${r1} ${r2} ${t1} ${t3} ${t2} ${t4} \
			ILLUMINACLIP:${adapter_path}/${adapters}.fa:2:30:10:4:TRUE \
			MINLEN:30
		echo "Running fastqc"
		mkdir ${path2}/fastqc
		fastqc -t ${threads} -o ${path2}/fastqc/ ${t1} ${t2} ${r1} ${r2}
		#Note that I am currently ignoring any unpaired reads
		fastq="${t1} ${t2}"
	elif [ ${PE} = "FALSE" ]
	then
		echo "Running trimmomatic SE"
		trimmomatic SE \
			-threads ${threads} \
			-phred33 \
			-trimlog ${path2}/trim_log.txt \
			-summary ${path2}/trim_summary.txt \
			${r1} ${t1} \
			ILLUMINACLIP:${adapter_path}/${adapters}.fa:2:30:10:4:TRUE \
			MINLEN:30 
		echo "Running fastqc"
		mkdir ${path2}/fastqc
		fastqc -t ${threads} -o ${path2}/fastqc/ ${t1} ${r1}
		fastq=${t1}
	fi
fi

#Align with STAR
if [ ${PE} = "TRUE" ]
then
	echo "Running STAR for ${sample} against ${fasta}"
	STAR \
		--runThreadN ${threads} \
		--runMode alignReads \
		--genomeDir STAR \
		--readFilesIn ${t1} ${t2}\
		--readFilesCommand zcat \
		--alignIntronMin ${IntronMin} \
		--alignIntronMax ${IntronMax} \
		--outSAMtype BAM SortedByCoordinate \
		--limitBAMsortRAM ${BAMsortRAM}
elif [ ${PE} = "FALSE" ]
then
	echo "Running STAR for ${sample} against ${fasta}"
	STAR \
		--runThreadN ${threads} \
		--runMode alignReads \
		--genomeDir STAR \
		--readFilesIn ${t1} \
		--readFilesCommand zcat \
		--alignIntronMin ${IntronMin} \
		--alignIntronMax ${IntronMax} \
		--outSAMtype BAM SortedByCoordinate \
		--limitBAMsortRAM ${BAMsortRAM}
fi

echo "Done"
