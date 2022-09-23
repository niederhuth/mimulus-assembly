#!/bin/bash --login
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50GB
#SBATCH --job-name star_rnaseq
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/expression/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/expression/lib:${LD_LIBRARY_PATH}"

#Other variables, these should not have to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition=""
datatype="rna"

#Adapter fasta, set automatically from misc/samples.csv
adapter_path="${conda}/envs/expression/share/trimmomatic/adapters" #Path to trimmomatic fastas 
adapters=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $5 == d) print $8}' \
	${path1}/samples.csv)

#Fastq files, these should not have to be changed, but should set automatically
path2="fastq/${datatype}"
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

#Get list of genomes
genomes=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $5 == d) print $7}' \
	${path1}/samples.csv)

#Run Star
for i in ${genomes}
do
	echo "Running STAR for ${sample} against ${i}"
	path3=${datatype}_${i}
	mkdir ${path3}
	index="$(pwd | sed s/${species}.*/${species}/)/${genotype}/ref/STAR"
	STAR \
		--runThreadN ${threads} \
		--runMode alignReads \
		--genomeDir ${index} \
		--readFilesIn ${fastq} \
		--outFileNamePrefix ${path3}/ \
		--readFilesCommand zcat \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMstrandField intronMotif \
		--outFilterType BySJout \
		--outFilterMultimapNmax 10 \
		--alignSJoverhangMin 5 \
		--alignSJDBoverhangMin 3 \
		--alignIntronMin 20 \
		--alignIntronMax 0 \
		--outFilterScoreMinOverLread 0.33 \
		--outFilterMatchNminOverLread 0.33 \
		--outFilterMismatchNmax 10 \
		--outFilterMismatchNoverReadLmax 0.1 \
		--quantMode GeneCounts
done

echo "Done"


