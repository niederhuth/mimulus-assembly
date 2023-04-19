#!/bin/bash --login
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=40GB
#SBATCH --job-name align_atac
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20
genrich_mode=FALSE #If using Genrich as peak-caller, set to "TRUE" to configure for
				   #Genrich appropriate output. As of now, only recommend using bowtie2
mark_dups=TRUE #Mark duplicate reads
keep_tmp=FALSE #Keep temporary bam files during filtering steps: "TRUE" or "FALSE"
rm_organelle="TRUE" #Remove organellar reads: "TRUE" or "FALSE"
rm_dups="TRUE" #Remove duplicate reads: "TRUE" or "FALSE"
rm_multimappers="TRUE" #Remove multimapping reads: "TRUE" or "FALSE"
mito_name="chrM" #Change this to the name of mitochondrial chromosome
chloro_name="chrC" #Change this to the name of the chloroplast chromosome
java_options="-Xmx16G"

#Which aligner to use: "bowtie", "bowtie2", "bwa", "by_length"
#If set to "by_length", the read length will first be checked and
#bowtie used for reads <= 50 nt, bowtie2 for reads > 50 nt
#Otherwise, the aligner specified will be used, I in general don't recommend bwa for this.
#I have not yet worked out the proper way to handle or filter bwa mapped reads for ATAC
aligner="bowtie2"

#Mapping quality score for filtering multimapping reads
#Recommend not changing, but if you want to apply your own filter, can set here
#If = "auto", will set using following defaults based on aligner:
#bowtie2: 4 - See: https://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html
#bowtie: 255 - This is very likely too conservative, but going with it for now
#bwa: 0 - I have not figured out proper setting for bwa yet.
mapq="auto"

#In general dont change this, unless using a similar datatype, e.g. DNase-seq, ChIP, etc
#This should match the dataype in the misc/samples.csv file
#Recommended options: "atac" "dnase" "mnase" "faire" "chip"
datatype="atac" 

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/atac/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/atac/lib:${LD_LIBRARY_PATH}"
#Path to picard
picard="${conda}/envs/atac/share/picard-2.23.9-0/picard.jar"
#Path to trimmomatic fastas 
adapter_path="${conda}/envs/atac/share/trimmomatic/adapters"

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
		if ls ${t1} >/dev/null 2>&1
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path3}/*_R1_001.fastq.gz > $r1
			cat ${path3}/*_R2_001.fastq.gz > $r2
		fi
	else
		echo "Data is Single-end"
		echo "Single-end ${datatype}? If this is wrong, double check and restart"
		PE="FALSE"
		if ls ${t1} >/dev/null 2>&1
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
		if ls ${t1} >/dev/null 2>&1
		then
			echo "Trimmed reads found, skipping trimming"
		else
			cat ${path3}/*_1.fastq.gz > $r1
			cat ${path3}/*_2.fastq.gz > $r2
		fi
	else
		echo "Data is Single-end"
		echo "Single-end ${datatype}? If this is wrong, double check and restart"
		PE="FALSE"
		if ls ${t1} >/dev/null 2>&1
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
if ls ${t1} >/dev/null 2>&1
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
			MINLEN:30 
		echo "Running fastqc"
		mkdir ${path3}/fastqc
		fastqc -t ${threads} -o ${path3}/fastqc/ ${t1} ${r1}
	fi
fi

#Read length
length=$(expr $(zcat ${r1} | head -2 | tail -1 | wc -c) - 1)
echo "Read length is ${length}"
#Set aligner
if [ ${aligner} = "by_length" ]
then
	echo "Aligner set to by_length"
	if [ ${length} -gt 50 ]
	then
		echo "Read length >50 nt, using bowtie2"
		aligner="bowtie2"
	elif [ ${length} -le 50 ]
	then
		echo "Read length <50 nt, using bowtie"
		aligner="bowtie"
	fi
else
	echo "Using ${aligner}"
fi

#Set aligner options
if [ ${PE} = "TRUE" ]
then
	if [ ${genrich_mode} = "TRUE" ]
	then
		bowtie2_options="--very-sensitive -X 1000 -k 10 --no-mixed --dovetail"
		bowtie_options="-X 1000"
		bwa_options=""
	else
		bowtie2_options="--very-sensitive -X 1000 --no-mixed --dovetail"
		bowtie_options="-m 1 -v 2 --best --strata --allow-contain -X 1000"
		bwa_options=""
	fi
else
	if [ ${genrich_mode} = "TRUE" ]
	then
		bowtie2_options="--very-sensitive -k 10"
		bowtie_options=""
		bwa_options=""
	else
		bowtie2_options="--very-sensitive"
		bowtie_options="-m 1 -v 2 --best --strata"
		bwa_options=""
fi

#Set mapq filter
if [ ${mapq} = "auto" ]
then
	if [ ${aligner} = "bowtie2" ]
	then
		mapq=4 # See: https://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html
	elif [ ${aligner} = "bowtie" ]
	then
		mapq=255 #This seems too conservative
	elif [ ${aligner} = "bwa" ]
	then
		mapq=30
	fi
	echo "mapq set to ${mapq}"
	echo "Note, this will only apply if rm_multimappers = TRUE"
fi

#Remove combined read files
rm ${r1} ${r2}

#Define Read Group
ID=$(zcat ${t1} | head -1 | cut -d ':' -f 3,4 | tr ':' '.')
PU=$(zcat ${t1} | head -1 | cut -d ':' -f 3,4,10 | tr ':' '.')
SM=$(pwd | sed s/^.*\\///)
PL="ILLUMINA"
LB="lib1"

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
	bam="${path2}/${sample}-${i}.bam"
	#Align Data
	if ls ${bam} >/dev/null 2>&1
	then
		echo "Existing bam file found, skipping to mark duplicates"
		echo "To rerun this step, delete ${bam} and resubmit"
	else
		if [ ${PE} = "TRUE" ]
		then
			if [ ${aligner} = "bowtie2" ]
			then
				echo "Aligner set to ${aligner}"
				echo "Running bowtie2 for ${sample} to ${i} genome"
				ref="${path4}/ref/bowtie2/${i}-v${version}"
				bowtie2 \
					-p ${threads} \
					${bowtie2_options} \
					-x ${ref} \
					-1 ${t1} \
					-2 ${t2} | samtools view -@ 4 -bSh | samtools sort -@ 4 > ${bam}
			elif [ ${aligner} = "bowtie" ]
			then
				echo "Aligner set to ${aligner}"
				echo "Running bowtie for ${sample} to ${i} genome"
				ref="${path4}/ref/bowtie/${i}-v${version}"
				bowtie \
					-p ${threads} \
					${bowtie_options} \
					-S \
					-x ${ref} \
					-1 ${t1} \
					-2 ${t2} | samtools view -@ 4 -bSh | samtools sort -@ 4 > ${bam}
			elif [ ${aligner} = "bwa" ]
			then
				echo "Aligner set to ${aligner}"
				echo "Running bwa for ${sample} to ${i} genome"
				ref="${path4}/ref/bwa/${i}-v${version}.fa"
				bwa mem \
					-t ${threads} \
					-R "@RG\tID:${ID}\tLB:${LB}\tPL:${PL}\tSM:${SM}\tPU:${PU}" \
					-M ${ref} ${t1} ${t2} | samtools view -@ 4 -bSh | samtools sort -@ 4 > ${bam}
			fi
		elif [ ${PE} = "FALSE" ]
		then
			if [ ${aligner} = "bowtie2" ]
			then
				echo "Aligner set to ${aligner}"
				echo "Running bowtie2 for ${sample} to ${i} genome"
				ref="${path4}/ref/bowtie2/${i}-v${version}"
				bowtie2 \
					-p ${threads} \
					${bowtie2_options} \
					-x ${ref} \
					-U ${t1} | samtools view -@ 4 -bSh | samtools sort -@ 4 > ${bam}
			elif [ ${aligner} = "bowtie" ]
			then
				echo "Aligner set to ${aligner}"
				echo "Running bowtie for ${sample} to ${i} genome"
				ref="${path4}/ref/bowtie/${i}-v${version}"
				bowtie \
					-p ${threads} \
					${bowtie_options} \
					-S \
					-x ${ref} ${t1} | samtools view -@ 4 -bSh | samtools sort -@ 4 > ${bam}
			elif [ ${aligner} = "bwa" ]
			then
				echo "Aligner set to ${aligner}"
				echo "Running bwa for ${sample} to ${i} genome"
				ref="${path4}/ref/bwa/${i}-v${version}.fa"
				bwa mem \
					-t ${threads} \
					-R "@RG\tID:${ID}\tLB:${LB}\tPL:${PL}\tSM:${SM}\tPU:${PU}" \
					-M ${ref} ${t1} | samtools view -@ 4 -bSh | samtools sort -@ 4 > ${bam}
			fi
		fi
		echo "Indexing ${sample}-${i} bam file"	
		samtools index ${bam}
		echo "Getting ${sample}-${i} alignment stats"
		samtools flagstat ${bam} > ${bam}.flagstats
		samtools stat ${bam} > ${bam}.stat
	fi
	#Mark Duplicates
	if [ ${mark_dups} = "TRUE" ]
	then
		echo "mark_dups is set to TRUE"
		md="${path2}/${sample}-${i}-md.bam"
		metrics="${path2}/${sample}-${i}-md-metrics.txt"
		if ls ${md} >/dev/null 2>&1
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
	else
		echo "mark_dups is set to FALSE"
		md=${bam}
	fi
	#Filter Data
	if [ ${rm_dups} = "TRUE" ] || [ ${rm_multimappers} = "TRUE" ] || [ ${rm_organelle} = "TRUE" ]
	then
		echo "Filtering bam file"
		filtered="${path2}/${sample}-${i}"
		output_bam="${path2}/${sample}-${i}"
		if [ ${rm_organelle} = "TRUE" ]
		then
			filtered="${filtered}-rm${mito_name}${chloro_name}"
		fi
		if [ ${rm_dups} = "TRUE" ]
		then
			filtered="${filtered}-rmDups"
		fi
		if [ ${rm_multimappers} = "TRUE" ]
		then
			filtered="${filtered}-rmMulti"
		fi
		if ls ${filtered}.bam >dev/null 2>&1
		then
			echo "Filtered bam file found, skipping"
			echo "To rerun this step, delete ${filtered}.bam and resubmit"
		else
			input_bam=${md}
			if [ ${rm_organelle} = "TRUE" ]
			then
				echo "Filtering organellar reads"
				output_bam="${output_bam}-rm${mito_name}${chloro_name}"
				samtools view -@ ${threads} -h ${input_bam} | \
				grep -v ${mito_name} | grep -v ${chloro_name} | \
				samtools sort -@ ${threads} -O bam -o ${output_bam}.bam
				echo "Indexing ${output_bam}.bam"
				samtools index ${output_bam}.bam
				#Alignment Stats
				echo "Getting ${output_bam}.bam alignment stats"
				samtools flagstat ${output_bam}.bam > ${output_bam}.bam.flagstats
				samtools stats ${output_bam}.bam > ${output_bam}.bam.stats
				if [ ${input_bam} != ${md} ] && [ ${keep_tmp} = "FALSE" ]
				then
					rm ${input_bam}
				fi
				input_bam="${output_bam}.bam"
			fi
			if [ ${rm_dups} = "TRUE" ]
			then
				echo "Filtering out duplicate reads"
				output_bam="${output_bam}-rmDups"
				samtools view -h -b -F 1024 ${input_bam} | \
				samtools sort -@ ${threads} -O bam -o ${output_bam}.bam
				echo "Indexing ${output_bam}.bam"
				samtools index ${output_bam}.bam
				#Alignment Stats
				echo "Getting ${output_bam}.bam alignment stats"
				samtools flagstat ${output_bam}.bam > ${output_bam}.bam.flagstats
				samtools stats ${output_bam}.bam > ${output_bam}.bam.stats
				if [ ${input_bam} != ${md} ] && [ ${keep_tmp} = "FALSE" ]
				then
					rm ${input_bam}
				fi
				input_bam="${output_bam}.bam"
			fi
			if [ ${rm_multimappers} = "TRUE" ]
			then
				echo "Filtering out multimapping reads"
				output_bam="${output_bam}-rmMulti"
				samtools view -h -q ${mapq} ${input_bam} | \
				samtools sort -@ ${threads} -O bam -o ${output_bam}.bam
				echo "Indexing ${output_bam}.bam"
				samtools index ${output_bam}.bam
				#Alignment Stats
				echo "Getting ${output_bam}.bam alignment stats"
				samtools flagstat ${output_bam}.bam > ${output_bam}.bam.flagstats
				samtools stats ${output_bam}.bam > ${output_bam}.bam.stats
				if [ ${input_bam} != ${md} ] && [ ${keep_tmp} = "FALSE" ]
				then
					rm ${input_bam}
				fi
				input_bam="${output_bam}.bam"
			fi
			if [ ${input_bam} != ${md} ] && [ ${input_bam} != "${filtered}.bam" ] && [ ${keep_tmp} = "FALSE" ]
			then
				rm ${input_bam}
			fi
		fi
	else
		echo "No filter options selected, skipping filtering"
	fi
done

echo "Done"

