#!/bin/bash --login
#SBATCH --time=128:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100GB
#SBATCH --job-name genespace-inversion-breakpoints
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20
threads2=4
extend=100000
ref= #reference species_genotype to map to, e.g. Mguttatus_S1, if left blank, get from misc/samples.csv
genespace= #path to genespace syntenicBlocks.txt.gz, if left to blank look in data/comparative

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/variant-calling/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/variant-calling/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
path2=$(pwd | sed s/data.*/data/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="breakpoints"
datatype="genome"
path3="breakpoints"

#Check for and make/cd working directory
if [ -d ${path3} ]
then
	cd ${path3}
	mkdir regions indexes
else
	mkdir ${path3}
	cd ${path3}
	mkdir regions indexes
fi

#Set path to genespace
if [ -z ${genespace} ]
then
	genespace="${path2}/comparative/genespace/results/syntenicBlocks.txt.gz"
fi

#Get list of genomes
if [ -z ${ref} ]
then
	ref=$(awk -v FS="," \
		-v a=${species} \
		-v b=${genotype} \
		-v c=${sample} \
		-v d=${condition} \
		-v e=${datatype} \
		'{if ($1 == a && $2 == b && $3 == c && $4 == d && $5 == e) print $7}' \
		${path1}/samples.csv)
fi

#Make the breakpoints bed file for the genotype
zcat ${genespace} | awk -v a=${genotype} -v b=${ref/*_/} '$1==a && $2==b' | cut -f4 | sort | uniq | while read seq
do
	chr=$(echo ${seq} | cut -d ' ' -f4)
	start2=
	stop2=
	sign2=
	zcat ${genespace} | awk -v a=${genotype} -v b=${ref/*_/} -v c=${chr} '$1==a && $2==b && $3==c && $4==c' | \
	cut -f3,6,7,20 | sort -k2n | while read line
	do
		start=$(echo ${line} | cut -d ' ' -f2)
		stop=$(echo ${line} | cut -d ' ' -f3)
		sign=$(echo ${line} | cut -d ' ' -f4)
		if [ ${sign} == "-" ]
		then
			echo "${chr} ${start} ${stop} ${sign}" | tr ' ' '\t' >> inversions.bed
		fi
		if [ -z ${start2} ]
		then
			start2=${start}
			stop2=${stop}
			sign2=${sign}
		else
			if [ ${sign} == ${sign2} ]
			then
				start2=${start}
				stop2=${stop}
				sign2=${sign}
			else
				if [[ ${stop2} -gt ${start} ]]
				then
					echo "${chr} ${start2} ${start}" | \
					tr ' ' '\t' >> ${species}_${genotype}_complex_rearrangments.bed
				else
					echo "${chr} ${stop2} ${start}" | \
					tr ' ' '\t' >> ${species}_${genotype}_breakpoints.bed
				fi
				start2=${start}
				stop2=${stop}
				sign2=${sign}
			fi
		fi
	done
done

#Make the breakpoints bed file for the reference 
zcat ${genespace} | awk -v a=${ref/*_/} -v b=${genotype} '$1==a && $2==b' | cut -f4 | sort | uniq | while read seq
do
	chr=$(echo ${seq} | cut -d ' ' -f4)
	start2=
	stop2=
	sign2=
	zcat ${genespace} | awk -v a=${ref/*_/} -v b=${genotype} -v c=${chr} '$1==a && $2==b && $3==c && $4==c' | \
	cut -f3,6,7,20 | sort -k2n | while read line
	do
		start=$(echo ${line} | cut -d ' ' -f2)
		stop=$(echo ${line} | cut -d ' ' -f3)
		sign=$(echo ${line} | cut -d ' ' -f4)
		if [ -z ${start2} ]
		then
			start2=${start}
			stop2=${stop}
			sign2=${sign}
		else
			if [ ${sign} == ${sign2} ]
			then
				start2=${start}
				stop2=${stop}
				sign2=${sign}
			else
				if [[ ${stop2} -gt ${start} ]]
				then
					echo "${chr} ${start2} ${start}" | tr ' ' '\t' >> ${ref}_complex_rearrangments.bed
				else
					echo "${chr} ${stop2} ${start}" | tr ' ' '\t' >> ${ref}_breakpoints.bed
				fi
				start2=${start}
				stop2=${stop}
				sign2=${sign}
			fi
		fi
	done
done


#Set the input files
version=$(ls ${path2}/${species}/${genotype}/ref/${genotype}-v*.fa | sed s/.*\-v// | sed s/\.fa//)
fasta="${path2}/${species}/${genotype}/ref/${genotype}-v${version}.fa"
if [[ ! -f ${fai} ]]
then
	samtools index ${fasta}
fi
ref_version=$(ls ${path2}/${ref/_*/}/${ref/*_/}/ref/${ref/*_/}-v*.fa | sed s/.*\-v// | sed s/\.fa//)
ref_fa="${path2}/${ref/_*/}/${ref/*_/}/ref/${ref/*_/}-v${ref_version}.fa"
<<<<<<< HEAD
ref_fai="${path2}/${species}/${genotype}/ref/${genotype}-v${version}.fa.fai"
=======
ref_fai="${path2}/${ref/_*/}/${ref/*_/}/ref/${ref/*_/}-v${ref_version}.fa.fai"
if [[ ! -f ${ref_fai} ]]
then
	samtools index ${ref_fasta}
fi
>>>>>>> 8b58413c50a5cfb19aa0f8321a88e3fc6a01a013

#Make & align simulated reads from breakpoints to ref genome
#Loop over each region in the breakpoints.bed
cat ${ref}_breakpoints.bed | while read line
do
	#get the region to analyze
	region=$(echo ${line} | tr ' ' '_')
	#Create an index for that chromosome/sequence only to reduce multimapping to other chroms
	#If index already exists, skip
	if [ ! -d indexes/${region/_*/} ]
	then
		mkdir indexes/${region/_*/}
		samtools faidx ${fasta} ${region/_*/} > indexes/${region/_*/}/${region/_*/}.fa
		bwa index -p indexes/${region/_*/}/${region/_*/} indexes/${region/_*/}/${region/_*/}.fa
	fi
	#Create a directory for that region
	dir="regions/${region}"
	mkdir regions/${region}
	#Create the bed for that region
	echo ${line} | tr ' ' '\t' > ${dir}/${region}.bed
	#Extend the region by the amount in the extend variabl
	bedtools slop \
		-b ${extend} \
		-i ${dir}/${region}.bed \
		-g ${ref_fai} > ${dir}/${region}_ext.bed
	#Extract the fasta sequence
	bedtools getfasta \
		-bed ${dir}/${region}_ext.bed \
<<<<<<< HEAD
		-fi ${ref_fasta} \
=======
		-fi ${ref_fa} \
>>>>>>> 8b58413c50a5cfb19aa0f8321a88e3fc6a01a013
		-fo ${dir}/${region}_ext.fa
	#Make simulated reads using wgsim
	wgsim \
		-e 0 \
		-d 400 \
		-s 50 \
		-N 1000000 \
		-1 150 \
		-2 150 \
		-r 0 \
		-R 0 \
		-X 0 \
		-S 0 \
		${dir}/${region}_ext.fa ${dir}/${region}_1.fastq ${dir}/${region}_2.fastq
	#Align the simulated reads
	bwa mem \
		-t ${threads} \
		-R "@RG\tID:${region}\tLB:SIM\tPL:SIM\tSM:${region}\tPU:SIM" \
		-M indexes/${region/_*/}/${region/_*/} ${dir}/${region}_1.fastq ${dir}/${region}_2.fastq | \
		samtools view -@ ${threads2} -bSh | \
		samtools sort -@ ${threads2} > ${dir}/${region}.bam
done
#Merge the mapped reads for all the regions
samtools merge ${species}_${genotype}_breakpoints.bam regions/*/*bam
#Index the reads
samtools index ${species}_${genotype}_breakpoints.bam

echo "Done"
