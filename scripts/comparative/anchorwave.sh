#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=400GB
#SBATCH --job-name anchorwave
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
mode="genoAli" #genoAli or proali
ref_fasta= #reference genome, if left blank, will look for in the ref directory for that genotype
ref_gff= #reference gff, if left blank, will look for in the ref directory for that genotype
threads=10
ref_coverage=1 #For proali use only, set max coverage for reference genome
query_coverage=1 #For proali use only, set max coverage for query genome
seqs="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/wga/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/wga/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="genome_alignment"
datatype="genome"
path2=$(pwd | sed s/data.*/data/)
path3="anchorwave"

#Check for and make/cd working directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#find the reference fasta if not provided
if [ -z ${ref_fasta} ]
then
	#Get version number for genome fasta
	genome_ver=$(ls ${path2}/${species}/${genotype}/ref/${genotype}-v*.fa | head -1 | \
		sed s/.*\-v// | sed s/\.fa//)
	#Set the ref_fasta path
	ref_fasta="${path2}/${species}/${genotype}/ref/${genotype}-v${genome_ver}.fa"
	if [[ ! -z ${seqs} ]]
	then
		echo ${seqs} | tr ' ' '\n' > regions
		samtools faidx ${ref_fasta} regions > ${genotype}-v${genome_ver}.fa
		ref_fasta="$(pwd)/${genotype}-v${genome_ver}.fa"
	fi
fi 
#Find the reference gff if not provided
if [ -z ${ref_gff} ]
then
	#Get the version number for the gff file
	anno_ver=$(ls ${path2}/${species}/${genotype}/ref/annotations/${genotype}-v*.gff | head -1 | \
		sed s/.*\-v// | sed s/\.gff//)
	#Set the ref_gff path
	ref_gff="${path2}/${species}/${genotype}/ref/annotations/${genotype}-v${anno_ver}.gff"
	if [[ ! -z ${seqs} ]]
	then
		fgrep -f regions ${ref_gff} > ${genotype}-v${genome_ver}.gff
		ref_gff="$(pwd)/${genotype}-v${genome_ver}.gff"
	fi
fi

#Get list of genomes
genomes=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${condition} \
	-v e=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $4 == d && $5 == e) print $7}' \
	${path1}/samples.csv)

#Loop over genomes and run anchorwave
for i in ${genomes}
do
	mkdir ${i}_alignment
	cd ${i}_alignment
	#Set path to query genome
	query_ver=$(ls ${path2}/$(echo ${i} | sed s/_/\\//)/ref/${i/*_/}-v*.fa | head -1 | \
		sed s/.*\-v// | sed s/\.fa//)
	query_fasta="${path2}/$(echo ${i} | sed s/_/\\//)/ref/${i/*_/}-v${query_ver}.fa"
	if [[ ! -z ${seqs} ]]]
	then
		samtools faidx ${query_fasta} ../regions > ${i/*_/}-v${query_ver}.fa
		query_fasta="$(pwd)/${i/*_/}-v${query_ver}.fa"
	fi
	#Liftover annotations
	echo "Extract the CDS sequences"
	anchorwave gff2seq \
		-i ${ref_gff} \
		-r ${ref_fasta} \
		-o cds.fa
	#Map CDS sequences to reference genome
	echo "Aligning ${species}_${genotype} CDS to itself with minimap2"
	minimap2 \
		-x splice \
		-t ${threads} \
		-k 12 \
		-a \
		-p 0.4 \
		-N 20 \
		${ref_fasta} cds.fa > ${species}_${genotype}.sam
	#Map CDS sequences to query genome
	echo "Aligning ${species}_${genotype} CDS to ${i} with minimap2"
	minimap2 \
		-x splice \
		-t ${threads} \
		-k 12 \
		-a \
		-p 0.4 \
		-N 20 \
		${query_fasta} cds.fa > ${i}.sam
	#Run Genome Alignment
	if [ ${mode} = "genoAli" ]
	then
		#Run genoAli
		echo "Aligning ${i} to ${species}_${genotype} with genoAli"
		anchorwave genoAli \
			-t ${threads} \
			-IV \
			-i ${ref_gff} \
			-r ${ref_fasta} \
			-as cds.fa \
			-a ${i}.sam \
			-ar ${species}_${genotype}.sam \
			-s ${query_fasta} \
			-n ${i}.anchors \
			-o ${i}.maf \
			-f ${i}.fragmentation.maf
	elif [ ${mode} = "proali" ]
	then
		#Run proali
		echo "Aligning ${i} to ${species}_${genotype} with proali"
		anchorwave proali \
			-t ${threads} \
			-R ${ref_coverage} \
			-Q ${query_coverage} \
			-i ${ref_gff} \
			-r ${ref_fasta} \
			-as cds.fa \
			-a ${i}.sam \
			-ar ${species}_${genotype}.sam \
			-s ${query_fasta} \
			-n ${i}.anchors \
			-o ${i}.maf \
			-f ${i}.fragmentation.maf
	fi
	#Change out of that directory
	cd ../
done

echo "Done"
