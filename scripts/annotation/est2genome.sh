#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=400GB
#SBATCH --job-name est2genome
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=40
minintron=10
maxintron=3000
bestn=5
ryo="no" #">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/EDTA/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/EDTA/lib:${LD_LIBRARY_PATH}"

#Set temporary directories for large memory operations
export TMPDIR=$(pwd | sed s/data.*/data/)
export TMP=$(pwd | sed s/data.*/data/)
export TEMP=$(pwd | sed s/data.*/data/)

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
path2=$(pwd | sed s/data.*/data/)
path3=$(pwd | sed s/data.*/scripts/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=${genotype}
version=$(ls ${genotype}-v*.fa | sed s/.*\-v// | sed s/.fa//) 
path4="exonerate"

#Look for fasta file, there can only be one!
if ls *.fa >/dev/null 2>&1
then
	fasta=$(ls *fa | sed s/.*\ //)
	echo "Fasta file ${fasta} found"
elif ls *.fasta >/dev/null 2>&1
then
	fasta=$(ls *fasta | sed s/.*\ //)
	echo "Fasta file ${fasta} found"
elif ls *.fna >/dev/null 2>&1
then
	fasta=$(ls *fna | sed s/.*\ //)
	echo "Fasta file ${fasta} found"
else
	echo "No fasta file found, please check and restart"
fi

#Make output directory
if [ -d ${path4} ]
then
	cd ${path4}
else
	mkdir ${path4}
	cd ${path4}
fi

#Get transcript sources
transcripts=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	'{if ($1 == a && $2 == b && $3 == c) print $5}' \
	${path1}/annotation/annotation_sources.csv)

#Run EDTA
for i in ${transcripts}
do
	source_species=$(echo $i | cut -d '_' -f1)
	source_genotype=$(echo $i | cut -d '_' -f2)
	input="${path2}/${source_species}/${source_genotype}/ref/annotations/${source_genotype}-v*-transcript.fa"
	output1="${i}_transcript.output"
	output2="${i}_transcript.gff"
	if [ -s ${output2} ]
	then
		echo "${output2} already exists."
		echo "To rerun this step, delete ${output2} and resubmit"
	else
		echo "Running exonerate est2genome on ${fasta}"
		exonerate \
			--cores ${threads} \
			--model est2genome \
			--bestn ${bestn} \
			--minintron ${minintron} \
			--maxintron ${maxintron} \
			--query ${input} \
			--target ../${fasta} \
			--showtargetgff yes \
			--showalignment no \
			--showvulgar no \
			--ryo ${ryo} \
			--refine > ${output1}
		#Reformat exonerate output
		echo "Reformatting exonerate gff"
		perl ${path3}/annotation/reformat_exonerate_transcript_gff.pl \
			--input_gff ${output1} \
			--output_gff ${output2}
	fi
done

echo "Done"
