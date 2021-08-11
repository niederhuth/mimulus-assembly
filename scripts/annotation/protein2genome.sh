#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=400GB
#SBATCH --job-name protein2genome
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Protein source
protein="Mguttatus/IM62/ref/annotations/IM62-v2-protein.fa"

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
export TMPDIR=$(pwd | sed s/cre-genomes.*//)
export TMP=$(pwd | sed s/cre-genomes.*//)
export TEMP=$(pwd | sed s/cre-genomes.*//)

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/data/)
path2=$(pwd | sed s/data.*/scripts/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${genotype}\\/// | sed s/\\/.*//)
version=$(ls ${genotype}-v*.fa | sed s/.*\-v// | sed s/.fa//) 
path3="exonerate"
output1="protein2genome.output"
output2="protein2genome.gff"

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
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#Run EDTA
if [ -s ${output2} ]
then
	echo "${output2} already exists."
	echo "To rerun this step, delete ${output2} and resubmit"
else
	echo "Running exonerate protein2genome on ${fasta}"
	exonerate \
		--cores ${threads} \
		--model protein2genome \
		--bestn ${bestn} \
		--minintron ${minintron} \
		--maxintron ${maxintron} \
		--query ${path1}/${protein} \
		--target ../${fasta} \
		--showtargetgff yes \
		--showalignment no \
		--showvulgar no \
		--ryo ${ryo} \
		--refine > ${output1}
	#Reformat exonerate output
	echo "Reformatting exonerate gff"
	perl ${path2}/annotation/reformat_exonerate_protein_gff.pl \
		--input_gff ${output1} \
		--output_gff ${output2}
fi

echo "Done"