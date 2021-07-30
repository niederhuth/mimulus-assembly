#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=400GB
#SBATCH --job-name edta
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=40
protein="Mguttatus/IM62/ref/annotations/IM62-v2-protein.fa"

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
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${genotype}\\/// | sed s/\\/.*//)
version=$(ls ${genotype}-v*.fa | sed s/.*\-v// | sed s/.fa//) 
path2=""

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

#Run EDTA
if [ something ]
then
	exonerate \
		--cores ${threads} \
		--model protein2genome \
		--bestn 5 \
		--minintron 10 \
		--maxintron 3000 \
		--query ${path1}/${protein} \
		--target ${fasta} \
		--showtargetgff yes \
		--showalignment no \
		--ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" \
		--targetchunkid 1 \
		--targetchunktotal 20 > ${output}

		--refine

	echo "Converting gtf to gff3"
	gffread -E protein_exonerate.gtf -o- > protein_exonerate.gff
fi


echo "Done"
