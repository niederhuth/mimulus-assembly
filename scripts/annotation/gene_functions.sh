#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=50GB
#SBATCH --job-name=gene_functions
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=50
protein_fasta="../maker_round2/*.all.maker.proteins.fasta"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/gene-function/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/gene-function/lib:${LD_LIBRARY_PATH}"

#Set temporary directories for large memory operations
export TMPDIR=$(pwd)
export TMP=$(pwd)
export TEMP=$(pwd)

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
path2=$(pwd | sed s/data.*/scripts/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*${species}\\/${genotype}\\/// | sed s/\\/.*//)
path3="gene_functions"

#Make & cd to directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#Set temporary directories for large memory operations
export TMPDIR=$(pwd)
export TMP=$(pwd)
export TEMP=$(pwd)

#Set some more variables
proteins=$(ls ${protein_fasta} | sed s/.*\ //)
output=$(echo ${proteins} | sed s/.*\\/// | sed s/\\..*//)
echo ${proteins}

#Run interproscan
echo "Running interproscan"
${path2}/annotation/interproscan/interproscan.sh \
	-cpu ${threads} \
	-appl pfam \
	-goterms \
	-pa \
	-dp \
	-iprlookup \
	-t p \
	-f TSV \
	-i ${proteins} \
	-o ${output}.iprscan

#Download Arabidopsis genes and create diamond DB
echo "Downloading Arabidopsis TAIR10 proteins"
wget https://www.arabidopsis.org/download_files/Proteins/TAIR10_protein_lists/TAIR10_pep_20110103_representative_gene_model
echo "Making diamond blast DB for "
diamond makedb \
	--threads ${threads} \
	--in TAIR10_pep_20110103_representative_gene_model \
	--db TAIR10.dmnd

#Run diamond blastp against Arabidopsis 
echo "Running diamond blastp on "
diamond blastp \
	--threads ${threads} \
	--db TAIR10.dmnd \
	--query ${proteins} \
	--out ${output}_TAIR10_blast.out \
	--evalue 1e-6 \
	--max-hsps 1 \
	--max-target-seqs 5 \
	--outfmt 0

#Download and format Arabidopsis TAIR10 functional descriptions
echo "Downloading and formatting Arabidopsis TAIR10 functional descriptions"
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_functional_descriptions
perl -e  'while (my $line = <>){ my @elems = split "\t", $line; if($elems[2] ne "") {print "$elems[0]\t$elems[2]\n"}}' \
TAIR10_functional_descriptions > TAIR10_short_functional_descriptions.txt

#Format annotations 
echo "Creating functional annotation file"
perl ${path2}/annotation/pl/create_functional_annotation_file.pl \
	--protein_fasta ${proteins} \
	--model_annot TAIR10_short_functional_descriptions.txt \
	--model_blast ${output}_TAIR10_blast.out \
	--pfam_results_file ${output}.iprscan \
	--max_hits 5 \
	--output prot_func_desc_list.txt

echo "Done"
