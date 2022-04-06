#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=50GB
#SBATCH --job-name filter_genes
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=50
maker_dir="../maker_round2"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/maker/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/maker/lib:${LD_LIBRARY_PATH}"

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
path3="gene_filtering"

#Look for fasta file, there can only be one!
if [ -z ${fasta} ]
then
	echo "No input fasta provided, looking for fasta"
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
else
	echo "Input fasta: ${fasta}"
fi

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

#Set various datasets
transcripts=${maker_dir}/maker_transcripts.fa
proteins=${maker_dir}/maker_proteins.fa
gff=${maker_dir}/maker.gff
gene_ids=maker_standard_gene_list.txt

#Check for Pfam A
echo "Downloading Pfam-A"
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz

#Prepare Pfam A hmm-database
echo "Preparing Pfam A hmm database"
hmmpress Pfam-A.hmm

#Search Pfam A hmm domains
echo "Searching against Pfam A hmm profiles"
hmmscan \
	--domE 1e-5 \
	-E 1e-5 \
	--cpu ${threads} \
	-o pfam_alignments.out \
	--tblout prot_domains.out \
	Pfam-A.hmm \
	${proteins}

#Generate maker standard gene list
echo "Generating Maker Standard Gene List"
perl ${path2}/annotation/generate_maker_standard_gene_list.pl \
	--input_gff ${gff} \
	--pfam_results prot_domains.out \
	--pfam_cutoff 1e-10 \
	--output_file maker_standard_gene_list.txt

#Make Maker Standard Files
echo "Making Maker Standard Transcripts fasta"
perl ${path2}/annotation/get_subset_of_fastas.pl \
	-l ${gene_ids} \
	-f ${transcripts} \
	-o maker_standard_transcripts.fa
	
echo "Making Maker Standard Protein fasta"
perl ${path2}/annotation/get_subset_of_fastas.pl \
	-l ${gene_ids} \
	-f ${proteins} \
	-o maker_standard_proteins.fa

echo "Making Maker Standard GFF file"
perl ${path2}/annotation/create_maker_standard_gff.pl \
	--input_gff ${gff} \
	--output_gff maker_standard.gff \
	--maker_standard_gene_list ${gene_ids}

#Download and make Transposase blast DB
echo "Downloading Tpases020812"
wget http://www.hrt.msu.edu/uploads/535/78637/Tpases020812.gz
gunzip Tpases020812.gz
#remove trailing white space in the Tpases020812 file
sed -i 's/[ \t]*$//' Tpases020812

echo "Making Transposase diamond blast DB"
diamond makedb \
	--threads ${threads} \
	--in Tpases020812 \
	--db Tpases020812.dmnd

#Run diamond blastp against Transposases 
echo "Running diamond blastp on Transposases"
diamond blastp \
	--threads ${threads} \
	--db Tpases020812.dmnd \
	--query ${proteins} \
	--out TE_blast.out\
	--evalue 1e-10 \
	--outfmt 6

#Download Gypsy DB hmm files and format the hmm database
echo "Downloading GyDB_collection"
wget https://gydb.org/extensions/Collection/collection/db/GyDB_collection.zip
unzip GyDB_collection.zip
echo "Combining GyDB hmm profiles"
cat GyDB_collection/profiles/*hmm > all_gypsy.hmm
echo "Formatting all_gypsy.hmm database"
hmmpress all_gypsy.hmm

#Search gypsy hmm profiles
echo "Searching against gypsy hmm profiles"
hmmscan \
	--domE 1e-5 \
	-E 1e-5 \
	--cpu ${threads} \
	-o gypsy_alignments.out \
	--tblout gypsyHMM_analysis.out \
	all_gypsy.hmm \
	${proteins}

#Create a genelist with no TEs
echo "Creating gene list with TEs removed"
python ${path2}/annotation/create_no_TE_genelist.py \
	--input_file_TEpfam ${path1}/annotation/TE_Pfam_domains.txt \
	--input_file_maxPfam prot_domains.out \
	--input_file_geneList_toKeep maker_standard_gene_list.txt \
	--input_file_TEhmm gypsyHMM_analysis.out \
	--input_file_TEblast TE_blast.out \
	--output_file noTE_maker_standard_gene_list.txt

#Generate noTE files
echo "Making no TE Maker Standard Transcripts fasta"
perl ${path2}/annotation/get_subset_of_fastas.pl \
	-l noTE_maker_standard_gene_list.txt \
	-f maker_standard_transcripts.fa \
	-o ${fasta/.fa/}_maker_standard_transcripts_noTE.fa

echo "Making no TE Maker Standard proteins fasta"
perl ${path2}/annotation/get_subset_of_fastas.pl \
	-l noTE_maker_standard_gene_list.txt \
	-f maker_standard_proteins.fa \
	-o ${fasta/.fa/}_maker_standard_proteins_noTE.fa

echo "Making no TE Maker Standard gff"
perl ${path2}/annotation/create_maker_standard_gff.pl \
	--input_gff maker_standard.gff \
	--output_gff ${fasta/.fa/}_maker_standard_transcripts_noTE.gff \
	--maker_standard_gene_list noTE_maker_standard_gene_list.txt

echo "Done"

