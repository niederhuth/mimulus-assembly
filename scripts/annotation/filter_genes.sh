#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200GB
#SBATCH --job-name filter_genes
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20

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
path3=""

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
if [[ -z Pfam-A.hmm || -z Pfam-A.hmm.gz ]]
then
	echo Pfam-A.hmm found
	if [ -z Pfam-A.hmm.gz ]
	then
		gunzip Pfam-A.hmm.gz
	fi
else
	echo "Downloading Pfam-A"
	wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
	gunzip Pfam-A.hmm.gz
fi

#Prepare Pfam A hmm-database
echo "Preparing Pfam A hmm database"
hmmpress Pfam-A.hmm

#Search Pfam A hmm domains
echo "Searching against Pfam A hmm profiles"
hmmscan \
	--domE 1e-5 \
	-E 1e-5 \
	--cpu ${threads} \
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
	--output_gff maker_standard.gff3 \
	--maker_standard_gene_list ${gene_ids}

#Make Transposase BLAST DB
echo "Making Transposase BLAST DB"
makeblastdb \
	-in Tpases020812 \
	-title "Tpases020812DB" \
	-parse_seqids \
	-out Tpases020812DB \
	-dbtype prot


#Run BLAST against Transposases 
echo "Running BLAST on Transposases"
blastp \
	-db Tpases020812DB \
	-query ${proteins} \
	-out blast.out \
	-evalue 1e-10 \
	-outfmt 6 \
	-num_threads ${threads}

#Download Gypsy DB hmm files and format the hmm database
if [ -z all_gypsy.hmm ]
then
	echo "all_gypsy.hmm found"
else
	echo "Downloading GyDB_collection"
	wget https://gydb.org/extensions/Collection/collection/db/GyDB_collection.zip
	unzip GyDB_collection.zip
	echo "Combining GyDB hmm profiles"
	cat GyDB_collection/profiles/*hmm > all_gypsy.hmm
fi
echo "Formatting all_gypsy.hmm database"
hmmpress all_gypsy.hmm

#Search gypsy hmm profiles
echo "Searching against gypsy hmm profiles"
hmmscan \
	--domE 1e-5 \
	-E 1e-5 \
	--cpu ${threads} \
	--tblout gypsyHMM_analysis.out \
	all_gypsy.hmm \
	${proteins}

#
python ${path2}/annotation/create_no_TE_genelist.py \
	--input_file_TEpfam TE_Pfam_domains.txt \
	--input_file_maxPfam Sreb_prot_domains.out \
	--input_file_geneList_toKeep Sreb_maker2_standard_gene_list.txt \
	--input_file_TEhmm gypsyHMM_analysis.out \
	--input_file_TEblast Sreb_blast.out \
	--output_file noTE_Sreb_v1.0_maker2standard_gene_list.txt

#
gene_ids=noTE_Sreb_v1.0_maker2standard_gene_list.txt
transcripts=Sreb_maker_standard_transcripts.fasta
proteins=Sreb_maker_standard_proteins.fasta
input_gff=Sred_maker2_standard.gff3

perl ${path2}/annotation/get_subset_of_fastas.pl \
	-l $gene_ids \
	-f $transcripts \
	-o Sreb_v1.0_maker_std_transcripts_noTE.fasta

perl ${path2}/annotation/get_subset_of_fastas.pl \
	-l $gene_ids \
	-f $proteins \
	-o Sreb_v1.0_maker_std_proteins_noTE.fasta

perl ${path2}/annotation/create_maker_standard_gff.pl \
	--input_gff $input_gff \
	--output_gff Sreb_v1.0_maker_std_noTE.gff3 \
	--maker_standard_gene_list noTE_Sreb_v1.0_maker2standard_gene_list.txt

echo "Done"

