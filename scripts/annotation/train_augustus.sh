#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200GB
#SBATCH --job-name train_augustus
#SBATCH --output=%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
fasta= #input fasta, if left blank, will look for it in current directory
input_gff="../maker_round1/maker_round1.gff"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/maker/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/maker/lib:${LD_LIBRARY_PATH}"
#Export path to agusutus config files
#export ZOE="${conda}/envs/maker" #Need to check
export AUGUSTUS_CONFIG_PATH="${conda}/envs/maker/config/"
#export REPEATMASKER_LIB_DIR=
#export REPEATMASKER_MATRICES_DIR=

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
path3="augustus_training"
path4="../maker_round1"

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

#maker2zff
#For determining which genes are High Confidence for Retraining, there are 7 criteria.
#-c fraction  The fraction of splice sites confirmed by an EST alignment, default 0.5
#-e fraction  The fraction of exons that overlap an EST alignment, default 0.5
#-o fraction  The fraction of exons that overlap any evidence (EST or Protein), default 0.5
#-a fraction  The fraction of splice sites confirmed by an ab-initio prediction, default 0
#-t fraction  The fraction of exons that overlap an ab-initio prediction, default 0
#-l number    The min length of the protein sequence produced by the mRNA
#-x number    Max AED to allow 0.5 is default
#-n           No filtering.  Accept all.
echo ""
echo "maker2zff -x 0.2 -l 200 $MAKER_GFF_FILE_W_FASTA"
maker2zff \
	-c 0.5 \
	-e 0.5 \
	-o 0.5 \
	-a 0 \
	-t 0 \
	-l 200 \
	-x 0.2 \
	${input_gff}

#fathom
#-validate [-quiet]
#-gene-stats [-errors-ok -nucleotide -dinucleotide]
#-categorize <int>
#-export <int> [-plus -errors-ok]
#-extract <feature> -length <int> -offset <int>
#-exon-intron
#-split <-number <int> | -training <float> | -GC <float> | -repeat <float>>
#-ace-format <-gene-method <string> [-dna -extra <string>]>
#-compare-genes <predictions> [-details]
#-score-genes <hmm> [-errors-ok]
#-filter-genes <hmm> -min-score <float> -min-length <int>
echo ""
echo "fathom genome.ann genome.dna -categorize 1000"
fathom \
	genome.ann \
	genome.dna \
	-categorize 1000

NUMFOUND="`grep -c x'>' uni.ann`"

if [ ${NUMFOUND} -gt 499 ]
then
	NUMFOUND=500
fi

TEMPSPLIT=$((NUMFOUND/2))
NUMSPLIT=${TEMPSPLIT/.*}

echo "number found after fathom: $NUMFOUND"
echo "number after split: $NUMSPLIT"

#Convert the uni.ann and uni.dna output from fathom into a genbank formatted file.
#fathom_to_genbank.pl is available on github: https://github.com/Childs-Lab/GC_specific_MAKER.
#Modify the path here, or ensure that fathom_to_genbank.pl is in your $PATH.
echo "fathom_to_genbank.pl --annotation_file uni.ann --dna_file uni.dna --genbank_file augustus.gb --number ${NUMFOUND}"
fathom_to_genbank.pl \
	--annotation_file uni.ann \
	--dna_file uni.dna \
	--genbank_file augustus.gb \
	--number ${NUMFOUND}

#To get the subset of fastas that correspond to the genes in the genbank file, follow these steps.
#get_subset_of_fastas.pl is available on github: https://github.com/Childs-Lab/GC_specific_MAKER.
#Modify the path here, or ensure that get_subset_of_fastas.pl is in your $PATH.
perl -e  'while (my $line = <>){ if ($line =~ /^LOCUS\s+(\S+)/) { print "$1\n"; } }' ${WORKING_DIR}/augustus.gb > ${WORKING_DIR}/genbank_gene_list.txt

echo "get_subset_of_fastas.pl -l ${WORKING_DIR}/genbank_gene_list.txt -f ${WORKING_DIR}/uni.dna -o ${WORKING_DIR}/genbank_gene_seqs.fasta"
perl ${path2}/annotation/get_subset_of_fastas.pl \
	-l  ${WORKING_DIR}/genbank_gene_list.txt \
	-f ${WORKING_DIR}/uni.dna \
	-o ${WORKING_DIR}/genbank_gene_seqs.fasta

#Split the known genes into test and training files.
#randomSplit.pl is a script provided by AUGUSTUS.
#Modify the path here, or ensure that randomSplit.pl is in your $PATH.
echo "randomSplit.pl ${WORKING_DIR}/augustus.gb ${NUMSPLIT}"
randomSplit.pl ${WORKING_DIR}/augustus.gb ${NUMSPLIT}

#We will use autoAug.pl for training because we have transcript alignments that will be used as hints.
#The etraining and optimize_augustus.pl scripts from AUGUSTUS will not be used as they do not allow for the use of hints.
#Run autoAug.pl.
#This script is provided in the AUGUSTUS installation.
#Modify the path here, or ensure that autoAug.pl is in your $PATH.
echo "autoAug.pl --species=$AUGUSTUS_SPECIES_NAME --genome=${WORKING_DIR}/genbank_gene_seqs.fasta --trainingset=${WORKING_DIR}/augustus.gb.train --cdna=$CDNA_FASTA --noutr"
autoAug.pl \
	--species=${AUGUSTUS_SPECIES_NAME} \
	--genome=${WORKING_DIR}/genbank_gene_seqs.fasta \
	--trainingset=${WORKING_DIR}/augustus.gb.train \
	--cdna=$CDNA_FASTA \
	--noutr

#Run the batch scripts generated by the previous command
cd "${WORKING_DIR}/autoAug/autoAugPred_abinitio/shells"
cd ${WORKING_DIR}/autoAug/autoAugPred_abinitio/shells

# The number of ./aug# scripts is variable.  Run until the new file does not exist.
x=1
while [ -e ./aug${x} ]
do
    echo "A.  $x"
    ./aug${x}
    let x=x+1
done

#Run the next command as indicated by autoAug.pl in step 5.
#When above jobs are finished, continue by running the command autoAug.pl
echo "cd $WORKING_DIR"
cd $WORKING_DIR
echo "autoAug.pl --species=$AUGUSTUS_SPECIES_NAME --genome=${WORKING_DIR}/genbank_gene_seqs.fasta --useexisting --hints=${WORKING_DIR}/autoAug/hints/hints.E.gff  -v -v -v  --index=1"
autoAug.pl \
	--species=${AUGUSTUS_SPECIES_NAME} \
	--genome=${WORKING_DIR}/genbank_gene_seqs.fasta \
	--useexisting \
	--hints=${WORKING_DIR}/autoAug/hints/hints.E.gff  \
	-v -v -v --index=1

#Run the batch scripts as indicated by autoAug.pl in step 7.
echo "cd ${WORKING_DIR}/autoAug/autoAugPred_hints/shells/"
cd ${WORKING_DIR}/autoAug/autoAugPred_hints/shells/

#The number of ./aug# scripts is variable.  Run until the new file does not exist.
let x=1
while [ -e ./aug${x} ]
do
    echo "B.  $x"
    ./aug${x}
    let x=x+1
done

#Checked sensitivity and specificity of the newly trained AUGUSTUS HMM by using the test data.
cd ${WORKING_DIR}
echo "augustus --species=$AUGUSTUS_SPECIES_NAME augustus.gb.test"
augustus --species=$AUGUSTUS_SPECIES_NAME augustus.gb.test

echo "Done"

