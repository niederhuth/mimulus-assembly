#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=100GB
#SBATCH --job-name LTR_retriever
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
fasta= #genome fasta, if left blank will look in the current directory
run_LTR_harvest=TRUE #TRUE/FALSE, run LTR harvest to find LTRs, will ignore inharvest and EDTA
inharvest= #results from previous LTR_harvest run, ignored if run_LTR_harvest=TRUE, mutually exclusive to EDTA
EDTA="edta" #directory for EDTA results, ignored if run_LTR_harvest=TRUE, mutually exclusive to inharvest
threads=40

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
path1=$(pwd | sed s/data.*/scripts/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*/${species}\\/${genotype}\\/// | sed s/\\/.*//)
version=$(ls ${sample}-v*.fa | sed s/.*\-v// | sed s/.fa//)
path2=$(pwd)
path3="LTR_retriever"

#Look for fasta file, there can only be one!
if [ -z ${fasta} ]
then
	echo "No fasta fasta provided, looking for fasta"
	if ls *.fa >/dev/null 2>&1
	then
		fasta_name=$(ls *fa | sed s/.*\ //)
		fasta=${path2}/$(ls *fa | sed s/.*\ //)
		echo "Fasta file ${fasta} found"
	elif ls *.fasta >/dev/null 2>&1
	then
		fasta_name=$(ls *fasta | sed s/.*\ //)
		fasta=${path2}/$(ls *fasta | sed s/.*\ //)
		echo "Fasta file ${fasta} found"
	elif ls *.fna >/dev/null 2>&1
	then
		fasta_name=$(ls *fna | sed s/.*\ //)
		fasta=${path2}/$(ls *fna | sed s/.*\ //)
		echo "Fasta file ${fasta} found"
	else
		echo "No fasta file found, please check and restart"
	fi
else
	echo "input fasta: ${fasta}"
fi

#Set path EDTA results
if [ ${run_LTR_harvest} = FALSE ]
then
	if [ ${inharvest} & -z ${EDTA} ]
	then
		if [ -f ${inharvest} ]
		then
			echo "LTR Harvest output scn found."
			scn=${inharvest}
		else
			echo "LTR Harvest output scn not found."
			echo "Please check that all paths are correct and restart"
		fi
	elif [ -z ${inharvest} & ${EDTA} ]
	then
		if [ -d ${EDTA} ]
		then
			echo "EDTA output directory found"
			echo "Checking for LTR Harvest output scn file"
			if [ -f ${EDTA}/${fasta_name}.mod.EDTA.raw/LTR/${fasta_name}.mod.rawLTR.scn ]
			then
				echo "LTR Harvest output scn file found"
				scn="${path2}/${EDTA}/${fasta_name}.mod.EDTA.raw/LTR/${fasta_name}.mod.rawLTR.scn"
			else
				echo "EDTA LTR Harvest output scn file not found"
				echo "Please check that all paths are correct and restart"
			fi
		else
			echo "EDTA output directory not found"
			echo "Please check that all paths are correct and restart"
		fi
	else
		echo "No input file found and run_LTR_harvest = FALSE"
		echo "inharvest & EDTA are mutually exclusive."
		echo "You must provide either an LTR Harvest output scn file or the path to a previous EDTA run."
		echo "Check your settings and restart."
	fi
fi

#Make & cd to output directory
if [ -d ${path3} ]
then
	cd ${path3}
	cp ${fasta} ./
	fasta=${fasta_name}
else
	mkdir ${path3}
	cd ${path3}
	cp ${fasta} ./
	fasta=${fasta_name}
fi

#Run LTR Harvest if run_LTR_harvest=TRUE
if [ ${run_LTR_harvest} = TRUE ]
then
	#Run gt suffixerator
	echo "Running genometools suffixerator"
	gt suffixerator \
		-db ${fasta} \
		-indexname ${fasta} \
		-tis -suf -lcp -des -ssp -sds -dna
	#Run LTR harvest
	echo "Running LTR Harvest"
	gt ltrharvest \
		-index ${fasta} \
		-minlenltr 100 \
		-maxlenltr 7000 \
		-mintsd 4 \
		-maxtsd 6 \
		-motif TGCA \
		-motifmis 1 \
		-similar 85 \
		-vic 10 \
		-seed 20 \
		-seqids yes > ${fasta}.harvest.scn
	#Run LTR_Finder
	echo "Running LTR_Finder"
	LTR_FINDER_parallel \
		-seq ${fasta} \
		-threads ${threads} \
		-harvest_out \
		-size 1000000 \
		-time 300
	#Combine results into rawLTR.scn
	cat ${fasta}.harvest.scn ${fasta}.finder.combine.scn > ${fasta}.rawLTR.scn
	#Set scn
	scn=${fasta}.rawLTR.scn
fi

#Run LTR_retriever using EDTA results
echo "Running LTR_retriever"
${path1}/annotation/LTR_retriever/LTR_retriever \
	-threads ${threads}	\
	-genome ${fasta} \
	-inharvest ${scn}

echo "Done"
