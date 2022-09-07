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
inharvest= #if left blank, will look for edta output in current directory
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

#If LTR Harvest not provided, look for EDTA output
if [ -z ${inharvest} ]
then
	echo "No LTR Harvest file provided, looking for EDTA output directory"
	if [ -d edta ]
	then
		echo "EDTA output directory found"
		echo "Checking for LTR Harvest output scn file"
		if [ -f edta/${fasta_name}.mod.EDTA.raw/LTR/${fasta_name}.mod.rawLTR.scn ]
		then
			echo "LTR Harvest output scn file found"
			inharvest="${path2}/edta/${fasta_name}.mod.EDTA.raw/LTR/${fasta_name}.mod.rawLTR.scn"
		else
			echo "LTR Harvest output scn file found"
			echo "Please provide a LTR Harvest scn file and restart"
		fi
	else
		echo "EDTA output directory not found"
		echo "Please provide a LTR Harvest scn file and restart"
	fi
fi

#Make & cd to output directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#Run LTR_retriever using EDTA results
echo "Running LTR_retriever"
${path1}/annotation/LTR_retriever/LTR_retriever \
	-threads ${threads}	\
	-genome ${fasta} \
	-inharvest ${inharvest}

echo "Done"
