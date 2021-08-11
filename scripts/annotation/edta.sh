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
alt_cds=""
TElibrary="" #Curated library

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/EDTA/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/EDTA/lib:${LD_LIBRARY_PATH}"

#Set temporary directories for large memory operations
export TMPDIR=$PBS_O_WORKDIR
export TMP=$PBS_O_WORKDIR
export TEMP=$PBS_O_WORKDIR

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${genotype}\\/// | sed s/\\/.*//)
version=$(ls ${genotype}-v*.fa | sed s/.*\-v// | sed s/.fa//) 
path2="edta"

#Set species for TIR-learner
if [ ${species} == "Zmays" ] 
then 
	TIRspecies="Maize"
elif [ ${species} == "Osativa" ]
then 
	TIRspecies="Rice"
else 
	TIRspecies="others"
fi

#Check 
if [ -z ${TElibrary} ]
then
	if [ ${species} == "Zmays" ]
	then
		TElibrary="${conda}/envs/EDTA/share/EDTA/database/maizeTE11122019"
		echo "Using currated library ${TElibrary}"
	elif [ ${species} == "Osativa" ]
	then
		TElibrary="${conda}/envs/EDTA/share/EDTA/database/rice6.9.5.liban"
		echo "Using currated library ${TElibrary}"
	else
		echo "No currated library provided, proceeding with denovo detection"
		args=""
	fi
else
	args="--curatedlib ${TElibrary}"
	echo "Using currated library ${TElibrary}"
fi

#Get CDS sequences
#This will adjust depending on whether or not these are available
if [ -f annotations/${genotype}-v${version}-cds.fa ]
then
	cds_seqs="annotations/${genotype}-v${version}-cds.fa"
	echo "CDS sequences found"
	echo "Using ${cds_seqs}"
	args="${args} --cds ${cds_seqs}"
else
	cds_seqs=${alt_cds}
	echo "No CDS sequences found"
	echo "Is this an un-annotated genome?"
	if [ -s ${alt_cds} ]
	then
		echo "Using ${alt_cds} as cds."
		echo "If this is a mistake, check if CDS sequences are in proper place and resubmit."
		args="${args} --cds ${alt_cds}"
	else
		echo "No CDS sequences provided"
	fi
fi

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
if [ -d ${path2} ]
then
	echo "Previous EDTA results found."
	echo "Restarting from previous run."
	echo "To repeat this analysis, delete ${path2} and resubmit."
	echo "Rusuming EDTA for ${fasta}"
	cd ${path2}
	EDTA.pl ${args} \
		--genome ${fasta} \
		--species ${TIRspecies} \
		--step all \
		--overwrite 0 \
		--sensitive 1 \
		--anno 1 \
		--evaluate 0 \
		--force 1 \
		--threads ${threads}
else
	mkdir ${path2}
	cp ${cds_seqs} ${path2}/cds.fa
	cp ${fasta} ${path2}/${fasta}
	cd ${path2}
	echo "Running EDTA for ${fasta}"
	EDTA.pl ${currated} \
		--genome ${fasta} \
		--species ${TIRspecies} \
		--step all \
		--overwrite 1 \
		--sensitive 1 \
		--anno 1 \
		--evaluate 0 \
		--force 1 \
		--threads ${threads}
fi

echo "Done"
