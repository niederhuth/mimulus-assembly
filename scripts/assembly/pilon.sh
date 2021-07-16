#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=500GB
#SBATCH --job-name pilon
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=40
rounds=4
datatype="wgs"
input="" #Can set to empty and script will find fasta in directory submitted

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/polishing/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/polishing/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${genotype}\\/// | sed s/\\/.*//)
assembly=$(pwd | sed s/^.*\\///)

#Extract reads from assembly job report
reads="$(grep reads: ../job_reports/${assembly}-*.SLURMout | head -1 | cut -d ' ' -f2)"

#Look for fasta file, there can only be one!
if [ -z ${input} ]
then
	echo "No input fasta provided, looking for fasta"
	if ls *.fa >/dev/null 2>&1
	then
		input=$(ls *fa | sed s/.*\ //)
		echo "Fasta file ${input} found"
	elif ls *.fasta >/dev/null 2>&1
	then
		input=$(ls *fasta | sed s/.*\ //)
		echo "Fasta file ${input} found"
	elif ls *.fna >/dev/null 2>&1
	then
		input=$(ls *fna | sed s/.*\ //)
		echo "Fasta file ${input} found"
	else
		echo "No fasta file found, please check and restart"
	fi
else
	echo "Input fasta: ${input}"
fi

#Loop over designated number of rounds for polishing
ref=../${input}
a=0
until [ ${a} -eq ${rounds} ]
do
	a=$(expr ${a} + 1)
	path2="pilon_${a}"
	echo "Round ${a} of polishing" 
	if [ -d ${path2} ]							
	then
		echo "Directory ${path2} exists."
		echo "Checking files."
		cd ${path2}
	else
		mkdir ${path2}
		cd ${path2}
	fi
	#Align data with minimap2
	if [ -f round_${a}.bam ]
	then
		echo "Round ${a} alignment found"
		echo "To rerun this step, delete ${path2}/round_${a}.paf and resubmit"
	else
		echo "Aligning with minimap2"
		bwa \
			-t ${threads} \
			-x ${preset} \
			${ref} \
			../../${reads} > round_${a}.paf
	fi
	#Polish with Pilon
	if [ -f pilon_${a}.fasta ]
	then
		echo "Round ${a} polishing found"
		echo "To rerun this step, delete ${path2}/pilon_${a}.fasta and resubmit"
	else
		echo "Polishing data with Pilon"
		pilon \
			--threads ${threads} \
			--genome ${ref} \
			--frags round_${a}.bam --unpaired unpaired.bam \
			--diploid \
			--fix all \
			--output pilon_${a}
	fi
	ref="../${path2}/pilon_${a}.fasta"
	cd ../
	echo "Round ${a} of polishing complete"
done

echo "Done"

