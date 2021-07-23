#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200GB
#SBATCH --job-name lrscaf
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20
datatype="ont"
reads="fastq/ont/clean.fastq.gz"
output="lrscaf"
input="" #Input assembly, if left blank, will search for fasta file

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/scaffolding/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/scaffolding/lib:$LD_LIBRARY_PATH"
#Path to LRScaf
LRScaf="${conda}/envs/scaffolding/bin/LRScaf-1.1.11.jar"

#The following shouldn't need to be changed, but should set automatically
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${genotype}\\/// | sed s/\\/.*//)
path1=$(pwd | sed s/${genotype}.*/${genotype}/)

#Change preset based on datatype
if [ ${datatype} = "ont" ]
then
	preset="map-ont"
elif [ ${datatype} = "ont-cor" ]
then
	preset="map-ont"
elif [ ${datatype} = "pac" ]
then
	preset="map-pb"
elif [ ${datatype} = "pac-cor" ]
then
	preset="map-pb"
elif [ ${datatype} = "hifi" ]
then
	preset="map-pb"
else
	echo "Do not recognize ${datatype}"
	echo "Please check and resubmit"
fi

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

#Align reads to assembly
if [ -s reads.paf.gz ]
then
	echo "Previous alignment found, proceeding to LRScaf"
	echo "To rerun this step, delete reads.paf.gz and resubmit"
else
	echo "Mapping reads to assembly"
	minimap2 \
		-t ${threads} \
		-x ${preset} \
		${input} \
		${path1}/${reads} > ${output}/reads.paf
fi

if [ -d lrscaf ]
then
	echo "Output directory ${output} found"
	echo "To rerun this step, delete and resubmit"
else
	echo "Scaffolding with LRScaf"
	java -jar ${LRScaf} \
		-c ${input} \
		-a ${output}/reads.paf \
		-o ${output} \
		-t mm \
		-p ${threads}
fi

echo "Done"
