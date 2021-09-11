#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200GB
#SBATCH --job-name purge_dups
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20
datatype="ont"
reads="fastq/ont/clean.fastq.gz"
asm="flye2/medaka_racon_2/consensus.fasta" #Input assembly

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/assembly/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/assembly/lib:$LD_LIBRARY_PATH"
#Path to LRScaf
LRScaf="${conda}/envs/assembly/bin/LRScaf-1.1.11.jar"

#The following shouldn't need to be changed, but should set automatically
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
path1=$(pwd | sed s/${genotype}\\/${sample}.*/${genotype}\\/${sample}/)
path2="lrscaf"

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

#Make output directory
if [ -d ${path2} ]
then
	cd ${path2}
else
	mkdir ${path2}
	cd ${path2}
fi

#Align reads to assembly
if [ -f reads.paf.gz ]
then
	echo "Previous alignment found, proceeding to LRScaf"
	echo "To rerun this step, delete reads.paf.gz and resubmit"
else
	echo "Mapping reads to assembly"
	minimap2 \
		-t ${threads} \
		-x ${preset} \
		${path1}/${asm} \
		${path1}/${reads} > reads.paf
fi

if [ -f test.txt ]
then
	echo ""
	echo "to rerun this step, delete and resubmit"
else
	echo "Scaffolding with LRScaf"
	java -jar ${LRScaf} \
		-c ${path1}/${asm} \
		-a reads.paf \
		-o ./ \
		-t mm \
		-p ${threads}
fi

echo "Done"
