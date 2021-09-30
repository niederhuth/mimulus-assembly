#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=500GB
#SBATCH --job-name sspace
#SBATCH --output=%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=40
extension=0 #extend the contigs of -s using paired reads in -l (-x 1=extension, -x 0=no extension, default -x 0)
insert_size=379 #mean insert size
insert_error=0.23 #insert size error, I take SD/mean insert size
read_orientation=FR

#In general dont change this, unless using a similar datatype
#This should match the dataype in the misc/samples.csv file
datatype="wgs"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/test/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/test/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
assembly=$(pwd | sed s/^.*\\///)
path2=$(pwd | sed s/${genotype}\\/${sample}\\/.*/${genotype}\\/${sample}/)
path3=sspace

#Fastq files, these should not have to be changed, but should set automatically
path4="${path2}/fastq/${datatype}"
t1="${path4}/trimmed.1.fastq.gz"
t2="${path4}/trimmed.2.fastq.gz"

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

#Make and cd to output directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#Creat library file
echo "lib1 ${t1} ${t2} ${insert_size} ${insert_error} ${read_orientation}" | tr ' ' '\t' > libraries.tsv

#Run SSPACE
echo ""
perl $HOME/Dev/sspace_basic/SSPACE_Basic.pl \
	-l libraries.tsv \
	-s ../${input} \
	-x ${extension} \
	-T ${threads}



