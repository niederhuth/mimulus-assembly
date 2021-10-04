#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=200GB
#SBATCH --job-name sspace
#SBATCH --output=%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
insert_size=379 #mean insert size
insert_error=0.23 #insert size error, I take SD/mean insert size
read_orientation=FR 
threads=40
extension=0 #extend the contigs of -s using paired reads (1=extension, 0=no extension, default -x 0)
dot_file=1 #1=make file, 0=dont make file
min_contig=500 #Minimum contig length used for scaffolding. Filters out contigs below this value (default -z 0)
min_lengths=5 #Minimum number of links (read pairs) to compute scaffold (default -k 5)
max_link_ratio=0.7  #Max link ratio between two contig pairs. Higher values less accurate (default -a 0.7)
min_overlap=15 #Minimum overlap required between contigs to merge adjacent contigs in a scaffold (default -n 15)
bowtie_gaps=2 #Maximum number of allowed gaps during mapping with Bowtie (default -g 0)

#In general dont change this, unless using a similar datatype
#This should match the dataype in the misc/samples.csv file
datatype="wgs"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/scaffolding/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/scaffolding/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/scripts/)
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

#Uncompress reads
echo "Uncompressing fastq files"
zcat ${t1} > r1.fastq
zcat ${t2} > r2.fastq
 
#Creat library file
echo "lib1 r1.fastq r2.fastq ${insert_size} ${insert_error} ${read_orientation}" | tr ' ' '\t' > libraries.tsv

#Run SSPACE
echo "Running SSPACE_Basic"
#perl ${path1}/assembly/sspace_basic/SSPACE_Basic.pl \
perl ${HOME}/Dev/sspace_basic/SSPACE_Basic.pl \
	-l libraries.tsv \
	-s ../${input} \
	-x ${extension} \
	-z ${min_contig} \
	-k ${min_lengths} \
	-a ${max_link_ratio} \
	-n ${min_overlap} \
	-g ${bowtie_gaps} \
	-T ${threads} \
	-p ${dot_file}

echo "Done"

