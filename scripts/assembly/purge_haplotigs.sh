#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200GB
#SBATCH --job-name purge_dups
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20
datatype="ont"
depth=200
minimum_length="10k"
asm=""

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/scaffolding/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/scaffolding/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${genotype}\\/// | sed s/\\/.*//)
assembly=$(pwd | sed s/^.*\\///)
path1=$(pwd | sed s/${genotype}.*/${genotype}/)
path2="purge_haplotigs"

#Output location
echo "Purging Duplicates for ${species} ${genotype} ${sample} ${assembly}"

#Extract reads from assembly job report
reads="$(grep reads: ${path1}/job_reports/${sample}-*.SLURMout | head -1 | cut -d ' ' -f2)"

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
if [ -z ${asm} ]
then
	echo "No input fasta provided, looking for fasta"
	if ls *.fa >/dev/null 2>&1
	then
		asm=$(ls *fa | sed s/.*\ //)
		echo "Fasta file ${asm} found"
	elif ls *.fasta >/dev/null 2>&1
	then
		asm=$(ls *fasta | sed s/.*\ //)
		echo "Fasta file ${asm} found"
	elif ls *.fna >/dev/null 2>&1
	then
		asm=$(ls *fna | sed s/.*\ //)
		echo "Fasta file ${asm} found"
	else
		echo "No fasta file found, please check and restart"
	fi
else
	echo "Input fasta: ${asm}"
fi

#Create output directory and change directory
if [ -d ${path2} ]
then
	cd ${path2}
else
	mkdir ${path2}
	cd ${path2}
fi

#Align reads to assembly
if [ -s aligned.bam ]
then
	echo "Aligned reads found, proceeding to coverage statistics."
	echo "To repeat this step, delete ${path2}/aligned.bam and resubmit."
else
	echo "Aligning reads to assembly"
	minimap2 \
		-t ${threads} \
		-x ${preset} \
		../${asm} \
		${path1}/${reads} | gzip -c - > aligned.bam
fi

#Generate read-depth histogram
if [ -s aligned.bam.genecov ]
then
	echo "Aligned reads found, proceeding to read-depth histogram."
	echo "To repeat this step, delete ${path2}/aligned.bam.genecov and resubmit."
else
	echo "Generating read-depth histogram"
	purge_haplotigs hist \
		-bam aligned.bame \
		-genome ${asm} \
		-threads ${threads} \
		-depth ${depth}
fi

#Get coverage stats
if [ -s coverage_stats.csv ]
then
	echo "Aligned reads found, proceeding to coverage statistics."
	echo "To repeat this step, delete ${path2}/coverage_stats.csv and resubmit."
else
	echo "Generating coverage statistics"
	purge_haplotigs cov \
		-in aligned.bam.genecov \
		-low ${low} \
		-high ${high} \
		-mid ${mid} \
		-out coverage_stats.csv \
		-junk \
		-suspect
fi

#Purge duplicates
if [ -s dups.bed ]
then
	echo "Aligned reads found, proceeding to coverage statistics."
	echo "To repeat this step, delete ${path2}/dups.bed and resubmit."
else
	echo "Purging haplotigs"
	purge_haplotigs purge \
		-genome ${asm}
		-coverage coverage_stats.csv \
		-threads ${threads} \
		-outprefix curated \
		-repeats ${repeats} \ BED-format file of repeats to ignore during analysis.
		-dotplots \ Generate dotplots for manual inspection.
		-bam \ Samtools-indexed bam file of aligned and sorted reads/subreads to the reference, required for generating dotplots.
		-align_cov \ Percent cutoff for identifying a contig as a haplotig. DEFAULT = 70
		-max_match \     Percent cutoff for identifying repetitive contigs. Ignored when using repeat annotations (-repeats). DEFAULT = 250
		-I                  Minimap2 indexing, drop minimisers every N bases, DEFAULT = 4G
		-v / -verbose       Print EVERYTHING.
		-limit_io           Limit for I/O intensive jobs. DEFAULT = -threads
		-wind_min           Min window size for BED coverage plots (for dotplots). DEFAULT = 5000
		-wind_nmax          Max windows per contig for BED coverage plots (for dotplots). DEFAULT = 200

fi

echo "Done"


