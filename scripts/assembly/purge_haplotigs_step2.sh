#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200GB
#SBATCH --job-name purge_haplotigs_step1
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20
datatype="ont"
depth=200
minimum_length="10k"
fasta=""

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/purge_haplotigs/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/purge_haplotigs/lib:$LD_LIBRARY_PATH"

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
if [ -z ${fasta} ]
then
	echo "No input fasta provided, looking for fasta"
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
else
	echo "Input fasta: ${fasta}"
fi

#Get coverage stats
if [ -s coverage_stats.csv ]
then
	echo "Coverage stats found, proceeding to haplotig purging."
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
		-genome ${fasta}
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
