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
path2="purge_dups"

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
if [ -s reads.paf.gz ]
then
	echo "Aligned reads found, proceeding to coverage statistics."
	echo "To repeat this step, delete ${path2}/reads.paf.gz and resubmit."
else
	echo "Aligning reads to assembly"
	minimap2 \
		-t ${threads} \
		-x ${preset} \
		../${asm} \
		${path1}/${reads} | gzip -c - > reads.paf.gz
fi

#Generate coverage stats
if [ -s PB.base.cov ]
then
	echo "Coverage stats found, proceeding to cutoff generation."
	echo "To repeat this step, delete ${path2}/PB.base.cov & ${path2}/PB.stat and resubmit."
else
	echo "Generating coverage statistics"
	pbcstat *.paf.gz #produces *.base.cov and *.stat files
fi

#Set cutoffs based on coverage
if [ -s cutoffs ]
then
	echo "Cutoffs file found, proceeding to fasta splitting."
	echo "To repeat this step, delete ${path2}/cutoffs and resubmit."
else
	echo "Generating cutoffs file"
	calcuts PB.stat > cutoffs 2>calcults.log
fi

#Split fasta sequence on 'Ns'
if [ -s fasta.split ]
then
	echo "Split fasta found, proceeding to self-alignment."
	echo "To repeat this step, delete ${path2}/fasta.split and resubmit."
else
	echo "Splitting fasta"
	split_fa ../${asm} > fasta.split
fi

#Align genome to itself
if [ -s fasta.split.self.paf.gz ]
then
	echo "Self-alignment found, proceding to duplicate purging."
	echo "To repeat this step, delete ${path2}/fasta.split.self.paf.gz and resubmit."
else
	echo "Aligning genome to itself"
	minimap2 \
		-t ${threads} \
		-x asm5 \
		-DP fasta.split \
		fasta.split | gzip -c - > fasta.split.self.paf.gz
fi

#Purge duplicates
if [ -s dups.bed ]
then
	echo "Duplicates bed found, proceeding to retrieve sequences."
	echo "To repeat this step, delete ${path2}/dups.bed and resubmit."
else
	echo "Running purge_dups"
	purge_dups \
		-2 \
		-T cutoffs \
		-c PB.base.cov \
		fasta.split.self.paf.gz > dups.bed 2> purge_dups.log
fi

#Retrieve sequences
if [ -s ${asm}.purge.fa ]
then
	echo "Purged fasta found, skipping."
	echo "To repeat this step, delete ${path2}/${asm}.purge.fa and resubmit."
else
	echo "Getting purged duplicate sequences"
	get_seqs \
		-l ${minimum_length} \
		-e dups.bed ../${asm} 
	mkdir hap
	mv hap.fa hap/
fi

echo "Done"


