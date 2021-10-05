#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100GB
#SBATCH --job-name ragtag_scaffold
#SBATCH --output=%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=10
aligner="nucmer"
min_len=1000
merge_dist=100000
min_gap=100
max_gap=100000
infer_gaps=TRUE #TRUE or FALSE
min_orient_score=0.0
min_location_score=0.0
min_grouping_score=0.2

#Set to whatever you used for genome assembly
datatype="ont"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/scaffolding/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/scaffolding/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
assembly=$(pwd | sed s/^.*\\///)
path2=$(pwd | sed s/${genotype}\\/${sample}.*/${genotype}\\/${sample}/)

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

#Change settings based on if infer_gaps is set to yes
if [ ${infer_gaps} = "TRUE" ]
then
	gaps="-r -g ${min_gap} -m ${max_gap} "
else
	gaps=
fi

#Get list of genomes
genomes=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${condition} \
	-v e=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $4 == d && $5 == e) print $7}' \
	${path1}/samples.csv)

#Run ragtag correct
for i in ${genomes}
do
	echo "Running ragtag scaffold with ${i} as reference genome"
	path3=$(pwd | sed s/${species}\\/.*/${species}\\/${i}/)
	version=$(ls ${path3}/ref/${i}-v*.fa | sed s/.*\-v// | sed s/.fa//)
	ref="${path3}/ref/${i}-v${version}.fa"
	#Run ragtag scaffold
	echo "Running ragtag scaffold"
	ragtag.py scaffold \
		-o ${i}_ragtag_scaffold \
		--aligner ${aligner} \
		-f ${min_len} \
		-d ${merge_dist} \
		-i ${min_grouping_score} \
		-a ${min_location_score} \
		-s ${min_orient_score} ${gaps}\
		-u \
		${ref} ${input}
done

echo "Done"


