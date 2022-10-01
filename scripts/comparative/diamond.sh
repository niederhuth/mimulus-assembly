#!/bin/bash --login
#SBATCH --time=128:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=50GB
#SBATCH --job-name diamond_blastp
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
datatype="proteins-primary" #Query seqs: proteins (all), proteins-primary (primary transcript)
genomes= #List targets as SPECIES_GENOTYPE, e.g. Mguttatus_S1; if left blank, will search in samples.csv
threads=50 #Number of threads to use
evalue=0.00001 #e-value cutoff
max_target_seqs=5 #max number of hits to retain
options="--unal 0 --more-sensitive" #additional diamond options
orthogroups="$(pwd | sed s/data.*/data/)/comparative/orthofinder/Results_Sep27/Phylogenetic_Hierarchical_Orthogroups/N0.tsv" #N0.tsv file, if blank will not run this step

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/orthofinder/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/orthofinder/lib:${LD_LIBRARY_PATH}"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
path2=$(pwd | sed s/data.*/data/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample="comparative"
condition="blastp"
path3="diamond_blastp"

#Make output directory and change to it
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#Get list of genomes
if [ -z ${genomes} ]
then
	genomes=$(awk -v FS="," \
		-v a=${species} \
		-v b=${genotype} \
		-v c=${sample} \
		-v d=${condition} \
		-v e=${datatype} \
		'{if ($1 == a && $2 == b && $3 == c && $4 == d && $5 == e) print $7}' \
		${path1}/samples.csv)
fi

#Get query sequences
query="${species}_${genotype}"
version=$(ls ${path2}/${species}/${genotype}/ref/annotations/${genotype}-v*-${datatype}.fa | \
	sed s/.*\-v// | sed s/\-${datatype}.fa//)
query_seqs="${path2}/${species}/${genotype}/ref/annotations/${genotype}-v${version}-${datatype}.fa"

#Run diamond
for i in ${genomes}
do
	mkdir ${i}
	cd ${i}
	#Find target sequences
	path4="${path2}/${i/_*/}/${i/*_/}"
	version2=$(ls ${path4}/ref/annotations/${i/*_/}-v*-${datatype}.fa | \
		sed s/.*\-v// | sed s/\-${datatype}.fa//)
	target_seqs="${path4}/ref/annotations/${i/*_/}-v${version2}-${datatype}.fa"
	dmnd_db="${path4}/ref/annotations/${i/*_/}-v${version2}-${datatype}.dmnd"
	echo "Target is ${i}-v${version2}"

	#Check for diamond database
	echo "Checking for ${i} diamond database"
	if [ -f ${path4}/ref/annotations/${i/*_/}-v${version2}-${datatype}.dmnd ]
	then
		echo "Found diamond database"
	else
		echo "No diamond database found"
		echo "Makind diamond database"
		diamond makedb \
			--in ${target_seqs} \
			--db ${dmnd_db}
	fi
	
	#Run diamond blast
	if [ -f ${query}-${i}.m8 ]
	then
		echo "${query}-${i}.m8 already exists, skipping"
		echo "To rerun diamond blastp on ${query}-${i}.m8, then delete ${query}-${i}.m8 and resubmit"
	else
		options="--threads ${threads} ${options}"
		echo "Running diamond blastp against ${i}-v${version2}"
		diamond blastp \
			${options} \
			--db ${dmnd_db} \
			--query ${query_seqs} \
			--out ${query}-${i}.m8 \
			--un ${query}-${i}-un.fa \
			--evalue ${evalue} \
			--max-target-seqs ${max_target_seqs}
	fi
	if [[ ! -z ${orthogroups} ]]
	then
		echo "Filtering results for othogroups"
		cat ${query}-${i}.m8 | while read line
		do
			a=$(echo ${line} | cut -d ' ' -f1)
			b=$(echo ${line} | cut -d ' ' -f2)
			if [[ $(grep ${a} ${orthogroups} | cut -f1) = $(grep ${b} ${orthogroups} | cut -f1) ]]
			then
				echo ${line} | tr ' ' '\t' >> ${query}-${i}_orthogroup_filtered.m8
			fi
		done
	fi
	cd ../
done

