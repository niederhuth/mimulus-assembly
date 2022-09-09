#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200GB
#SBATCH --job-name nucmer
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
ref= #reference genome, if left blank, will look for in the ref directory for that genotype
threads=20
breaklen=500 #distance to attempt to extend poor scoring regions before giving up, default 200
mincluster=100 #min length of a cluster of matches, default 65
minmatch=50 #min length for single match, default 20

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/wga/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/wga/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="genome_alignment"
datatype="genome"
path2=$(pwd | sed s/data.*/data/)
path3="nucmer"

#Check and look for ref genome
if [ -z ${ref} ]
then
	echo "No reference genome provided, looking for reference genome"
	ref="${path2}/${species}/${genotype}/ref/${genotype}-v*.fa"
	echo "Reference genome ${ref} found"
fi

#Check for and make/cd working directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
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

#Run nucmer alignments
for i in ${genomes}
do
	mkdir ${i}_alignment
	cd ${i}_alignment
	#Set path to query genome
	query="${path2}/$(echo ${i} | sed s/_/\\//)/ref/${i/*_/}-v*.fa"
	#Run nucmer
	echo "Aligning ${i} against ${species}_${genotype} with nucmer"
	nucmer \
		--threads ${threads} \
		--maxmatch \
		-b ${breaklen} \
		-c ${mincluster} \
		-l ${minmatch} \
		-p ${i}.delta \	
		${ref} \
		${query}
	cd ../
done

echo "Done"
