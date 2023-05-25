#!/bin/bash --login
#SBATCH --time=128:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100GB
#SBATCH --job-name genespace
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
genomes= #Space separated list of Species_Genotype, if left blank look for in misc/samples.csv
ploidy= #Space separated list of ploidy level (e.g. 1 1 2 2), if left blank look for in misc/samples.csv

#Change to current directory
#cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/orthofinder/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/orthofinder/lib:${LD_LIBRARY_PATH}"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed 's/data.*/misc/')
path2=$(pwd | sed 's/data.*/data/')
path3=$(pwd | sed 's/data.*/scripts\/comparative\/R/')
species=$(pwd | sed 's/^.*\/data\///' | sed 's/\/.*//')
genotype="genespace"
sample="genespace"
condition="genespace"
datatype="proteins-primary"
path4="genespace"

#Make output directory and change to it
if [ -d ${path4} ]
then
	cd ${path4}
else
	mkdir ${path4}
	cd ${path4}
fi

#Get list of genomes
if [ -z ${genomes} ]
then
	read -a genomes <<< $(awk -v FS="," \
		-v a=${species} \
		-v b=${genotype} \
		-v c=${sample} \
		-v d=${condition} \
		-v e=${datatype} \
		'{if ($1 == a && $2 == b && $3 == c && $4 == d && $5 == e) print $7}' \
		${path1}/samples.csv)
fi

#Get list of genomes
if [ -z ${ploidy} ]
then
	read -a ploidy <<< $(awk -v FS="," \
		-v a=${species} \
		-v b=${genotype} \
		-v c=${sample} \
		-v d=${condition} \
		-v e=${datatype} \
		'{if ($1 == a && $2 == b && $3 == c && $4 == d && $5 == e) print $10}' \
		${path1}/samples.csv)
fi

#Prepare input files
mkdir rawGenomes
echo "Species,Genotype,Version,Ploidy" > genomes.csv
for i in ${!genomes[@]}
do
	echo "Copying files for ${genomes[i]}"
	path4="${path2}/${genomes[i]/_*/}/${genomes[i]/*_/}/ref/annotations/"
	path5="rawGenomes/${genomes[i]/_*/}/${genomes[i]/*_/}/annotation"
	version=$(ls ${path4}/${genomes[i]/*_/}-v*${datatype}.fa | sed 's/.*\-v//' | sed s/\-${datatype}.fa//)
	mkdir rawGenomes/${genomes[i]/_*/} rawGenomes/${genomes[i]/_*/}/${genomes[i]/*_/} ${path5}
	cp ${path4}/${genomes[i]/*_/}-v${version}.gff ${path5}/
	cp ${path4}/${genomes[i]/*_/}-v${version}-${datatype}.fa ${path5}/
	echo "${genomes[i]/_*/},${genomes[i]/*_/},${genomes[i]/*_/},${ploidy[i]}" >> genomes.csv
done

#Run genespace
echo "Running genespace"
Rscript ${path3}/genespace.R

#Classify genespace orthologs
echo "Classifying genespace orthologs"
Rscript ${path3}/classify_genespace.R


echo "Done"

