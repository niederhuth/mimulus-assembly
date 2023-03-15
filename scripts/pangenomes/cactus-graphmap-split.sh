#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100GB
#SBATCH --job-name cactus-graphmap-split
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=10
name=Mguttatus
seqFile=/data/misc/Mguttatus-pangenome-seqs.txt
reference=

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/pangenome/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/pangenome/lib:$LD_LIBRARY_PATH"
#Export path to UDOCKER_DIR. All images will be downloaded and installed here
export UDOCKER_DIR=${conda}/envs/pangenome/udocker
#Export path to UDOCKER_CONTAINERS. All containers will be saved there
mkdir containers
export UDOCKER_CONTAINERS=$(pwd)/containers

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
#Set output directory
path2=cactus_pg
#Create output directory
if [[ ! -d ${path2} ]]
then
	mkdir ${path2}
fi

#Copy over the seqFile
if [[ ! -f ${path2}/seqFile.txt ]]
	cp ${path1}/${seqFile} ${path2}/seqFile.txt
fi

#Get the reference genome
if [[ ! -f ${path2}/reference.txt ]]
then
	if [ -z ${reference} ]
	then
		echo $(grep -A 1 Haploid ${path2}/seqFile.txt | tail -n 1 | cut -f1) > ${path2}/reference.txt
	else
		echo ${reference} > ${path2}/reference.txt
	fi
fi

#Set the name
if [[ ! -f ${path2}/name.txt ]]
then
	echo ${name} > ${path2}/name.txt
fi

#Check for docker container, if it doesnt exist, create it
if [[ ! -d containers/cactus_pg ]]
then
	udocker create --name=cactus_pg quay.io/comparative-genomics-toolkit/cactus:v2.4.3
fi
#Run docker container
udocker run --volume=$(PBS_O_WORKDIR):/data cactus_pg

#Set files & variables
path2=cactus_pg
seqFile=${path2}/seqFile.txt
reference=$(head -1 ${path2}/reference.txt)
name=$(head -1 ${path2}/reference.txt)
inputGFA=${path2}/${name}-pg.gfa.gz
inputPAF=${path2}/${name}-pg.paf
outDir=${path2}/split
jobstore=${path2}/jobstore
logFile={path2}/${name}-pg-graphmap-split.log

#Run cactus-minigraph
echo "Running cactus-graphmap-split"
cactus-graphmap-split ${jobstore} ${seqFile} ${inputGFA} ${inputPAF} \
	--reference ${reference} \
	--outDir ${outDir}

echo "Done"
