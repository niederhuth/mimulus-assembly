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
seqFile=/data/misc/Mguttatus-pangenome-seqs.txt
reference=

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/pangenomes/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/pangenomes/lib:$LD_LIBRARY_PATH"
#Export path to UDOCKER_DIR. All images will be downloaded and installed here
export UDOCKER_DIR=${conda}/envs/pangenome/udocker
#Export path to UDOCKER_CONTAINERS. All containers will be saved there
mkdir containers
export UDOCKER_CONTAINERS=$(pwd)/containers

#The following shouldn't need to be changed, but should set automatically
#Set new path for within container
path1=$(echo ${PBS_O_WORKDIR} | sed s/.*data/data/)
#cd to that path
cd ${path1}
#Set rest of variables
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
path2=cactus_pg

#Create docker container
udocker create --name=cactus_pg quay.io/comparative-genomics-toolkit/cactus:v2.4.3
#Run docker container
udocker run --volume=$(pwd | sed s/data.*//) cactus_pg
#Change directories within container
cd /data/${path1}

#Create output directory
if [[ ! -d ${path2} ]]
then
	mkdir ${path2}
fi

#Set output files
inputGFA=${path2}/${species}-pg.gfa.gz
inputPAF=${path2}/${species}-pg.paf
outDir=${path2}/split
logFile={path2}/${species}-pg-graphmap-split.log

#Get the reference genome
if [ -z ${reference} ]
then
	reference=$(grep -A 1 Haploid ${seqFile} | tail -n 1 | cut -f1)
fi

#Run cactus-minigraph
echo "Running cactus-graphmap-split"
cactus-graphmap-split ${path2}/jobstore ${seqFile} ${inputGFA} ${inputPAF} \
	--reference ${reference} \
	--outDir ${outDir} \
	--mapCores ${threads} 

echo "Done"
