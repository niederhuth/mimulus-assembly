#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200GB
#SBATCH --job-name cactus-minigraph
#SBATCH --job_reports/output=%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
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
if [[ ! -d containers ]]
then
	mkdir containers
fi
export UDOCKER_CONTAINERS=$(pwd)/containers

#The following shouldn't need to be changed, but should set automatically
#Set new path for within container
path1=$(echo ${PBS_O_WORKDIR} | sed s/.*data/data/)
#cd to that path
cd ${path1}
#Set rest of variables
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
path2=cactus_pg_gpu

#Check for docker container, if it doesnt exist, create it
if [[ ! -d containers/minigraph_cactus ]]
then
	udocker create --name=cactus_pg_gpu quay.io/comparative-genomics-toolkit/cactus:v2.4.3-gpu
fi
#Run docker container
udocker run --volume=$(pwd | sed s/data.*//) cactus_pg_gpu
#Change directories within container
cd /data/${path1}

#Create output directory
if [[ ! -d ${path2} ]]
then
	mkdir ${path2}
fi

#Set output files
outputGFA=${path2}/${species}-pg.gfa.gz
logFile={path2}/${species}-pg-minigraph.log

#Get the reference genome
if [ -z ${reference} ]
then
	reference=$(grep -A 1 Haploid ${seqFile} | tail -n 1 | cut -f1)
fi

#Run cactus-minigraph
echo "Running cactus-minigraph"
cactus-minigraph ${path2} ${seqFile} ${outputGFA} \
	--reference ${reference} \
	--logFile ${logFile} \
	--mapCores ${threads}

echo "Done"
