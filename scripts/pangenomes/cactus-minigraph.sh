#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100GB
#SBATCH --job-name cactus-minigraph
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=10
name=Mguttatus
seqFile=Mguttatus-pangenome-seqs.txt
reference=

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/pangenome/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/pangenome/lib:$LD_LIBRARY_PATH"
#Export path to UDOCKER_DIR. All images will be downloaded and installed here
export UDOCKER_DIR=${conda}/envs/pangenome/udocker
#Export path to UDOCKER_CONTAINERS. All containers will be saved there
if [[ ! -d containers ]]
then
	mkdir containers
fi
export UDOCKER_CONTAINERS=$(pwd)/containers

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*//)
path2=$(pwd | sed s/data.*/misc/)
path3=$(pwd)/cactus_pg
seqFile=${path2}/${seqFile}

#Check for docker container, if it doesnt exist, create it
if [[ ! -d containers/cactus_pg ]]
then
	udocker create --name=cactus_pg quay.io/comparative-genomics-toolkit/cactus:v2.4.3
fi

#Create output directory
if [[ ! -d ${path3} ]]
then
        mkdir ${path3}
fi

#Check for seqFile and copy over is not already
#if [[ ! -f ${path3}/${seqFile} ]]
#then
 #       cp ${path1}/${seqFile} ${path3}/${seqFile}
 #       seqFile=${path3}/${seqFile}
#else
 #       seqFile=${path3}/${seqFile}
#fi

#Get the reference genome
if [ -z ${reference} ]
then
        reference=$(grep -A 1 Haploid ${seqFile} | tail -n 1 | cut -f1)
fi

#Set files & variables
outputGFA=${path3}/${name}-pg.gfa.gz
jobstore=${path3}/jobstore
logFile=${path3}/${name}-pg-minigraph.log

#Run cactus-minigraph
echo "Running cactus-minigraph"
udocker run \
	--volume=${path1}:/data \
	cactus_pg \
	cactus-minigraph ${jobstore} ${seqFile} ${outputGFA} \
		--reference ${reference} \
		--logFile ${logFile} \
		--mapCores ${threads}

echo "Done"
