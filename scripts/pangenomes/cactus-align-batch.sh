#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200GB
#SBATCH --job-name cactus-align-batch
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
path1=/data/$(pwd | sed s/.*data/data/)
path2=/data/misc
path3=cactus_pg

#Check for docker container, if it doesnt exist, create it
if [[ ! -d containers/cactus_pg ]]
then
	udocker create --name=cactus_pg quay.io/comparative-genomics-toolkit/cactus:v2.4.3
fi

#Run docker container
udocker run \
	--env="path1=${path1}" \
	--env="path3=${path2}" \
	--env="path3=${path3}" \
	--env="name=${name}" \
	--env="reference=${reference}" \
	--env="seqFile=${seqFile}" \
	--env="threads=${threads}" \
	--volume=$(pwd | sed s/data.*//):/data \
	cactus_pg

#Change to working directory
cd ${path1}

#Create output directory
if [[ ! -d ${path3} ]]
then
	mkdir ${path3}
fi

#Check for seqFile and copy over is not already
if [[ ! -f ${path3}/${seqFile} ]]
then
	cp ${path2}/${seqFile} ${path3}/${seqFile}
	seqFile=${path3}/${seqFile}
else
	seqFile=${path3}/${seqFile}
fi

#Set files & variables
chromfile=
splitDir=${path3}/split
logFile=${path3}/${species}-pg-graphmap-split.log

#Get the reference genome
if [ -z ${reference} ]
then
	reference=$(grep -A 1 Haploid ${seqFile} | tail -n 1 | cut -f1)
fi

#Run cactus-minigraph
echo "Running cactus-graphmap"
cactus-align-batch ${path3}/jobstore ${chromfile} ${splitDir}


 --alignCores 16 --realTimeLogging --alignOptions "--pangenome --maxLen 10000 --reference Glycine_max_v4_0 --outVG" --logFile soybean-pg-${VERSION}-align.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 1000 --maxNodes 20 --betaInertia 0 --targetTime 1

echo "Done"
