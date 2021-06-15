#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=500GB
#SBATCH --job-name nextDenovo
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
#NA

#This should match the dataype in the misc/samples.csv file
#Options include:
#"ont" = Raw Nanopore
#"ont-cor" = Corrected Nanopore
#"pac" = raw PacBio
#"pac-cor" = Corrected PacBio
#"hifi" = PacBio HiFi
datatype="ont" 

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/assembly/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/assembly/lib:${LD_LIBRARY_PATH}"

#Set temporary directories for large memory operations
export TMPDIR=$(pwd | sed s/cre-genomes.*//)
export TMP=$(pwd | sed s/cre-genomes.*//)
export TEMP=$(pwd | sed s/cre-genomes.*//)

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/^.*\\///)
path2="nextDenovo"
reads="fastq/${datatype}/filtered.fastq.gz"

#Declare reads
echo "reads: ${reads}"

#Run nextDenovo
if ls ${path2} >/dev/null 2>&1 
then
	cd ${path2}
else
	mkdir ${path2}
	cd ${path2}
	cp ${path1}/nextDenovo/${species}_${genotype}.cfg ./
	cp ${path1}/nextDenovo/${species}_${genotype}.fofn ./
fi

nextDenovo ${species}_${genotype}.cfg

echo "Done"
