#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=50GB
#SBATCH --job-name est_genotype_error
#SBATCH --output=%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables


#In general dont change this, unless using a similar datatype
#This should match the dataype in the misc/samples.csv file
datatype="radseq"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/scaffolding/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/scaffolding/lib:${LD_LIBRARY_PATH}"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
assembly=$(pwd | sed s/^.*\\///)
path2=$(pwd | sed s/${genotype}\\/${sample}\\/.*/${genotype}\\/${sample}/)
path3=${path1}/genetic_map/${datatype}

#Create bad_marks.txt
touch bad_marks.txt

#
for i in 
do
	python ${path}/dev/GOOGA/hmm_scaff_likelihood.py ${i/.txt//}
done

