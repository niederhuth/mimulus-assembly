#!/bin/bash --login
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=150GB
#SBATCH --job-name convert_gam
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=10
index="$(pwd | sed s/Mguttatus.*/Mguttatus/)/pangenome/giraffe/index.giraffe.gbz"
output_format=BAM #BAM, SAM, CRAM, GAMP, GAF
datatype="poolseq"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/pangenome/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/pangenome/lib:$LD_LIBRARY_PATH"

#Other variables, these should not have to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*data\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*data\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)

#Set output format
arguments="-t ${threads} -x ${index}"
#Set output format
if [ ${output_format} = "BAM" ]
then
	arguments="${arguments} -b"
elif [ ${output_format} = "SAM" ]
then
	arguments="${arguments} -s"
elif [ ${output_format} = "CRAM" ]
then
	arguments="${arguments} c"
elif [ ${output_format} = "GAMP" ]
then
	arguments="${arguments} -m"
elif [ ${output_format} = "GAF" ]
then
	arguments="${arguments} -G"
else
	echo "Output format not specified!"
	exit 1
fi

#Set input
input="giraffe/${sample}_${datatype}.gam"
#Set output
output="giraffe/${sample}_${datatype}.${output_format}"

# Compute the read support from the gam
echo "Converting to ${output_format}"
vg surject ${arguments} ${input} > ${output}

echo "Done"
