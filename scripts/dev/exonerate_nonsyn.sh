#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500GB
#SBATCH --job-name S1_L1_nonsynt_hits
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
model=protein2genome
seqs="L1_nonsynt-protein.fa"
genome="S1-v1.fa"
minintron=10
maxintron=5000
bestn=5
ryo=">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/maker/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/maker/lib:${LD_LIBRARY_PATH}"

#Set temporary directories for large memory operations
export TMPDIR=$(pwd | sed s/data.*/data/)
export TMP=$(pwd | sed s/data.*/data/)
export TEMP=$(pwd | sed s/data.*/data/)

#The following shouldn't need to be changed, but should set automatically
path2="nonsyntenic"

#Make & cd to directory
if [ -d ${path2} ]
then
	cd ${path2}
else
	mkdir ${path2}
	cd ${path2}
fi

#Run Exonerate
echo "Running exonerate"
exonerate \
	--model ${model} \
	--bestn 5 \
	--minintron 10 \
	--maxintron 5000 \
	--query ${seq} \
	--target ${genome} \
	--showtargetgff yes \
	--showalignment no \
	--showvulgar no \
	--ryo "${ryo}" > target_${genome/-*/}_query_${seqs/.fa}.output

echo "Done"

