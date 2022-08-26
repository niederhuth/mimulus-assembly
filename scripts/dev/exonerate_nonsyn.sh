#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500GB
#SBATCH --job-name protein2genome-Athaliana
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
fasta=$(ls -l repeatmasker/*.fa.masked | sed s/.*\ //)
proteins=$(ls -l proteins/*.fa | sed s/.*\ //)

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
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*${species}\\/${genotype}\\/// | sed s/\\/.*//)
path2="nonsyntenic"

#Make & cd to directory
if [ -d ${path2} ]
then
	cd ${path2}
else
	mkdir ${path2}
	cd ${path2}
fi

#Extract sequences for GOI
seqtk subseq ${} ${genelist} > ${}

#Search genome
echo ""
exonerate \
	--model est2genome \
	--bestn 5 \
	--minintron 10 \
	--maxintron 5000 \
	--query L1_nonsyn-cds.fa \
	--target ../S1/pseudomolecule/S1-v1.fa \
	--showtargetgff yes \
	--showalignment no \
	--showvulgar no \
	--ryo "${ryo}" > S1_L1_nonsyn-cds.output

echo "Done"

