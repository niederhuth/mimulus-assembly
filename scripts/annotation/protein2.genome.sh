#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500GB
#SBATCH --job-name protein2genome
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
fasta="repeatmasker/*.fa.masked"
proteins=$(ls -l proteins/*fa | sed s/.*\ //)
chunks=20
minintron=10
maxintron=3000
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
version=$(ls ${sample}-v*.fa | sed s/.*\-v// | sed s/.fa//) 
path2="exonerate"

#Make & cd to directory
if [ -d ${path2} ]
then
	cd ${path2}
else
	mkdir ${path2}
	cd ${path2}
fi

#Loop over protein sources and run exonerate
for i in ${proteins}
do
	outdir=$(echo ${i} | sed s/.*\\/// | sed s/\.fa//)
	a=1
	if [ -d ${outdir} ]
	then
		echo "${outdir} directory present, checking files"
	else
		mkdir ${outdir}
	fi
	
	while [ ${a} -le ${chunks} ]
	do
		if [ -z ${outdir}/chunk_${a} ]
		then
			echo "${outdir} chunk_${a} already complete, skipping to next chunk"
			a=$(expr ${a} + 1)
		else
			echo "Running exonerate protein2genome chunk ${a} on ${fasta}"
			exonerate \
			--model protein2genome \
			--bestn ${bestn} \
			--minintron ${minintron} \
			--maxintron ${maxintron} \
			--targetchunkid ${a} \
			--targetchunktotal ${chunks} \
			--query ../${i} \
			--target ../${fasta} \
			--showtargetgff yes \
			--showalignment no \
			--showvulgar no \
			--ryo ${ryo} > ${outdir}/chunk_${a}
		fi
		a=$(expr ${a} + 1)
	done
done

echo "Done"
