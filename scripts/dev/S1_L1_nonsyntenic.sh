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
query_chunks=10
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

#Run Exonerate
a=1
while [ ${a} -le ${query_chunks} ]
do
	if [ -s target_${genome/-*/}_query_${seqs/.fa}_chunk_${a}.output ]
	then
		if [ $(wc -l target_${genome/-*/}_query_${seqs/.fa}_chunk_${a}.output | cut -d ' ' -f1) -gt 20 ]
		then
			echo "target_${genome/-*/}_query_${seqs/.fa}_chunk_${a}.output already complete" 
			echo "Skipping to next chunk"
			a=$(expr ${a} + 1)
		else
			rm target_${genome/-*/}_query_${seqs/.fa}_chunk_${a}.output
		fi
	fi
	if [ ! -s target_${genome/-*/}_query_${seqs/.fa}_chunk_${a}.output ]
	then
		if [ ${a} -le ${query_chunks} ]
		then
			echo "Aligning chunk_${a}"
			exonerate \
				--model ${model} \
				--bestn 5 \
				--minintron 10 \
				--maxintron 5000 \
				--query ${seq} \
				--target ${genome} \
				--querychunkid ${a} \
				--querychunktotal ${query_chunks} \
				--showtargetgff yes \
				--showalignment no \
				--showvulgar no \
				--ryo "${ryo}" > target_${genome/-*/}_query_${seqs/.fa}_chunk_${a}.output
			a=$(expr ${a} + 1)
		fi
	fi
done

echo "Done"

