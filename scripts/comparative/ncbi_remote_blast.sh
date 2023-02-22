#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB
#SBATCH --job-name ncbi_remote_blast
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
fasta="$(pwd)/seqs/prot.fa"
email="niederhu@msu.edu"
blast_type="blastp"
evalue=0.0001
result_number=5000

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/phylo/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/phylo/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
path2=$(pwd | sed s/data.*/scripts/)
path3=blast

#cd to blast directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#Loop over the fasta file and run blast for each sequence
grep \> ${fasta} | sed s/\>// > blast_list
cat blast_list | while read seq
do
	if [ -d ${seq} ]
	then
		echo "Previous blast results for ${seq} found. Skipping."
		echo "To rerun this blast, delete the directory and resubmit."
	else
		echo "Running ${blast_type} on ${seq}"
		mkdir ${seq}
		cd ${seq}
		samtools faidx ${fasta} ${seq} > ${seq}.fa
		python ${path2}/comparative/py/NCBI_BLAST.py \
			--bp ${blast_type} \
			--evalue ${evalue} \
			--rn ${result_number} \
			--email ${email} \
			--output-prefix ${seq}-blastp \
			--input ${seq}.fa
		cd ..
	fi
done

echo "Done"
