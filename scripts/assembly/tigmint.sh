#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=500GB
#SBATCH --job-name tigmint
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20 #doesn't seem to want to use more than 6
datatype="ont"
input= #input fasta, if left blank, will look for it in current directory

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/scaffolding/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/scaffolding/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
assembly=$(pwd | sed s/^.*\\///)
path2=$(pwd | sed s/${genotype}\\/${sample}.*/${genotype}\\/${sample}/)
path3="tigmint"

#Get genome size estimate
genomeSize=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${condition} \
	-v e=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $4 == d && $5 == e) print $9}' \
	${path1}/samples.csv)
genomeSize2=$(python -c "print(${genomeSize/g/} * 1000000000)")

#Look for fasta file, there can only be one!
if [ -z ${input} ]
then
	echo "No input fasta provided, looking for fasta"
	if ls *.fa >/dev/null 2>&1
	then
		input=$(ls *fa | sed s/.*\ //)
		echo "Fasta file ${input} found"
	elif ls *.fasta >/dev/null 2>&1
	then
		input=$(ls *fasta | sed s/.*\ //)
		echo "Fasta file ${input} found"
	elif ls *.fna >/dev/null 2>&1
	then
		input=$(ls *fna | sed s/.*\ //)
		echo "Fasta file ${input} found"
	else
		echo "No fasta file found, please check and restart"
	fi
else
	echo "Input fasta: ${input}"
fi

#Make and cd to output directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#Copy and rename files...because of stupid eccentricities of some code
cp ../${input} input.fa
cp ${path2}/fastq/${datatype}/clean.fastq.gz reads.fq.gz

#Run tigmint
echo "Running tigmint"
tigmint-make tigmint-long \
	draft=input \
	reads=reads \
	span=auto \
	G=${genomeSize2/.0/} \
	dist=auto \
	t=${threads}

#Clean some stuff up for downstream analyses
unlink input.cut500.tigmint.fa
rm input.fa input.fa.fai reads.fq.gz

echo "Done"
