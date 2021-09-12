#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500GB
#SBATCH --job-name protein2genome_test
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=1
chunks=20
minintron=10
maxintron=3000
bestn=5
ryo="no" #">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n"

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
path2=$(pwd | sed s/data.*/data/)
path3=$(pwd | sed s/data.*/scripts/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*${species}\\/${genotype}\\/// | sed s/\\/.*//)
version=$(ls ${sample}-v*.fa | sed s/.*\-v// | sed s/.fa//) 
path4="exonerate_test"

#Look for fasta file, there can only be one!
if ls *.fa >/dev/null 2>&1
then
	fasta=$(ls *fa | sed s/.*\ //)
	echo "Fasta file ${fasta} found"
elif ls *.fasta >/dev/null 2>&1
then
	fasta=$(ls *fasta | sed s/.*\ //)
	echo "Fasta file ${fasta} found"
elif ls *.fna >/dev/null 2>&1
then
	fasta=$(ls *fna | sed s/.*\ //)
	echo "Fasta file ${fasta} found"
else
	echo "No fasta file found, please check and restart"
fi

#Make output directory
if [ -d ${path4} ]
then
	cd ${path4}
else
	mkdir ${path4}
	cd ${path4}
fi

#Get protein sources
proteins=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	'{if ($1 == a && $2 == b && $3 == c) print $4}' \
	${path1}/annotation/annotation_sources.csv)

#Copy & combined protein sources
for i in ${proteins}
do
	a=$(echo $i | cut -d '_' -f1)
	b=$(echo $i | cut -d '_' -f2)
	c=$(echo $i | cut -d '_' -f3)
	input="${path2}/${a}/${b}/${c}/ref/annotations/${c}-v*-protein.fa"
	cat ${input} >> protein_seqs.fa
done

output1="protein2genome.output"
output2="protein2genome.gff"

#Check for previous runs
if ls ${output1}_chunk_* >/dev/null 2>&1
then
	a=$(ls test_* | tail -1 | sed s/.*_//)
	echo "Previous chunks found, picking up starting with chunk ${a}"
else
	a=1
	echo "No previous chunks found, starting from chunk 1"
fi

#Run exonerate over for each chunk
while [ ${a} -le ${chunks} ]
do
	echo "Running exonerate protein2genome chunk ${a} on ${fasta}"
	exonerate \
		--cores ${threads} \
		--model protein2genome \
		--bestn ${bestn} \
		--minintron ${minintron} \
		--maxintron ${maxintron} \
		--querychunkid ${a} \
		--querychunktotal ${chunks} \
		--query protein_seqs.fa \
		--target ../${fasta} \
		--showtargetgff yes \
		--showalignment no \
		--showvulgar no \
		--ryo ${ryo} > ${output1}_chunk_${a}
	a=$(expr ${a} + 1)
done

#Combine chunks
echo "Combining exonerate output chunks"
cat ${output1}_chunk_* > ${output1}

#Reformat gff files
echo "Reformatting exonerate gff"
perl ${path3}/annotation/reformat_exonerate_protein_gff.pl \
	--input_gff ${output1} \
	--output_gff ${output2}

echo "Done"
