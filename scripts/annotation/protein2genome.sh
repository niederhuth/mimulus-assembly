#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=500GB
#SBATCH --job-name protein2genome
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20
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
sample=${genotype}
version=$(ls ${genotype}-v*.fa | sed s/.*\-v// | sed s/.fa//) 
path4="exonerate"

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

#Run EDTA
for i in ${proteins}
do
	source_species=$(echo $i | cut -d '_' -f1)
	source_genotype=$(echo $i | cut -d '_' -f2)
	input="${path2}/${source_species}/${source_genotype}/ref/annotations/${source_genotype}-v*-protein.fa"
	cat ${input} >> protein_seqs.fa
done

output1="protein2genome.output"
output2="protein2genome.gff"

a=1
while [ ${a} -le ${chunks} ]
do
	echo "Running exonerate protein2genome chunk ${a} on ${fasta}"
	exonerate \
		--cores ${threads} \
		--model protein2genome \
		--bestn ${bestn} \
		--minintron ${minintron} \
		--maxintron ${maxintron} \
		--targetchunkid ${a} \
		--targetchunktotal ${chunks} \
		--query protein_seqs.fa \
		--target ../${fasta} \
		--showtargetgff yes \
		--showalignment no \
		--showvulgar no \
		--ryo ${ryo} > ${output1}_chunk_${a}
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
