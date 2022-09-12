#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=300GB
#SBATCH --job-name orthofinder
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=100 #sequence search threads
threads2=10 #analysis threads
inflation=1.3 #Inflation parameter, default 1.5
seq_searc_program="diamond_ultra_sens" #blast/diamond/diamond_ultra_sens/blast_gz/mmseqs/blast_nucl sequence search program 
msa_program="mafft" #mafft/muscle
msa=TRUE #TRUE/FALSE use multiple sequence alignment
tree_program="fasttree" #fasttree/raxml/raxml-ng/iqtree only applies if msa=TRUE
tree= #path to a user provided tree, if left blank, orthofinder will generate its own
is_DNA=FALSE #TRUE/FALSE, sequences are DNA

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/orthofinder/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/orthofinder/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species="comparative"
genotype="orthofinder"
sample="orthofinder"
condition="orthofinder"
datatype="proteins-primary"
path2=$(pwd | sed s/data.*/data/)
path3="orthofinder"

#Get list of genomes
genomes=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${condition} \
	-v e=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $4 == d && $5 == e) print $7}' \
	${path1}/samples.csv)

#Copy over input sequences
echo "Copying sequence files"
if [[ ! -d seqs ]]
then
	mkdir ${datatype}
fi
for i in ${genomes}
do
	#Find input sequences
	species2=$(echo ${i} | sed s/_.*//)
	genotype2=$(echo ${i} | sed s/${species2}_// | sed s/_.*//)
	path4="${path2}/${species2}/${genotype2}/ref/annotations"
	version=$(ls ${path4}/${genotype2}-v*-${datatype}.fa | sed s/.*\-v// | sed s/.fa//)
	input_seqs="${path4}/${genotype2}-v${version2}-${datatype}.fa" ${datatype}/${species2}_${genotype2}.fa
done

#Set msa options
if [ ${msa} = TRUE ]
then
	settings="-M msa -A ${msa_program} -T ${tree_program}"
elif [ ${msa} = FALSE ]
then
	settings="-M dendroblast -A ${msa_program}"
else
	echo "msa must be set to either TRUE or FALSE"
fi
#Set user provided tree
if [[ ! -z ${tree} ]]
then
	echo "User provided rooted tree supplied"
	settings="${settings} -s ${tree}"
fi
#Are sequences DNA?
if [ ${is_DNA} = TRUE ]
then
	settings="${settings} -d"
fi

#Run Orthofinder
echo "Running OrthoFinder"
orthofinder \
	-t ${threads} \
	-a ${threads2} \
	${settings} \
	-S ${seq_searc_program} \
	-I ${inflation}\
	-y \
	-o ${path2} \
	-f ${datatype}

echo "Done"

