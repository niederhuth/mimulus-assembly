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
threads2=20 #analysis threads
datatype="proteins" #proteins, proteins-primary, etc
input_seqs=proteins #Directory of input sequences, if left blank will copy over from list of genomes in misc/samples.csv
inflation=1.3 #Inflation parameter, default 1.5
seq_search_program=diamond #blast/diamond/diamond_ultra_sens/blast_gz/mmseqs/blast_nucl sequence search program 
msa=TRUE #TRUE/FALSE use multiple sequence alignment
msa_program=mafft #mafft/muscle
tree_program=fasttree #fasttree/raxml/raxml-ng/iqtree only applies if msa=TRUE
tree= #path to a user provided tree, if left blank, orthofinder will generate its own
split_HOGs=TRUE #TRUE/FALSE Split paralogous orthogroups: -y argument in orthofinder
is_DNA=FALSE #TRUE/FALSE, sequences are DNA
start_from=f #f/b/fg/ft : f=full pipeline, b=from previous blast, fg=from previous orthogroup, ft=from previous gene tree  
previous_results= #directory of previous results for restarting orthofinder at later stages

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
path2=$(pwd | sed s/data.*/data/)
path3="orthofinder"

#Get list of genomes
if [ -z ${input_seqs} ]
then
	genomes=$(awk -v FS="," \
		-v a=${species} \
		-v b=${genotype} \
		-v c=${sample} \
		-v d=${condition} \
		-v e=${datatype} \
		'{if ($1 == a && $2 == b && $3 == c && $4 == d && $5 == e) print $7}' \
		${path1}/samples.csv)
	echo "Genomes: ${genomes}"

fi

#Copy over input sequences
if [ -z ${input_seqs} ]
then
	if [ ${start_from} = "f" ]
	then
		echo "Copying sequence files"
		if [[ ! -d ${datatype} ]]
		then
			mkdir ${datatype}
		fi
		for i in ${genomes}
		do
			#Find input sequences
			species2=$(echo ${i} | sed s/_.*//)
			genotype2=$(echo ${i} | sed s/${species2}_//)
			path4="${path2}/${species2}/${genotype2}/ref/annotations"
			version=$(ls ${path4}/${genotype2}-v*-${datatype}.fa | sed s/.*\-v// | sed s/\-${datatype}.fa//)
			cp ${path4}/${genotype2}-v${version}-${datatype}.fa ${datatype}/${species2}_${genotype2}.fa
		done
	fi
else
	echo 'Input sequences already provided in: ${input_seqs}'
fi

#Set msa options
if [ ${msa} = TRUE ]
then
	settings="-M msa -A ${msa_program} -T ${tree_program}"
elif [ ${msa} = FALSE ]
then
	settings="-M dendroblast"
else
	echo "msa must be set to either TRUE or FALSE"
fi
#Set user provided tree
if [[ ! -z ${tree} ]]
then
	echo "User provided rooted tree supplied"
	settings="${settings} -s ${tree}"
fi
#Split HOGs
if [ ${split_HOGs} = TRUE ]
then
	settings="${settings} -y"
fi
#Are sequences DNA?
if [ ${is_DNA} = TRUE ]
then
	settings="${settings} -d"
fi
#Set starting point
if [ ${start_from} = "f" ]
then
	if [ -z ${input_seqs} ]
	then
		settings="${settings} -f ${datatype} -o ${path3}"
	else
		settings="${settings} -f ${input_seqs} -o ${path3}"
	fi
elif [ ${start_from} = "b" ]
then
	settings="${settings} -b ${previous_results}"	
elif [ ${start_from} = "fg" ]
then
	settings="${settings} -fg ${previous_results}"
elif [ ${start_from} = "ft" ]
then
	settings="${settings} -ft ${previous_results}"
fi

#Run Orthofinder
echo "Running OrthoFinder"
orthofinder \
	-t ${threads} \
	-a ${threads2} \
	-S ${seq_search_program} \
	-I ${inflation}\
	${settings} 	

echo "Done"

