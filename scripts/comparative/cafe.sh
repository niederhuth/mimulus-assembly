#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100GB
#SBATCH --job-name cafe
#SBATCH --output=%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
orthogroup_counts="../orthofinder/orthogroup_counts.tsv"
species_tree="../orthofinder/Results_Sep27/Species_Tree/SpeciesTree_rooted_node_labels.txt"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/phylo/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/phylo/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*data\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*data\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
version=$(ls ${genotype}*-v*.fa | sed s/.*\-v// | sed s/\.fa//)

#
#head -1 ${orthogroup_counts} | cut -f 2- | awk -v OFS="\t" '{print "Desc","Family ID",$0}' > orthogroup_counts.tsv
#sed '1d' ${orthogroup_counts} | awk -v OFS="\t" '{print "(null)",$0}' >> orthogroup_counts.tsv

#
cafe5 \
	-i orthogroup_counts.tsv \
	-t ${species_tree}

echo "Done"