#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=300GB
#SBATCH --job-name orthofinder
#SBATCH --output=%x-%j.SLURMout

cd ${PBS_O_WORKDIR}
export PATH="${HOME}/miniconda3/envs/orthofinder/bin:${PATH}"
export LD_LIBRARY_PATH="${HOME}/miniconda3/envs/orthofinder/lib:${LD_LIBRARY_PATH}"

#Set Variables
species=$(cut -d ',' -f1 ../../misc/genomes.csv | sed '1d' | tr '\n' ' ')
threads=100
threads2=6

#Run OrthoFinder
echo "Copying sequence files"
mkdir seqs
for i in ${species}
do
	cp ../${i}/ref/mcscanx/${i}-protein.fa  seqs/${i}.fa
done

echo "Running OrthoFinder"
orthofinder \
	-t ${threads} \
	-a ${threads2} \
	-M dendroblast \
	-S diamond_ultra_sens \
	-I 1.3 \
	-y \
	-s ../../misc/SpeciesTree_rooted.txt \
	-o orthofinder \
	-f seqs/

echo "Create orthogroup list"
if [ -f orthogroup_list.tsv ]
then
	rm orthogroup_list.tsv
fi
sed '1d' orthofinder/*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv | while read line
do
	og=$(echo ${line} | cut -d ' ' -f1 | sed s/\\:$//) 
	echo ${line} | cut -d ' ' -f4- | tr -d $'\r' | sed s/,//g | tr ' ' ',' | sed s/,$// | tr ',' '\n' | awk -v OFS='\t' -v a=${og} '{print $0,a}' | sed s/gene_BVRB/gene:BVRB/ >> orthogroup_list.tsv 
done

echo "Get orthogroups for each species & gene"
for i in $species
do
	echo ${i}
	cut -f2 ../${i}/ref/mcscanx/${i}.gff > tmp
	fgrep -f tmp orthogroup_list.tsv > ../${i}/ref/mcscanx/${i}_orthogroups.tsv
	cut -f1 ../${i}/ref/mcscanx/${i}_orthogroups.tsv > tmp2
	a=1
	fgrep -w -v -f tmp2 tmp | while read line
	do
		echo $line ${i}_${a} | tr ' ' '\t' >> ../${i}/ref/mcscanx/${i}_orthogroups.tsv
		a=$(expr ${a} + 1)
	done
done
rm tmp tmp2

echo "Done"

