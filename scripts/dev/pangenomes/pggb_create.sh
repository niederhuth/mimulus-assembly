#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20GB
#SBATCH --job-name mafft
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/pangenomes/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/pangenomes/lib:$LD_LIBRARY_PATH"

#
wfmash Mguttatus.fa.gz -p 90 -n 6 -t 4 -m > Mguttatus.mapping.paf
#
paf2net.py -p Mguttatus.mapping.paf
#
net2communities.py \
	-e Mguttatus.mapping.paf.edges.list.txt \
	-w Mguttatus.mapping.paf.edges.weights.txt \
	-n Mguttatus.mapping.paf.vertices.id2name.txt
#
seq 0 14 | while read i; do
    chromosomes=$(cat scerevisiae7.mapping.paf.edges.weights.txt.community.$i.txt | cut -f 3 -d '#' | sort | uniq | tr '\n' ' ');
    echo "community $i --> $chromosomes";
done