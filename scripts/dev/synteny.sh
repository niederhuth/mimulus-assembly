#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --job-name busco2allmaps_synt
#SBATCH --output=%x-%j.SLURMout

#
i="IM62"

#Change to current directory
cd ${PBS_O_WORKDIR}

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
assembly=$(pwd | sed s/^.*\\///)
path2=$(pwd | sed s/${genotype}\\/${sample}\\/.*/${genotype}\\/${sample}/)
path3=allmaps


#Make and cd to output directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

grep -v \# ../../busco/*/*sol*/full_table.tsv | \
awk -v OFS="\t" '{if ($2=="Complete") print $3,$4,$5,$1,100,$6}' > S1.bed

path4=$(pwd | sed s/${genotype}\\/${sample}\\/.*/${i}/)
grep -v \# ${path4}/ref/busco/*/*sol*/full_table.tsv | \
awk -v OFS="\t" '{if ($2=="Complete") print $3,$4,$5,$1,100,$6}' | grep  > REF.bed

cut -f4 S1.bed > tmp3
cut -f4 REF.bed > tmp4
fgrep -w -f tmp3 tmp4 | awk -v OFS="\t" '{print $1,$1,100}' > REF.S1.1x1.anchors

python -m jcvi.assembly.syntenypath bed --switch --scale=100000 REF.S1.1x1.anchors -o synteny.bed

#python -m jcvi.assembly.allmaps mergebed synteny.bed -o allmaps.bed

python -m jcvi.assembly.allmaps mergebed ../primers.bed synteny.bed -o allmaps.bed

python -m jcvi.assembly.allmaps path allmaps.bed ../fasta_other/input.fa --notex