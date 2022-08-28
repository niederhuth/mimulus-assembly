#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500GB
#SBATCH --job-name liftoff
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables


#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/liftoff/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/liftoff/lib:${LD_LIBRARY_PATH}"

#Set temporary directories for large memory operations
export TMPDIR=$(pwd | sed s/data.*/data/)
export TMP=$(pwd | sed s/data.*/data/)
export TEMP=$(pwd | sed s/data.*/data/)

#Run liftoff
echo ""
liftoff \
	-polish \
	-o mapped \
	-u unmapped \
	-g ${gff} \
	-chroms ${chroms} \
	${target} \
	${reference}

#Get genes with a valid ORF otherwise classify as pseudogenes
grep valid_ORFs=1 mapped_polished | cut -f9 | sed s/\;.*// | sed s/ID\=// > validgenes
grep valid_ORFs=0 mapped_polished | cut -f9 | sed s/\;.*// | sed s/ID\=// > pseudogenes

#Subset gff of valid genes
fgrep -f validgenes mapped_polished > validgenes.gff

#Subset gff of pseudogenes
fgrep -f pseudogenes mapped_polished > pseudogenes.gff

#Classify annotations with gffcompare
gffcompare -r ${ref_gff} validgenes.gff

#= - Discard - complete, exact match of intron chain
#c - Prob Discard - contained in reference (intron compatible)
#e - Prob Discard - single exon transfrag partially covering an intron, possible pre-mRNA fragment
#i - Prop Keep - fully contained within a reference intron
#j - Discard - multi-exon with at least one junction match
#k - Prob Discard - containment of reference (reverse containment)
#m - Prob Discard - retained intron(s), all introns matched or retained
#n - Prob Discard - retained intron(s), not all introns matched/covered
#o - Prob Discard - other same strand overlap with reference exons
#p - Keep - possible polymerase run-on (no actual overlap)
#s - Keep - intron match on the opposite strand (likely amapping error)
#x - Keep - exonic overlap on the opposite strand (like o or e but on the opposite strand)
#y - Keep - contains a reference within its intron(s)
#u - Keep - none of the above (unknown, intergenic)
#Keep p,s,u,x,y
awk '$4=="p" || $4=="s" || $4 =="u" || $4=="x" || $4=="y"' gffcmp.tracking | cut -f5 | sed s/q1\:// | sed s/\|.*// > newgenes

#Subset gff to keep new_genes
fgrep -f newgenes mapped_polished | bedtools sort > newgenes.gff

#
awk '$3=="gene"' newgenes.gff > tmp
#Combine the files
cat S1-genes tmp | bedtools sort | cut -f9 | sed s/\;.*// | sed s/ID\=// > combined_list





rm tmp2
count=0
cat combined_list | while read line
do
	genome=$(echo ${line} | sed s/Mg// | sed s/_.*//)
	if [ ${genome} == "S1" ]
	then
		gene=${line}
		chrom=$(echo ${line} | sed s/.*_// | sed s/g.*//)
		count=0
	elif [ ${genome} == "L1" ]
	then
		chrom2=$(echo ${line} | sed s/.*_// | sed s/g.*//)
		if [ ${chrom2} == ${chrom} ]
		then
			count=$(expr ${count} + 1)
			echo ${line} | awk -v a=${gene} -v b=${count} -v OFS="\t" '{print $0',a,b}'' >> tmp2
		fi	
	fi
done

for i in $(cut -f2 tmp2 | sort | uniq)
do
	grep ${i} tmp2 > tmp3
	
done












