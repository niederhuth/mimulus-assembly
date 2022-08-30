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
gene_list=""
chroms="../chroms.txt"
target_fa="../S1-v1.fa"
ref_fa="../L1-v1.fa"
target_gff="../S1-v1.gff"
ref_gff="../L1-v1.gff"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/liftoff/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/liftoff/lib:${LD_LIBRARY_PATH}"

#Set temporary directories for large memory operations
export TMPDIR=$(pwd | sed s/data.*/data/)
export TMP=$(pwd | sed s/data.*/data/)
export TEMP=$(pwd | sed s/data.*/data/)

#Subset gff
if [ ${gene_list} ]
then
	fgrep -f ${gene_list} ${ref_gff} > search.gff
else
	cat ${ref_gff} > search.gff
fi

#Run liftoff
liftoff \
	-cds \
	-polish \
	-o mapped.gff \
	-g search.gff \
	-chroms ${chroms} \
	${target_fa} \
	${ref_fa}

#Get genes with a valid ORF otherwise classify as pseudogenes
grep valid_ORFs=1 mapped.gff_polished | cut -f9 | sed 's/\;.*//' | sed 's/ID\=//' > validgenes
grep valid_ORFs=0 mapped.gff_polished | cut -f9 | sed 's/\;.*//' | sed 's/ID\=//' > pseudogenes

#Subset gff of valid genes
fgrep -f validgenes mapped.gff_polished > validgenes.gff

#Subset gff of pseudogenes
fgrep -f pseudogenes mapped.gff_polished > pseudogenes.gff

#Classify annotations with gffcompare
gffcompare -r ${target_gff} validgenes.gff

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
awk '$4=="p" || $4=="s" || $4 =="u" || $4=="x" || $4=="y"' gffcmp.tracking | \
cut -f5 | sed 's/q1\://' | sed 's/|.*//' > newgenes

#Subset gff to keep new_genes
fgrep -f newgenes mapped.gff_polished | bedtools sort > newgenes.gff

#
cat ${target_gff} newgenes.gff | bedtools sort > combined.gff

#
awk '$3=="gene"' combined.gff | cut -f1,2,4,5,9 | sed 's/ID\=//' | sed 's/\;.*//' > combined_order








