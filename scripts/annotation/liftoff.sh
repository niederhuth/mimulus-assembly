#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --job-name liftoff
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
target_fa="/mnt/gs21/scratch/niederhu/mimulus-assembly/data/Mguttatus/S1/ref/S1-v1.fa" #target genome fasta to map gff files to, if left blank, look in current directory
target_gff="/mnt/gs21/scratch/niederhu/mimulus-assembly/data/Mguttatus/S1/ref/annotations/S1-v1.gff" #gff file for target genome, only necessary if gffcompare=TRUE

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/liftoff/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/liftoff/lib:${LD_LIBRARY_PATH}"

#Set temporary directories for large memory operations
export TMPDIR=$(pwd | sed s/data.*/data/)
export TMP=$(pwd | sed s/data.*/data/)
export TEMP=$(pwd | sed s/data.*/data/)

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="annotation"
datatype="liftoff"
path2=$(pwd | sed s/data.*/data/)
path3="liftoff"

#Look for fasta file, there can only be one!
if [ -z ${target_fa} ]
then
	echo "No input fasta provided, looking for fasta"
	if ls *.fa >/dev/null 2>&1
	then
		fasta=$(ls *fa | sed s/.*\ //)
		echo "Fasta file ${target_fa} found"
	elif ls *.fasta >/dev/null 2>&1
	then
		fasta=$(ls *fasta | sed s/.*\ //)
		echo "Fasta file ${target_fa} found"
	elif ls *.fna >/dev/null 2>&1
	then
		fasta=$(ls *fna | sed s/.*\ //)
		echo "Fasta file ${target_fa} found"
	else
		echo "No fasta file found, please check and restart"
	fi
else
	echo "Input fasta: ${target_fa}"
fi

#Check for and make/cd working directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#Get list of genomes
genomes=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${condition} \
	-v e=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $4 == d && $5 == e) print $7}' \
	${path1}/samples.csv)

#Loop over genomes and liftover & process the annotations
for i in ${genomes}
do
	mkdir ${i}
	cd ${i}
	version=$(ls ${path2}/${i/_*/}/${i/*_/}/ref/${i/*_/}-v*.fa | sed s/.*\-v// | sed s/\.fa//)
	ref_fa=$(ls ${path2}/${i/_*/}/${i/*_/}/ref/${i/*_/}-v${version}.fa)
	ref_gff=$(ls ${path2}/${i/_*/}/${i/*_/}/ref/annotations/${i/*_/}-v${version}.gff)
	chroms="${path1}/annotation/${i}_${species}_${genotype}_chr_mapping.csv"
	ignore_transcripts="${path1}/annotation/${i}_ignore_transcripts.txt"

	#Run liftoff
	mkdir liftoff
	cd liftoff
	#Filter the ref_gff
	fgrep -v -f ${ignore_transcripts} ${ref_gff} > ${i}.gff
	echo "Lifting over annotations from ${i/_*/} ${i/*_/} to ${species} ${genotype} with liftoff"
	liftoff \
		-cds \
		-polish \
		-o mapped.gff \
		-g ${i}.gff \
		-chroms ${chroms} \
		${target_fa} \
		${ref_fa}

	#Change out of directory
	cd ../

	#Classify annotations with gffcompare
	if [ ! -z ${target_gff} ]
	then
		echo "Running gffcompare"
		mkdir gffcompare
		cd gffcompare
		gffcompare -r ${target_gff} ../liftoff/mapped.gff_polished

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
		cut -f5 | sed 's/q1\://' | sed 's/|.*//' > ../newgenes
		#Get the transcript name
		awk '$4=="p" || $4=="s" || $4 =="u" || $4=="x" || $4=="y"' gffcmp.tracking | \
		cut -f5 | tr '|' '\t' | cut -f2 >../newtranscripts
		cd ../

		#Subset gff to keep new_genes
		echo "Subsetting the new annotations"
		cat newgenes newtranscripts | sort | uniq > newannots
		fgrep -f newannots liftoff/mapped.gff_polished | sort | uniq | bedtools sort > newannots.gff
	else
		cp liftoff/mapped.gff_polished newannots.gff
	fi

	#Get genes with a valid ORF otherwise classify as new_pseudogenes
	grep valid_ORFs=1 newannots.gff | cut -f9 | sed 's/\;.*//' | sed 's/ID\=//' > new_valid_genes
	grep valid_ORF=True newannots.gff | cut -f9 | sed 's/\;.*//' | sed 's/ID\=//' > new_valid_transcripts
	cat new_valid_genes new_valid_transcripts > new_valid_annots
	grep valid_ORFs=0 newannots.gff | cut -f9 | sed 's/\;.*//' | sed 's/ID\=//' > new_pseudogenes
	grep valid_ORF=False newannots.gff | cut -f9 | sed 's/\;.*//' | sed 's/ID\=//' > new_pseudogenes_transcripts
	cat new_pseudogenes new_pseudogenes_transcripts > new_pseudogenes_annots
	echo "$(wc -l new_valid_genes | sed s/\ new_valid_genes//) putative genes with intact ORFs found"
	echo "$(wc -l new_pseudogenes | sed s/\ new_pseudogenes//) putative pseudogenes found"

	#Subset gff of valid genes
	fgrep -f new_valid_annots newannots.gff | sort | uniq | bedtools sort > new_valid_genes.gff

	#Subset gff of new_pseudogenes
	#Change column 3 from gene to pseudogene
	#Add attribute pseudo_gene=TRUE
	fgrep -f new_pseudogenes_annots newannots.gff | sort | uniq | bedtools sort \
	awk -v OFS="\t" '{if ($3=="gene") print $1,$2,"pseudogene",$4,$5,$6,$7,$8,$9";putative_pseudogene=TRUE"; 
					else print$0";putative_pseudogene=TRUE"}' > new_pseudogenes.gff

	#Combine the gff files
	if [ ! -z ${target_gff} ]
	then
		cat ${target_gff} newannots.gff | bedtools sort > all_annot.gff
		#If it is first genome, we will combine with all the genes from the target gff
		#Otherwise, we will use only the final valid genes from the previous round
		if [ ${i} = ${genomes/\ */} ]
		then
			cat ${target_gff} new_valid_genes.gff | bedtools sort > all_valid.gff
		else
			cat ${final_valid} new_valid_genes.gff | bedtools sort > all_valid.gff
		fi
		cat ${final_pseudo} new_pseudogenes.gff | bedtools sort > all_pseudogenes.gff
	else
		cat newannots.gff | bedtools sort > all_annot.gff
		cat new_valid_genes.gff | bedtools sort > all_valid.gff
		cat new_pseudogenes.gff | bedtools sort > all_pseudogenes.gff
	fi

	#Get the order of genes for renaming
	awk '$3=="gene"' all_valid.gff | cut -f1,2,4,5,9 | sed 's/ID\=//' | sed 's/\;.*//' > valid_gene_order.tsv
	awk '$3=="mRNA"' all_valid.gff | cut -f1,2,4,5,9 | sed 's/ID\=//' | sed 's/\;.*//' > valid_mRNA_order.tsv

	#Change the target_gff
	target_gff="$(pwd)/all_annot.gff"
	#Set final outputs
	final_valid="$(pwd)/all_valid.gff"
	final_gene_order="$(pwd)/valid_gene_order.tsv"
	final_mRNA_order="$(pwd)/valid_mRNA_order.tsv"
	final_pseudo="$(pwd)/all_pseudogenes.gff"

	#Change out of the directory
	cd ../
done

#Copy the final files
cp ${final_valid} final_valid_genes.gff
cp ${final_gene_order} final_valid_gene_order.tsv
cp ${final_mRNA_order} final_valid_mRNA_order.tsv
cp ${final_pseudo} final_pseudogenes.gff

echo "Done"






