#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --job-name rename
#SBATCH --output=%x-%j.SLURMout

species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)

#Hack Script for renaming
mkdir tmp
cd tmp
#Rename chr genes
cut -f1 ../final_valid_gene_order.tsv | sort | uniq | grep "^chr" | while read line
do
	echo ${line}
	count=10
	if [[ $(echo ${line} | wc -c) -gt 5 ]]
	then
		chr=$(echo ${line} | sed s/chr//)
	else 
		chr=$(echo ${line} | sed s/chr/0/)
	fi
	awk -v a=${line} -v OFS="\t" '$1==a' ../final_valid_gene_order.tsv | while read genes
	do
		if [[ $(expr 6 - $(echo ${count} | wc -c)) = 0 ]]
		then
			zeros=
		elif [[ $(expr 6 - $(echo ${count} | wc -c)) = 1 ]]
		then
			zeros="0"
		elif [[ $(expr 6 - $(echo ${count} | wc -c)) = 2 ]]
		then 
			zeros="00"
		elif [[ $(expr 6 - $(echo ${count} | wc -c)) = 3 ]]
		then 
			zeros="000"
		elif [[ $(expr 6 - $(echo ${count} | wc -c)) = 4 ]]
		then 
			zeros="0000"
		fi
		name="Mg${genotype}_${chr}g${zeros}${count}"
		count=$(expr ${count} + 10)
		echo "${genes} ${name}" | tr ' ' '\t' >> ${line}_rename_genes.tsv
	done
done
#Rename contig genes
count=10
chr="UN"
grep -v "^chr" ../final_valid_gene_order.tsv | while read genes
do
	if [[ $(expr 6 - $(echo ${count} | wc -c)) = 0 ]]
	then
		zeros=
	elif [[ $(expr 6 - $(echo ${count} | wc -c)) = 1 ]]
	then
		zeros="0"
	elif [[ $(expr 6 - $(echo ${count} | wc -c)) = 2 ]]
	then 
		zeros="00"
	elif [[ $(expr 6 - $(echo ${count} | wc -c)) = 3 ]]
	then 
		zeros="000"
	elif [[ $(expr 6 - $(echo ${count} | wc -c)) = 4 ]]
	then 
		zeros="0000"
	fi
	name="Mg${genotype}_${chr}g${zeros}${count}"
	count=$(expr ${count} + 10)
	echo "${genes} ${name}" | tr ' ' '\t' >> chrUN_rename_genes.tsv
done
#Rename chr transcripts
cut -f1 ../final_valid_mRNA_order.tsv | sort | uniq | grep "^chr" | while read line
do
	echo ${line}
	count=10
	if [[ $(echo ${line} | wc -c) -gt 5 ]]
	then
		chr=$(echo ${line} | sed s/chr//)
	else 
		chr=$(echo ${line} | sed s/chr/0/)
	fi
	awk -v a=${line} -v OFS="\t" '$1==a' ../final_valid_mRNA_order.tsv | while read genes
	do
		if [[ $(expr 6 - $(echo ${count} | wc -c)) = 0 ]]
		then
			zeros=
		elif [[ $(expr 6 - $(echo ${count} | wc -c)) = 1 ]]
		then
			zeros="0"
		elif [[ $(expr 6 - $(echo ${count} | wc -c)) = 2 ]]
		then 
			zeros="00"
		elif [[ $(expr 6 - $(echo ${count} | wc -c)) = 3 ]]
		then 
			zeros="000"
		elif [[ $(expr 6 - $(echo ${count} | wc -c)) = 4 ]]
		then 
			zeros="0000"
		fi
		name="Mg${genotype}_${chr}g${zeros}${count}.1"
		count=$(expr ${count} + 10)
		echo "${genes} ${name}" | tr ' ' '\t' >> ${line}_rename_transcripts.tsv
	done
done
#Rename contig transcripts
count=10
chr="UN"
grep -v "^chr" ../final_valid_mRNA_order.tsv | while read genes
do
	if [[ $(expr 6 - $(echo ${count} | wc -c)) = 0 ]]
	then
		zeros=
	elif [[ $(expr 6 - $(echo ${count} | wc -c)) = 1 ]]
	then
		zeros="0"
	elif [[ $(expr 6 - $(echo ${count} | wc -c)) = 2 ]]
	then 
		zeros="00"
	elif [[ $(expr 6 - $(echo ${count} | wc -c)) = 3 ]]
	then 
		zeros="000"
	elif [[ $(expr 6 - $(echo ${count} | wc -c)) = 4 ]]
	then 
		zeros="0000"
	fi
	name="Mg${genotype}_${chr}g${zeros}${count}.1"
	count=$(expr ${count} + 10)
	echo "${genes} ${name}" | tr ' ' '\t' >> chrUN_rename_transcripts.tsv
done
cat chr*_rename_genes.tsv chr*_rename_transcripts.tsv | cut -f5,6 | sort > ../rename.map
cd ..

#Rename the gff
cat final_valid_genes.gff | while read line
do
	mol=$(echo ${line} | cut -d ' ' -f3)
	id=$(echo ${line} | cut -d ' ' -f9 | sed 's/\;.*//' | sed 's/ID\=//')
	echo ${id}
	if [[ ${mol} == "gene" ]]
	then
		new=$(awk -v a=${id} '{if ($1==a) print $2}' rename.map)
		echo ${line} | awk -v OFS="\t" -v a=${new} '{print $1,$2,$3,$4,$5,$6,$7,$8,"ID="a";Name="a}' >> renamed.gff
	elif [[ ${mol} == "mRNA" ]]
	then
		new=$(awk -v a=${id} '{if ($1==a) print $2}' rename.map)
		parent=$(echo ${line} | cut -d ' ' -f9 | sed 's/.*Parent\=//' | sed 's/\;.*//')
		new_parent=$(awk -v a=${parent} '{if ($1==a) print $2}' rename.map)
		echo ${line} | awk -v OFS="\t" -v a=${new} -v b=${new_parent} '{print $1,$2,$3,$4,$5,$6,$7,$8,"ID="a";Name="a";Parent="b}' >> renamed.gff
	elif [[ ${mol} == "exon" ]]
	then
		id2=$(echo ${id} | sed 's/\.exon.*//' | sed 's/\:.*//')
		new=$(awk -v a=${id2} '{if ($1==a) print $2}' rename.map)
		new2=$(echo ${id} | sed s/${id2}/${new}/ | sed 's/\:/\.exon\./')
		parent=$(echo ${line} | cut -d ' ' -f9 | sed 's/.*Parent\=//' | sed 's/\;.*//')
		new_parent=$(awk -v a=${parent} '{if ($1==a) print $2}' rename.map)
		echo ${line} | awk -v OFS="\t" -v a=${new2} -v b=${new_parent} '{print $1,$2,$3,$4,$5,$6,$7,$8,"ID="a";Parent="b}' >> renamed.gff
	elif [[ ${mol} == "CDS" ]]
	then
		id2=$(echo ${id} | sed 's/\.CDS.*//' | sed 's/\:.*//')
		new=$(awk -v a=${id2} '{if ($1==a) print $2}' rename.map | sed 's/$/\.CDS/')
		parent=$(echo ${line} | cut -d ' ' -f9 | sed 's/.*Parent\=//' | sed 's/\;.*//')
		new_parent=$(awk -v a=${parent} '{if ($1==a) print $2}' rename.map)
		echo ${line} | awk -v OFS="\t" -v a=${new} -v b=${new_parent} '{print $1,$2,$3,$4,$5,$6,$7,$8,"ID="a";Parent="b}' >> renamed.gff
	elif [[ ${mol} == "five_prime_UTR" ]]
	then
		id2=$(echo ${id} | sed 's/\.five_prime_UTR.*//' | sed 's/\:.*//')
		new=$(awk -v a=${id2} '{if ($1==a) print $2}' rename.map | sed 's/$/\.five_prime_UTR/')
		parent=$(echo ${line} | cut -d ' ' -f9 | sed 's/.*Parent\=//' | sed 's/\;.*//')
		new_parent=$(awk -v a=${parent} '{if ($1==a) print $2}' rename.map)
		echo ${line} | awk -v OFS="\t" -v a=${new} -v b=${new_parent} '{print $1,$2,$3,$4,$5,$6,$7,$8,"ID="a";Parent="b}' >> renamed.gff
	elif [[ ${mol} == "three_prime_UTR" ]]
	then
		id2=$(echo ${id} | sed 's/\.three_prime_UTR.*//' | sed 's/\:.*//')
		new=$(grep -w ${id2} rename.map  | cut -f2 | sed 's/$/\.three_prime_UTR/')
		parent=$(echo ${line} | cut -d ' ' -f9 | sed 's/.*Parent\=//' | sed 's/\;.*//')
		new_parent=$(awk -v a=${parent} '{if ($1==a) print $2}' rename.map)
		echo ${line} | awk -v OFS="\t" -v a=${new} -v b=${new_parent} '{print $1,$2,$3,$4,$5,$6,$7,$8,"ID="a";Parent="b}' >> renamed.gff
	fi
done


