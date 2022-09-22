cat final_pseudogenes.gff | sed s/mRNA/pseudo_mRNA/ > tmp
cat renamed.gff tmp | bedtools sort > valid_pseudo_genes.gff
rm tmp
#
awk '$3=="gene" || $3=="pseudogene"' valid_pseudo_genes.gff | cut -f1,3,4,5,9 | \
sed 's/ID\=//' | sed 's/\;.*//' > pseudogene_order.tsv
#
awk '$3=="mRNA" || $3=="pseudo_mRNA"' valid_pseudo_genes.gff | cut -f1,3,4,5,9 | sed 's/ID\=//' | \
sed 's/\;.*//' > pseudogene_mRNA_order.tsv

while read line
do
a=$(echo $line | cut -f1)
b=$(echo $line | cut -f2)
m=$(echo $a | sed 's/.*v2\.0/v2\.0/')
echo $m
if [[ ${m} == "v2.0" ]]
then
c=$(echo $a | sed s/\.v2\.0/\.1\.v2\.0/)
echo "${c} ${b}.1" | tr ' ' '\t' >> tmp
else
echo "${a}.1 ${b}.1" | tr ' ' '\t' >> tmp
fi
done < rename_pseudogene.map

mv tmp rename_pseudogene_mRNA.map
cat rename_pseudogene_mRNA.map rename_pseudogene.map | sort > rename.map


cat final_pseudogenes.gff | while read line
do
	mol=$(echo ${line} | cut -d ' ' -f3)
	id=$(echo ${line} | cut -d ' ' -f9 | sed 's/\;.*//' | sed 's/ID\=//')
	echo ${id}
	if [[ ${mol} == "pseudogene" ]]
	then
		new=$(awk -v a=${id} '{if ($1==a) print $2}' rename.map)
		echo ${line} | awk -v OFS="\t" -v a=${new} -v c=${id} '{print $1,$2,$3,$4,$5,$6,$7,$8,"ID="a";Name="a";putative_pseudogene=TRUE;pseudogene_model="c}' >> renamed.gff
	elif [[ ${mol} == "mRNA" ]]
	then
		new=$(awk -v a=${id} '{if ($1==a) print $2}' rename.map)
		parent=$(echo ${line} | cut -d ' ' -f9 | sed 's/.*Parent\=//' | sed 's/\;.*//')
		new_parent=$(awk -v a=${parent} '{if ($1==a) print $2}' rename.map)
		echo ${line} | awk -v OFS="\t" -v a=${new} -v b=${new_parent} -v c=${id} '{print $1,$2,"mRNA",$4,$5,$6,$7,$8,"ID="a";Name="a";Parent="b";putative_pseudogene=TRUE;pseudogene_model="c}' >> renamed.gff
	elif [[ ${mol} == "exon" ]]
	then
		id2=$(echo ${id} | sed 's/\.exon.*//' | sed 's/\:.*//')
		new=$(awk -v a=${id2} '{if ($1==a) print $2}' rename.map)
		new2=$(echo ${id} | sed s/${id2}/${new}/ | sed 's/\:/\.exon\./')
		parent=$(echo ${line} | cut -d ' ' -f9 | sed 's/.*Parent\=//' | sed 's/\;.*//')
		new_parent=$(awk -v a=${parent} '{if ($1==a) print $2}' rename.map)
		echo ${line} | awk -v OFS="\t" -v a=${new2} -v b=${new_parent} -v c=${id} '{print $1,$2,$3,$4,$5,$6,$7,$8,"ID="a";Parent="b";putative_pseudogene=TRUE;pseudogene_model="c}' >> renamed.gff
	elif [[ ${mol} == "CDS" ]]
	then
		id2=$(echo ${id} | sed 's/\.CDS.*//' | sed 's/\:.*//')
		new=$(awk -v a=${id2} '{if ($1==a) print $2}' rename.map | sed 's/$/\.CDS/')
		parent=$(echo ${line} | cut -d ' ' -f9 | sed 's/.*Parent\=//' | sed 's/\;.*//')
		new_parent=$(awk -v a=${parent} '{if ($1==a) print $2}' rename.map)
		echo ${line} | awk -v OFS="\t" -v a=${new} -v b=${new_parent} -v c=${id} '{print $1,$2,$3,$4,$5,$6,$7,$8,"ID="a";Parent="b";putative_pseudogene=TRUE;pseudogene_model="c}' >> renamed.gff
	elif [[ ${mol} == "five_prime_UTR" ]]
	then
		id2=$(echo ${id} | sed 's/\.five_prime_UTR.*//' | sed 's/\:.*//')
		new=$(awk -v a=${id2} '{if ($1==a) print $2}' rename.map | sed 's/$/\.five_prime_UTR/')
		parent=$(echo ${line} | cut -d ' ' -f9 | sed 's/.*Parent\=//' | sed 's/\;.*//')
		new_parent=$(awk -v a=${parent} '{if ($1==a) print $2}' rename.map)
		echo ${line} | awk -v OFS="\t" -v a=${new} -v b=${new_parent} -v c=${id} '{print $1,$2,$3,$4,$5,$6,$7,$8,"ID="a";Parent="b";putative_pseudogene=TRUE;pseudogene_model="c}' >> renamed.gff
	elif [[ ${mol} == "three_prime_UTR" ]]
	then
		id2=$(echo ${id} | sed 's/\.three_prime_UTR.*//' | sed 's/\:.*//')
		new=$(awk -v a=${id2} '{if ($1==a) print $2}' rename.map | sed 's/$/\.three_prime_UTR/')
		parent=$(echo ${line} | cut -d ' ' -f9 | sed 's/.*Parent\=//' | sed 's/\;.*//')
		new_parent=$(awk -v a=${parent} '{if ($1==a) print $2}' rename.map)
		echo ${line} | awk -v OFS="\t" -v a=${new} -v b=${new_parent} -v c=${id} '{print $1,$2,$3,$4,$5,$6,$7,$8,"ID="a";Parent="b";putative_pseudogene=TRUE;pseudogene_model="c}' >> renamed.gff
	fi
done