
#
awk '$38 != $39' orthogroup_counts.tsv > cnv_orthogroup_counts.tsv

#
fgrep -f cnv_orthogroup_list.txt orthogroups.csv | cut -d ',' -f1,38 > L1_cnv_orthogroups.csv
fgrep -f cnv_orthogroup_list.txt orthogroups.csv | cut -d ',' -f1,39 > S1_cnv_orthogroups.csv

#
cut -d ',' -f2 L1_cnv_orthogroups.csv | tr ' ' '\n' | sort | uniq | sed '/^$/d' > L1_cnv_orthologs.txt
cut -d ',' -f2 S1_cnv_orthogroups.csv | tr ' ' '\n' | sort | uniq | sed '/^$/d' > S1_cnv_orthologs.txt

#
fgrep -v -f ../genespace/results/L1_1x1_syntelogs.txt L1_cnv_orthologs.txt > tmp
mv tmp L1_cnv_orthologs.txt
fgrep -v -f ../genespace/results/S1_1x1_syntelogs.txt S1_cnv_orthologs.txt > tmp
mv tmp S1_cnv_orthologs.txt

#
awk '$2=="Mguttatus_L1"' Results_Sep27/Gene_Duplication_Events/Duplications.tsv | cut -f6,7 | tr '\t' '\n' | sed s/\ //g | tr ',' '\n' | sort | uniq | sed s/Mguttatus_L1_// | sed /^$/d > L1_gene_duplications.txt
awk '$2=="Mguttatus_S1"' Results_Sep27/Gene_Duplication_Events/Duplications.tsv | cut -f6,7 | tr '\t' '\n' | sed s/\ //g | tr ',' '\n' | sort | uniq | sed s/Mguttatus_S1_// | sed /^$/d > S1_gene_duplications.txt