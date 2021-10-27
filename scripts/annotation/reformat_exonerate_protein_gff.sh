#Combine and reformat exonerate data into gff

path1=$(pwd | sed s/data.*/scripts/)
a=$(pwd | sed s/.*\\///)

mkdir exonerate_output
mkdir tmp
for i in target_chunk_*_query_chunk_*
do
	perl ${path1}/annotation/reformat_exonerate_protein_gff.pl --input_gff ${i} --output_gff tmp/${i}.tmp
	sed '1,2d' tmp/${i}.tmp | grep -v "\-\-\ completed\ exonerate\ analysis" > tmp/${i}.tmp2
done
mv target_chunk_*_query_chunk_* exonerate_output
cat tmp/*tmp2 > ${a}.gff
rm -R tmp

echo "Done"
