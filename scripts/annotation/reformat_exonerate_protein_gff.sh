#Combine and reformat exonerate data into gff

path1=$(pwd | sed s/data.*/scripts/)
a=$(pwd | sed s/.*\\///)

if [ -f target_chunk_1_query_chunk_1 ]
then
	mkdir exonerate_output
	mv target_chunk_*_query_chunk_* exonerate_output
fi

mkdir tmp
cd exonerate_output

for i in target_chunk_*_query_chunk_*
do
	sed '1,2d' tmp/${i}.tmp | grep -v "\-\-\ completed\ exonerate\ analysis" > ../tmp/${i}.tmp
done
cd ..
cat tmp/*tmp > ${a}

perl ${path1}/annotation/reformat_exonerate_protein_gff.pl --input_gff ${a} --output_gff ${a}.gff

rm -R tmp

echo "Done"
