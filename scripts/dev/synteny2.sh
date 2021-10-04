
ref="IM62"
target="test"
transcripts="cds.fa"

for i in ${target} ${ref}
do
	#map transcripts to fasta
	minimap2 -x splice ${i}.fa ${transcripts} > ${i}.paf 

	#convert to bed file
	awk -v OFS="\t" '{print $6,$8,$9,$1,100,$5}' ${i}.paf > ${i}_aln.bed

	#retain only unique alignments
	cut -f4 ${i}_aln.bed | sort | uniq -c | sed 's/^ *//' | \
	awk -v FS=" " '{if ($1 == 1) print $2}' > ${i}_unique.bed
done

cut -f1 ${ref}.bed > tmp
fgrep -f tmp ${target}.bed | cut -f1 | awk -v OFS="\t" '{print $1,$1,100}' > ${ref}.${target}.1x1.anchors

