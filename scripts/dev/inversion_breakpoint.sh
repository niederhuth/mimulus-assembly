
zcat genespace/results/syntenicBlocks.txt.gz | awk '$1=="L1" && $2=="S1" && $3=="chr8" && $4=="chr8"' | cut -f4,5,6,7,20,21,22 | sort -k3n

zcat genespace/results/syntenicBlocks.txt.gz | awk '$1=="S1" && $2=="L1" && $3=="chr8" && $4=="chr8"' | cut -f4,5,6,7,20,21,22 | sort -k3n


#
extend=50000
fasta="../../../L1/ref/L1-v1.fa"
fai="../../../L1/ref/L1-v1.fa.fai"
ref="../../ref/S1-v1.fa"
threads=10
threads2=4


mkdir regions indexes
cat inversions.bed | while read line
do
	region=$(echo ${line} | tr ' ' '_')
	if [ ! -d indexes/${region/_*/} ]
	then
		mkdir indexes/${region/_*/}
		samtools faidx ${ref} ${region/_*/} > indexes/${region/_*/}/${region/_*/}.fa
		bwa index -p indexes/${region/_*/}/${region/_*/} indexes/${region/_*/}/${region/_*/}.fa
	fi
	dir="regions/${region}"
	mkdir regions/${region}
	echo ${line} | tr ' ' '\t' > ${dir}/${region}.bed
	bedtools slop -b ${extend} -i ${dir}/${region}.bed -g ${fai} > ${dir}/${region}_ext.bed
	bedtools getfasta -bed ${dir}/${region}_ext.bed -fi ${fasta} -fo ${dir}/${region}_ext.fa
	wgsim -e 0 -d 400 -s 50 -N 1000000 -1 150 -2 150 -r 0 -R 0 -X 0 -S 0 ${dir}/${region}_ext.fa ${dir}/${region}_1.fastq ${dir}/${region}_2.fastq
	bwa mem -t ${threads} \
				-R "@RG\tID:${region}\tLB:SIM\tPL:SIM\tSM:${region}\tPU:SIM" \
				-M indexes/${region/_*/}/${region/_*/} ${dir}/${region}_1.fastq ${dir}/${region}_2.fastq | \
				samtools view -@ ${threads2} -bSh | \
				samtools sort -@ ${threads2} > ${dir}/${region}.bam
done
samtools merge combined.bam regions/*/*bam
samtools index combined.bam
