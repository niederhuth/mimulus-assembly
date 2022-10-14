bedtools slop -b 10000 -i inversions.bed -g ../../../L1/ref/L1-v1.fa.fai > inversions_+10000.bed
bedtools getfasta -bed inversions_+3000.bed -fi ../../../L1/ref/L1-v1.fa -fo inversions_+10000.fa
wgsim -e 0 -d 300 -s 50 -N 30000 -1 100 -2 100 -r 0 -R 0 -X 0 -S 0 inversions_+10000.fa inversions_1.fastq inversions_2.fastq
bowtie2 -x S1 -1 inversions_1.fastq -2 inversions_2.fastq -S inversions.sam
samtools sort inversions.sam > inversions.bam
samtools index inversions.bam