#Set variables
threads=20 #threads for mpileup
threads2=5 #threads for call
max_alleles=2 #Max allele count for filtering, recommend 2 for now
maf=0.1 #Minor allele frequency, assuming this data is from F2s or F3s, so keeping high
max_non_ref_af=0.9 #Maximum non-ref allele frequencyi, use to eliminate allels without ref allele
max_missing=0.2 #Max number of samples can be missing calls

name=$1

vcftools \
	--vcf ${name}.vcf \
	--out ${name} \
	--recode \
	--remove-indels \
	--max-alleles ${max_alleles} \
	--maf ${maf} \
	--max-non-ref-af ${max_non_ref_af} \
	--max-missing ${max_missing}
#Modify genotypes
echo "Reformatting genotypes"
cat ${name}.recode.vcf | \
sed s/\\.\\/\\.\\:/NN\:/ | \
sed s/0\\/0\\:/AA\:/g | \
sed s/0\\/1\\:/AB\:/g | \
sed s/1\\/1\\:/BB\:/g > ${name}.modified.recode.vcf
if [ -s ${name}.modified.recode.vcf ]
then
	vcf=${name}.modified.recode.vcf
	#Make g_files directory if not present
	if [ ! -d g_files ]
	then
		mkdir g_files
	fi
	#Get sample number from vcf
	ncol=$(grep "#CHROM" ${name}.modified.recode.vcf | awk '{print NF; exit}')
	#we start at column 10, where first sample is
	a=10
	echo "Outputing individual sample genotypes"
	#Loop over each sample column and output a g_gile
	until [ ${a} -gt ${ncol} ]
	do
		#Get sample for that column
		sample=$(grep "#CHROM" ${vcf} | cut -f${a} | sed s/bam_files\\/// | sed s/.bam//)
		#Cut the propter columsn and output g_file
		cut -f1,2,${a} ${vcf} | grep -v \# | sed s/A:.*/A/ | sed s/B:.*/A/ | sed s/NN:.*/NN/ |\
		awk -v OFS="\t" -v x=${sample} '{print x,$0}' > g_files/${sample}_genotype.txt
		#Add 1 to the column
		a=$(expr ${a} + 1)
	done
fi
