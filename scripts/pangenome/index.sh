#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=60GB
#SBATCH --job-name index
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20
worflow="giraffe" #map, mpmap, rpvg, giraffe
ref_fasta="$(pwd | sed s/pangenome.*//)/L1/ref/L1-v1.fa"
vcf="data/S1_SVs.vcf.gz"
gff=
tmp_dir=

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/pangenome/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/pangenome/lib:$LD_LIBRARY_PATH"

#Make output directory
if [[ ! -d ${worflow} ]]
then
	mkdir ${worflow}
fi

#Check for vcf directory, make if it doesn't exist
if [[ ! -d vcf && ! -z ${vcf} ]]
then
	mkdir vcf
fi
#Check for gff directory, make if it doesn't exist
if [[ ! -d gff && ! -z ${gff} ]]
then
	mkdir gff
fi

#Merge vcf files
if [ -f vcf/merged.vcf.gz ]
then
	echo "vcf/merged.vcf.gz already exists"
	echo "To redo merging of vcf files please delete vcf/merged.vcf.gz and resubmit"
else
	if [ $(echo ${vcf} | tr ' ' 'n' | wc -l) = 1 ]
	then
		cp ${vcf} vcf/merged.vcf.gz
	else
		echo "Merging vcf files"
		bcftools merge \
			--threads ${threads} \
			-O z \
			-o vcf/merged.vcf.gz \
			-m none \
			-0 ${vcf}
	fi
	echo "Indexing vcf file"
	bcftools index vcf/merged.vcf.gz
fi

#Set arguments
arguments="-t ${threads} --workflow ${worflow} --prefix ${worflow}/index --ref-fasta ${ref_fasta}"
#Set vcf
if [[ ! -z ${vcf} ]]
then
	arguments="${arguments} --vcf vcf/merged.vcf.gz"
fi
#Add gff/gtf if given
if [[ ! -z ${gff} ]]
then
	arguments="${arguments} --tx-gff gff/merged.gff"
fi
#Set temporary directory
if [ -z ${tmp_dir} ]
then
	arguments="${arguments} --tmp-dir $(pwd)"
else
	arguments="${arguments} --tmp-dir ${tmp_dir}"
fi

#Run vg autoindex
echo "Running vg autoindex --workflow ${worflow}"
vg autoindex ${arguments} 

echo "Done"
