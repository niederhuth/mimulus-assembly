#!/bin/bash --login
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20GB
#SBATCH --job-name mafft
#SBATCH --output=%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
goi= #path to a csv of genes of interest, if blank, will look for goi.csv in misc dir
threads=20 #number of threads to use with mafft
maxiterate=1000 #maximum number of iterations
mode="linsi" #linsi/einsi/ginsi, how to run mafft alignments
datatype=proteins-primary
op=1.53 #Gap opening penalty, default: 1.53
ep=0.0 #Offset (works like gap extension penalty), default: 0.0

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/phylo/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/phylo/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
path2=$(pwd | sed s/data.*/scripts/)

#Set mafft settings
settings="--thread ${threads} --maxiterate ${maxiterate} --op ${op} --ep ${ep}"
if [ ${mode} = "linsi" ]
then
	echo "Running mafft in linsi mode"
	settings="${settings} --localpair"
elif [ ${mode} = "einsi" ]
then
	echo "Running mafft in einsi mode"
	settings="${settings} --genafpair"
elif [ ${mode} = "ginsi" ]
	echo "Running mafft in ginsi mode"
	settings="${settings} --globalpair"
fi

while read line
do
	name=$(echo ${line} | cut -d ' ' -f1)
	echo "Working on ${name}"
	#Run mafft
	echo "Running mafft on ${name}"
	mafft \
		${settings} \
		${i}/${i}-${datatype}.fa > ${i}/${i}-${datatype}.fas

	#convert to cds nucleotide alignment
	echo "Converting ${i}/${i}-${datatype}.fas to CDS alignment"
	pal2nal.pl \
		${i}/${i}-${datatype}.fas \
		${i}/${i}-${datatype/proteins/cds}.fa \
		-output fasta > ${i}/${i}-${datatype/proteins/cds}.fas
	pal2nal.pl \
		${i}/${i}-${datatype}.fas \
		${i}/${i}-${datatype/proteins/cds}.fa \
		-output paml > ${i}/${i}-${datatype/proteins/cds}.phy

	#gblocks
	Gblocks cds.fas -t=c -b3=8 -b4=10 -b5=h -s=y -p=t -e=.gb
	sed -i s/\ //g cds.fas.gb

	#axt format
	perl ${scripts}/phylognetics/pl/parseFastaIntoAXT.pl cds.fas
fi


