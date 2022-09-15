#!/bin/bash --login
#SBATCH --time=3:59:00
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
trim=FALSE #TRUE/FALSE Trim the alignments, if TRUE, will generate files for both trimmed & untrimmed
trim_method= #gblocks/trimal

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/phylo/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/phylo/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
path2=$(pwd | sed s/data.*/scripts/)

#Set path to gene of interest list
if [ -z ${goi} ]
then
	goi=${path1}/goi.csv
fi

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
then
	echo "Running mafft in ginsi mode"
	settings="${settings} --globalpair"
fi

sed '1d' ${path1}/goi.csv | while read line
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
	#Trim alignments
	if [ ${trim} = TRUE ]
	then
		#Trim with gblocks
		if [ ${trim_method} = "gblocks" ]
		then
			echo "Trimming the alignment with gblocks"
			Gblocks ${i}/${i}-${datatype/proteins/cds}.fas \
				-t=c -b3=8 -b4=10 -b5=h -s=y -p=t -e=.gb
			#sed -i s/\ //g cds.fas.gb
		#Trim with trimAL
		elif [ ${trim_method} = "trimal" ]
		then
			echo "Trimming the alignment with trimAL"
			trimal \
				-in ${i}/${i}-${datatype/proteins/cds}.fas \
				-out ${i}/${i}-${datatype/proteins/cds}_trimal.fas
		#Trim with
		elif [ ${trim_method} = "guidance2" ]
		then
			echo "Trimming the alignment with "
		#Trim with
		elif [ ${trim_method} = "something" ]
		then
			echo "Trimming the alignment with "
		fi
	fi	

	#axt format
	#perl ${scripts}/phylognetics/pl/parseFastaIntoAXT.pl cds.fas
	#perl $SCRIPTS/parseFastaIntoAXT.pl cds.fas.gb
done


