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
prefilter=prequal #prequal filter seqs prior to alignment using specified tool, if blank, no prefiltering
align_unfiltered=TRUE #TRUE/FALSE applies if prefiltered data, if TRUE, will align both filtered & unfiltered
maxiterate=1000 #mafft maximum number of iterations
mode="linsi" #linsi/einsi/ginsi, how to run mafft alignments
datatype=proteins-primary
op=1.53 #mafft gap opening penalty, default: 1.53
ep=0.0 #mafft offset (works like gap extension penalty), default: 0.0
leavegappyregion=FALSE #TRUE/FALSE use --leavegappyregion with mafft
prot2cds=TRUE #Convert the protein alignment to a CDS alignment
trim=all #gblocks/trimal/all trim alignments with X, all will use all trim methods, if blank, no trimming

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
if [ ${leavegappyregion} = TRUE ]
then
	settings="${settings} --leavegappyregion"
fi

#Loop over each gene in gene of interest and perform alignments
sed '1d' ${path1}/goi.csv | while read line
do
	name=$(echo ${line} | cut -d ',' -f1)
	echo "Working on ${name}"
	mkdir ${name}/alignments
	seqs="${name}/seqs/${name}-${datatype}"
	alns="${name}/alignments/${name}-${datatype}"
	#Remove stop "*" from seqs
	sed -i s/\\*$// ${seqs}.fa
	if [ ${prefilter} = "prequal" ]
	then
		echo "Prefiltering the sequences using PREQUAL"
		prequal ${seqs}.fa
		filter=".filtered"
	else
		filter=
	fi
	#Run mafft
	if [ ${prefilter} = "prequal" ]
	then
	echo "Running mafft on ${name} filtered sequences"
		mafft \
			${settings} \
			${seqs}.fa${filter} > ${alns}${filtered}.fas
	fi
	if [[ ! -z ${prefilter} || ${align_unfiltered} = TRUE ]]
	then
		echo "Running mafft on ${name} unfiltered sequences"
		mafft \
			${settings} \
			${seqs}.fa > ${alns}.fas
	fi
	#convert to cds nucleotide alignment
	if [ ${prot2cds} = TRUE ]
	then
		echo "Converting protein alignments to CDS alignments"
		if [ ${prefilter} = "prequal" ]
		then
			   \
				${alns}${filtered}.fas \
				${seqs/proteins/cds}.fa \
				-output fasta > ${alns/proteins/cds}${filtered}.fas

		fi
		if [[ ! -z ${prefilter} || ${align_unfiltered} = TRUE ]]
		then
			pal2nal.pl \
				${alns}.fas \
				${seqs/proteins/cds}.fa \
				-output fasta > ${alns/proteins/cds}.fas
		fi
	fi
	#Trim alignments
	#Trim with gblocks
	if [[ ${trim_method} = "gblocks" || ${trim_method} = "all" ]]
	then
		echo "Trimming the alignment with Gblocks"
		Gblocks ${alns/proteins/cds}.fas \
			-t=c -b3=8 -b4=10 -b5=h -s=y -p=t -e=.gblocks
		#sed -i s/\ //g cds.fas.gb
	fi
	#Trim with trimAL
	if [[ ${trim_method} = "trimal" || ${trim_method} = "all" ]]
	then
		echo "Trimming the alignment with trimAL"
		trimal \
			-in ${alns/proteins/cds}.fas \
			-out ${alns/proteins/cds}.fas.trimal
	fi	

	#axt format
	#perl ${scripts}/phylognetics/pl/parseFastaIntoAXT.pl cds.fas
	#perl $SCRIPTS/parseFastaIntoAXT.pl cds.fas.gb
done


