#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=50GB
#SBATCH --job-name stringtie
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=40
bam="../SRrna/Aligned.sortedByCoord.out.bam" #Can be space separated list
read_type="rf" #fr: fr-secondstrand, rf: fr-firststrand, lr: longread, mix: mix lr & shortread
#To mimic --conservative set min_multi_exon_reads=1.5, min_iso_frac=0.05, trim=FALSE
trim=TRUE #use coverage based trimming of transcript ends
min_multi_exon_reads="1.5" #min reads per bp cov for multi-exon transcript default: 1
min_single_exon_reads="4.75" #min reads per bp cov for single-exon transcript default: 4.75, 1.5 for long-reads
min_iso_frac="0.05" #minimum isoform fraction default: 0.01
max_gap="50" #maximum gap allowed between read mappings default: 50, set to 0 for long-reads
min_transcript_len="200" #minimum assembled transcript length default: 200
min_anchor_len="10" #minimum anchor length for junctions default: 10
min_junc_cov="1" #minimum junction coverage default: 1
frac_multi_hit="0.95" #fraction of bundle allowed to be covered by multi-hit reads default: 1
LR_splice_window="25" #window around possibly erroneous splice sites from long reads default: 25

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/transcript-assembly/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/transcript-assembly/lib:${LD_LIBRARY_PATH}"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="annotation"
assembly=$(pwd | sed s/^.*\\///)
path2="stringtie"

#Add various settings
settings="-v -p ${threads}"
if [ ${trim} = FALSE ]
then
	settings="${settings} -t"
fi
#Adjust setting based on read type
if [ ${read_type} = "rf" ]
then
	settings="${settings} --rf"
elif [ ${read_type} = "fr" ]
then
	settings="${settings} --fr"
elif [${read_type} = "lr" ]
then
	settings="${settings} -L -E ${LR_splice_window}"
elif [${read_type} = "mix" ]
then
	settings="${settings} --mix -E ${LR_splice_window}"
fi

#Make output dir and cd
if [ ! -d ${path2} ]
then
	mkdir ${path2}
fi
cd ${path2}

#Run stringtie
echo "Running stringtie"
stringtie ${bam} \
	${settings} \
	-c ${min_multi_exon_reads} \
	-s ${min_single_exon_reads} \
	-f ${min_iso_frac} \
	-g ${max_gap} \
	-m ${min_transcript_len} \
	-a ${min_anchor_len} \
	-j ${min_junc_cov} \
	-M ${frac_multi_hit} \
	-o ${path2} \
	-l stringtie 

echo "Done"

