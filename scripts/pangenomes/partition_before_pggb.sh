#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20GB
#SBATCH --job-name partition-before-pggb
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20
haplotype_number=2 #Number of haplotypes/genomes
pidentity=95 #Mercent identity for mapping/alignment [default: 90]
segment_length=1000  #Segment length for mapping [default: 5000]
min_match_length=19  #Filter exact matches below this length [default: 19]
vcf= #Make a VCF against "ref" decomposing variants, SPEC = REF:DELIM[:LEN]
stats=TRUE #TRUE/FALSE generate statistics of the seqwish and smoothxg graph
multiqc=TRUE #TRUE/FALSE generate MultiQC report of graphs' statistics and visualizations, runs odgi stats

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/pangenome/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/pangenome/lib:$LD_LIBRARY_PATH"
#Export path to UDOCKER_DIR. All images will be downloaded and installed here
export UDOCKER_DIR=${conda}/envs/pangenome/udocker
#Export path to UDOCKER_CONTAINERS. All containers will be saved there
if [[ ! -d containers ]]
then
	mkdir containers
fi
export UDOCKER_CONTAINERS=$(pwd)/containers

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
path2=$(pwd | sed s/data.*/data/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample="pangenome"
condition="pggb"
datatype="genome"
path3="pggb"

#Make & cd to working directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#Get list of genomes
genomes=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${condition} \
	-v e=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $4 == d && $5 == e) print $7}' \
	${path1}/samples.csv)

if [[ -d sequences ]]
then
	mkdir sequences
fi

#Change sequence names to PanSN-spec: https://github.com/pangenome/PanSN-spec
if [ -f sequences/${species}.fa.gz ]
then
	echo ""
else
	echo "Copying and renaming input fastas"
	for seq in ${genomes}
	do
		echo "${seq} ... Starting"
		#Get the version #
		version=$(ls ${path2}/${seq/_*/}/${seq/*_/}/ref/${seq/*_/}-v*.fa | sed s/.*\-v// | sed s/\.fa//)
		#Set input path & output path
		infa=${path2}/${seq/_*/}/${seq/*_/}/ref/${seq/*_/}-v${version}.fa
		outfa=sequences/${seq/*_/}-v${version}.fa.gz
		#Haplotype number...need to work on this for phased genomes
		haplotype=1
		#Modify sequence names and output
		sed s/\>/\>${seq/*/}\#${haplotype}\#/ ${infa} | bgzip --threads ${threads} > ${outfa}
	done
	#Combine the fasta files
	echo "Combining input fastas"
	cat sequences/*-v*fa.gz | bgzip --threads ${threads} > ${species}.fa.gz
	#Index the combined fasta
	samtools faidx ${species}.fa.gz
fi

#Set settings
settings="-t ${threads} -n ${haplotype_number} -i ${species}.fa.gz -o partitions"
if [[ ! -z ${pidentity} ]]
then
	settings="${settings} -p ${pidentity}"
fi
if [[ ! -z ${segment_length} ]]
then
	settings="${settings} -s ${segment_length}"
fi
if [[ ! -z ${min_match_length} ]]
then
	settings="${settings} -k ${min_match_length}"
fi
if [ ${stats} = "TRUE" ]
then
	settings="${settings} -S"
fi
if [ ${multiqc} = "TRUE" ]
then
	settings="${settings} --multiqc"
fi
if [[ ! -z ${vcf} ]]
then
	settings="${settings} -V ${vcf}"
fi

#Check for docker container, if it doesnt exist, create it
if [[ ! -d containers/pggb ]]
then
	udocker create --name=pggb ghcr.io/pangenome/pggb:20230321161006ffd010
fi

#Run partition-before-pggb
echo "Running partition-before-pggb"
echo "Settings: ${settings}"
udocker run \
	--volume=$(pwd):/data \
	partition-before-pggb ${settings}

#Make output directory for pggb runs for invidual paritions
if [[ ! -d partition_scripts ]]
then
	mkdir partition_scripts
fi

#Parse out 
echo "Parsing out pggb commands for each community"
cat partitions/*yml | while read line
do
	if [[ ${line} =~ ^"pggb -i".* ]]
	then
		community=$(echo ${line} | sed s/.*community/community/ | sed s/\.fa.*//)
		echo ${line} | sed s/\-o\ partitions/\-o\ community_graphs/ > partition_scripts/${community}.sh
	fi
done

echo "Done"
