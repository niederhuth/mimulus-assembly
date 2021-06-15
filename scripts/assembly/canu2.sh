#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100GB
#SBATCH --job-name canu2
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
correctedErrorRate=0.154 #increased slightly due to lower coverage
useGrid="true" #true or false, if run into issues, maybe try false
gridOptions="--time=168:00:00"
saveReads=True #Keep corrected & trimmed reads
corMinCoverage=0
ovsMemory=62
canuIterationMax=4

#This should match the dataype in the misc/samples.csv file
#Options include:
#"ont" = Raw Nanopore
#"ont-cor" = Corrected Nanopore
#"pac" = raw PacBio
#"pac-cor" = Corrected PacBio
#"hifi" = PacBio HiFi
datatype="ont" 

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/assembly/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/assembly/lib:${LD_LIBRARY_PATH}"

#Set temporary directories for large memory operations
export TMPDIR=$(pwd | sed s/cre-genomes.*//)
export TMP=$(pwd | sed s/cre-genomes.*//)
export TEMP=$(pwd | sed s/cre-genomes.*//)

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/^.*\\///)
path2="canu2"
reads="fastq/${datatype}/clean2.fastq.gz"

#Declare reads
echo "reads: ${reads}"

#Set canu options based on datatype
if [ ${datatype} = "ont" ]
then
	dt="-nanopore"
	echo "Raw Nanopore reads, using canu argument ${dt}"
elif [ ${datatype} = "ont-cor" ]
then
	dt="-nanopore"
	echo "Corrected Nanopore reads, using canu argument ${dt}"
elif [ ${datatype} = "pac" ]
then
	dt="-pacbio"
	echo "Raw PacBio reads, using canu argument ${dt}"
elif [ ${datatype} = "pac-cor" ]
then
	dt="-pacbio"
	echo "Corrected PacBio reads, using canue argument ${dt}"
elif [ ${datatype} = "hifi" ]
then
	dt="-pacbio-hifi"
	echo "PacBio Hifi reads, using canu argument ${dt}"
else
	echo "Do not recognize ${datatype}"
	echo "Please check and resubmit"
fi

#Get genome size estimate
genomeSize=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${sample} \
	-v d=${datatype} \
	'{if ($1 == a && $2 == b && $3 == c && $5 == d) print $9}' \
	${path1}/samples.csv)

#Run canu
if [ -d ${path2} ] 
then
	echo "Previous canu assembly detected, picking up from previous run"
	echo "If you wish to start over, please delete ${path2} and resubmit"
fi	
echo "Beginning canu assembly"
canu \
	-p ${genotype} \
	-d ${path2} \
	genomeSize=${genomeSize} \
	correctedErrorRate=${correctedErrorRate} \
	saveReads=${saveReads} \
	corMinCoverage=${corMinCoverage} \
	useGrid=${useGrid} \
	gridOptions=${gridOptions} \
	ovsMemory=${ovsMemory} \
	canuIterationMax=${canuIterationMax} \
	${dt} $reads

echo "Done"
