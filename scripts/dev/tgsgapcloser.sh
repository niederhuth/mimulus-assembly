#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=500GB
#SBATCH --job-name tigmint-long
#SBATCH --output=../../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=20 #doesn't seem to want to use more than 6
datatype="ont" #ont or pb


input= #input fasta, if left blank, will look for it in current directory, mutually exclusive with input_dir
input_dir="pilon" #common directory name, e.g. pilon or polca to look for assemblies, eclusive with "input"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/test/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/test/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
assembly=$(pwd | sed s/^.*\\///)
path2=$(pwd | sed s/${genotype}\\/${sample}.*/${genotype}\\/${sample}/)

#Look for fasta file, there can only be one!
if [ -z ${input} ]
then
	echo "No input fasta provided, looking for fasta"
	if ls *.fa >/dev/null 2>&1
	then
		input=$(ls *fa | sed s/.*\ //)
		echo "Fasta file ${input} found"
	elif ls *.fasta >/dev/null 2>&1
	then
		input=$(ls *fasta | sed s/.*\ //)
		echo "Fasta file ${input} found"
	elif ls *.fna >/dev/null 2>&1
	then
		input=$(ls *fna | sed s/.*\ //)
		echo "Fasta file ${input} found"
	else
		echo "No fasta file found, please check and restart"
	fi
else
	echo "Input fasta: ${input}"
fi

#Set options
options="--thread ${threads} --tgstype ${datatype}"
if [ ${racon} = "TRUE" ]
then
	options="${options} --racon"
fi
if [ ${pilon} = "TRUE" ]
then
	options="${options} --ngs ${ngs} --pilon"
fi
if [[ ${racon} = "FALSE" && ${pilon} = "FALSE" ]]
	options="--ne"
fi

#Run tgsgapcloser
echo "Running tgsgapcloser"
tgsgapcloser \
	--scaff ${input} \
	--reads ${} \
	--output tgsgapcloser \
	${options}

echo "Done"

          --ne                             do not error correct. error correct by default.
          or
          --racon     <racon>              the installed racon.
          or
          --ngs       <ngs_reads>          the ngs reads used for pilon
          --pilon     <pilon>              the installed pilon.
          --samtools  <samtools>           the installed samtools.
          --java      <java>               the installed java.
      optional:
          --min_idy   <min_idy>            min_idy for filter reads .
                                           0.3 for ont by default.
                                           0.2 for pb by default.
          --min_match <min_idy>            min match length for filter reads .
                                           300bp for ont by default.
                                           200bp for pb by default.
          --chunk     <chunk_num>          split candidate into chunks to error-correct separately.
          --pilon_mem <t_num>              memory used for pilon , 300G for default.
          --p_round   <pilon_round>        pilon error-correction round num . 3 by default.
          --r_round   <racon_round>        racon error-correction round num . 1 by default.
          --g_check                        gapsize diff check , none by default.

echo "Done"
