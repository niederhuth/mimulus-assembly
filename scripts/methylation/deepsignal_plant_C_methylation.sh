#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=200GB
#SBATCH --job-name deepsignal-plant-C-methylation
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=50
fast5_dir="fast5"
model="$(pwd | sed s/data.*/scripts/)/methylation/ont_models/model.dp2.CNN.arabnrice2-1_120m_R9.4plus_tem.bn13_sn16.both_bilstm.epoch6.ckpt" #Path to model for deepsignal-plant
output_format="tsv" #tsv or bed format

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/deepsignal/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/deepsignal/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
path2=$(pwd | sed s/data.*/data/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/data\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/data\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="methylation"
datatype="ont"
path3="ont_mods"

#Check for and make/cd working directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#Call modfifications with deepsignal-plant
echo "Calling methylation with deepsignal-plant"
CUDA_VISIBLE_DEVICES=-1 deepsignal_plant call_mods \
	--gzip \
	--input_path ${fast5_dir}/ \
	--model_path ${model_path} \
	--result_file ${sample}_ont_C_methylation.tsv.gz \
	--corrected_group RawGenomeCorrected_000 \
	--motifs C \
	--nproc ${threads}

#Call modification frequency and output in tsv format
if [ ${output_format} == "tsv" ]
then
	echo "Calling methylation frequency"
	echo "Outputing in tsv format"
	deepsignal_plant call_freq \
		--nproc ${threads} \
		--sort \
		--gzip \
		--input_path ${sample}_ont_C_methylation.tsv.gz \
		--result_file ${sample}_ont_C_methylation_freqs.tsv.gz
elif [ ${output_format} == "bed" ]
then
	echo "Calling methylation frequency"
	echo "Outputing in bed format"
	deepsignal_plant call_freq \
		--nproc ${threads} \
		--bed \
		--sort \
		--gzip \
		--input_path ${sample}_ont_C_methylation.bed.gz \
		--result_file ${sample}_ont_C_methylation_freqs.bed.gz
else
	echo "Error: output_format must be either 'tsv' or 'bed'"
fi

echo "Done"

