#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50GB
#SBATCH --job-name allmaps
#SBATCH --output=%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
primers=TRUE #Paired primer sequences for genetic markers
markers=FALSE #Have not implemented this option
synteny=FALSE #Have not implemented this option
optical=FALSE #Have not implemented this option
primer_max_dist=5000 #max distance for primers to be separated

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/allmaps/bin:$PATH"
export LD_LIBRARY_PATH="${conda}/envs/allmaps/lib:$LD_LIBRARY_PATH"

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
condition="assembly"
assembly=$(pwd | sed s/^.*\\///)
path2=$(pwd | sed s/${genotype}\\/${sample}\\/.*/${genotype}\\/${sample}/)
path3=allmaps

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

#Make and cd to output directory
if [ -d ${path3} ]
then
	cd ${path3}
	cp ../${input} input.fa
else
	mkdir ${path3}
	cd ${path3}
	cp ../${input} input.fa
fi

#Aligning primer data with bowtie
if [ ${primers} = "TRUE" ]
then
	#Make bowtie index
	if [ -s bowtie_index/input.rev.2.ebwt ]
	then
		echo "bowtie index found"
	else
		echo "Building bowtie index"
		mkdir bowtie_index
		bowtie-build input.fa bowtie_index/input
	fi

	#Align primers
	echo "Aligning primer data with bowtie"
	bowtie \
		--sam-nohead \
		--no-unal \
		-f \
		-X ${primer_max_dist} \
		-x bowtie_index/input \
		-1 ${path1}/genetic_map/Marker_forward_primers.fa \
		-2 ${path1}/genetic_map/Marker_reverse_primers.fa \
		-S primers.sam

	#Format primers.bed
	for i in $(sed '1d' ${path1}/genetic_map/genetic_map.csv)
	do
		M=$(echo ${i} | cut -d ',' -f1)
		LG=$(echo ${i} | cut -d ',' -f2)
		GP=$(echo ${i} | cut -d ',' -f3)
		F=$(grep ${M}_F primers.sam | cut -f3,4)
		Fchr=$(echo ${F} | cut -d ' ' -f1)
		Fpos=$(echo ${F} | cut -d ' ' -f2)
		R=$(grep ${M}_R primers.sam | cut -f3,4)
		Rchr=$(echo ${R} | cut -d ' ' -f1)
		Rpos=$(echo ${R} | cut -d ' ' -f2)
		if [[ ! -z ${F} && ! -z ${R} ]]
		then
			if [ ${Fchr} == ${Rchr} ]
			then
				echo "Primers ${M} properly paired"
				if [ ${Fpos} -gt ${Rpos} ]
				then
					echo "${Fchr},${Rpos},${Fpos},LG${LG}:${GP}" | tr ',' '\t' >> primers.bed
				elif [ ${Fpos} -lt ${Rpos} ]
				then
					echo "${Fchr},${Fpos},${Rpos},"LG"${LG}:${GP}" | tr ',' '\t' >> primers.bed
				fi
			else
				echo "Primers ${M} improperly paired, skipping"
			fi
		elif [ ! -z ${F} ]
		then
			echo "${M} forward primer only"
			echo "${Fchr},${Fpos},$(expr `${Fpos} + 1`),LG${LG}:${GP}" | tr ',' '\t' >> primers.bed
		elif [ ! -z ${R} ]
		then
			echo "${M} reverse primer only"
			echo "${Rchr},${Fpos},$(expr `${Rpos} + 1`),LG${LG}:${GP}" | tr ',' '\t' >> primers.bed
		else
			echo "Primers ${M} missing data" >> missing_primers.txt
		fi
	done
	position_data="primers.bed ${position_data} "
fi

#Merge files and create weights file
echo "Merging position data files"
python -m jcvi.assembly.allmaps mergebed ${position_data} -o allmaps.bed

#Run allmaps path
echo "Running allmaps path"
python -m jcvi.assembly.allmaps path allmaps.bed input.fa

#Cleanup
mkdir fasta_other
mv allmaps.chr.fasta fasta_other
mv allmaps.unplaced.fasta fasta_other
mv input.fa fasta_other

echo "Done"




