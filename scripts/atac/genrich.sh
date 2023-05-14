#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50GB
#SBATCH --job-name genrich
#SBATCH --output=job_reports/%x-%j.SLURMout

#Set variables
filter_chrs="chrM,chrC"
remove_PCR_dups=TRUE
atac_mode=TRUE
auc=10 #Minimum AUC for a peak (default 200.0)
max_dist=100 #Maximum distance between significant sites (default 100)
min_length=0 #Minimum length of a peak (def. 0)
sig_cutoff=0.01 #Maximum p-value/q-value (default 0.01)
qvalue=FALSE #Use q-value rather than p-value
skip_peak_calling=FALSE #Skip peak calling and generate log file only
recall_peaks=TRUE #Recall peaks from log file

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"
#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/atac/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/atac/lib:${LD_LIBRARY_PATH}"

#Other variables, these should not have to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*data\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*data\\/${species}\\/${genotype}\\/// | sed s/\\/.*//)
datatype="atac"
path2=$(pwd | sed s/data.*/data/)/${species}/${genotype}

#Creat output directory & cd to it
if [ -d ${datatype} ]
then 
	cd ${datatype}
else
	mkdir ${datatype}
	cd ${datatype}
fi

#Get a list of the condition groups
groups=$(awk -v FS="," \
	-v a=${species} \
	-v b=${genotype} \
	-v c=${datatype} \
	'{if ($1 == a && $2 == b && $5 == c) print $4}' \
	${path1}/samples.csv | sort | uniq | sed '/^$/d' | grep -v nDNA)

#Loop over each reference genome and run genrich
for condition in ${groups}
do
	#Make directory for that condition and cd to it
	if [ -d ${condition} ]
	then
		cd ${condition}
	else
		mkdir ${condition}
		cd ${condition}
	fi

	#Get replicate samples
	replicates=$(awk -v FS="," \
		-v a=${species} \
		-v b=${genotype} \
		-v c=${condition} \
		-v d=${datatype} \
		'{if ($1 == a && $2 == b && $4 == c && $5 == d) print $3}' \
		${path1}/samples.csv | sort | uniq | sed '/^$/d')

	#Get list of genomes
	genomes=$(awk -v FS="," \
		-v a=${species} \
		-v b=${genotype} \
		-v c=${condition} \
		-v d=${datatype} \
		'{if ($1 == a && $2 == b && $4 == c && $5 == d) print $7}' \
		${path1}/samples.csv | sort | uniq | sed '/^$/d')

	#Run genrich for each assembly mapping
	for assembly in ${genomes}
	do
		#Get ref genome info
		path3="$(pwd | sed s/data\\/.*/data/)/${assembly/_*}/${assembly/*_/}"
		version=$(ls ${path3}/ref/${assembly/*_/}-v*.fa | sed s/.*\-v// | sed s/.fa//)
		
		#Set output file
		narrowPeak="${condition}_ref_${assembly}-v${version}_auc_${auc}.narrowPeak"
		log_file="${condition}_ref_${assembly}-v${version}.log"
		reps_log="${condition}_ref_${assembly}-v${version}_reps.log"
		
		#Make a list of input bam files
		bam_files=""
		for sample in ${replicates}
		do
			bam_files="${path2}/${sample}/${datatype}/${sample}_ref_${assembly}-v${version}.bam,${bam_files}"
			nDNA="${path2}/nDNA/atac/nDNA_ref_${assembly}-v${version}.bam,${nDNA}"
		done

		#Set various options for genrich
		genrich_options="-z -v"
		#Recall peaks?
		if [ ${recall_peaks} = "TRUE" ]
		then
			#If recall_peaks = TRUE, set input files
			genrich_options="${genrich_options} -P -f ${log_file}.gz"
		else
			#If recall_peaks = FALSE, set input files
			genrich_options="${genrich_options} -t ${bam_files} -c ${nDNA} -f ${log_file} -k ${reps_log}"
		fi
		#Skip peak calling?
		if [ ${skip_peak_calling} = "TRUE" ]
		then
			genrich_options="${genrich_options} -X"
		else
			genrich_options="${genrich_options} -a ${auc} -g ${max_dist} -l ${min_length} -o ${narrowPeak}"
		fi
		#atac mode
		if [ ${atac_mode} = "TRUE" ] 
		then
			genrich_options="${genrich_options} -j"
		fi
		#Set q-value or p-value
		if [ ${qvalue} = "TRUE" ] 
		then
			genrich_options="${genrich_options} -q ${sig_cutoff}"
		else
			#If qvalue=FALSE, use p-value instead
			genrich_options="${genrich_options} -p ${sig_cutoff}"
		fi
		#Add chromosomes for filtering
		if [[ ! -f $[filter_chrs] ]]
		then
			genrich_options="${genrich_options} -e ${filter_chrs}"
		fi
		#Remove PCR duplicates
		if [ ${remove_PCR_dups} = "TRUE" ] 
		then
			pcr_dups_out="${condition}_ref_${assembly}-v${version}_pcr_duplicates.tsv"
			genrich_options="${genrich_options} -r -R ${pcr_dups_out}"
		fi

		#If log file already exists and recall_peaks is false, skip
		if [[ -f ${log_file}.gz && ${recall_peaks} = "FALSE" ]]
		then
			echo "Output files already exist. Skipping"
			echo "To rerun this step, delete ${log_file}.gz, ${narrowPeak}.gz, ${reps_log}.gz and resubmit"
		#Otherwise run Genrich pipeline
		else
			echo "Running Genrich for ${condition} against reference ${assembly}-v${version}"
			Genrich ${genrich_options}
		fi
	done
	cd ../
done	

echo "Done"
