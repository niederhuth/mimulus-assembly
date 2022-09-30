#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=50GB
#SBATCH --job-name=update-gene-functions
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
threads=50
old_version="1"
new_version="1.2"
old_proteins="../../final/pseudomolecule/annotations/${old_version}/*protiens.fa"
new_proteins="../../final/pseudomolecule/annotations/${new_version}/*protiens.fa"

#Change to current directory
cd ${PBS_O_WORKDIR}
#Export paths to conda
export PATH="${conda}/envs/gene-function/bin:${PATH}"
export LD_LIBRARY_PATH="${conda}/envs/gene-function/lib:${LD_LIBRARY_PATH}"

#Set temporary directories for large memory operations
export TMPDIR=$(pwd)
export TMP=$(pwd)
export TEMP=$(pwd)

#The following shouldn't need to be changed, but should set automatically
path1=$(pwd | sed s/data.*/misc/)
path2=$(pwd | sed s/data.*/scripts/)
species=$(pwd | sed s/^.*\\/data\\/// | sed s/\\/.*//)
genotype=$(pwd | sed s/.*\\/${species}\\/// | sed s/\\/.*//)
sample=$(pwd | sed s/.*${species}\\/${genotype}\\/// | sed s/\\/.*//)
output="${genotype}-v${new_version}"
path3="update_gene_functions"



#Make & cd to directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#Find the blast results
blast="$(pwd | sed s/data.*/data/)/${species}/${genotype}/comparative/diamond_blastp/Athaliana_Athaliana/${species}_${genotype}-Athaliana_Athaliana_orthogroup_filtered.m8"

#Copy over the proteins
cp ${old_proteins} ./
cp ${new_proteins} ./

#map new ids onto the v1 roteins
${conda}/envs/maker/bin/map_fasta_ids ../liftoff/rename.map ${genotype}-v${old_version}-proteins.fa

#Get gene names for old and new proteins
grep \> ${genotype}-v${old_version}-proteins.fa | sed s/\>// > old_proteins
grep \> ${genotype}-v${new_version}-proteins.fa | sed s/\>// > new_proteins 

#Run interproscan
echo "Running interproscan"
${path2}/annotation/interproscan/interproscan.sh \
	--cpu ${threads} \
	-appl pfam \
	-goterms \
	-pa \
	-dp \
	-iprlookup \
	-t p \
	-f TSV \
	-i ${output}-proteins.fa \
	-o ${output}.iprscan

#Download and format Arabidopsis TAIR10 functional descriptions
echo "Downloading and formatting Arabidopsis TAIR10 functional descriptions"
wget -q https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_functional_descriptions
perl -e  'while (my $line = <>){ my @elems = split "\t", $line; if($elems[2] ne "") {print "$elems[0]\t$elems[2]\n"}}' \
TAIR10_functional_descriptions > TAIR10_short_functional_descriptions.txt

#Download Arabidopsis GO terms
echo "Downloading Arabidopsis GO terms"
wget -q https://www.arabidopsis.org/download_files/GO_and_PO_Annotations/Gene_Ontology_Annotations/gene_association.tair.gz
gunzip gene_association.tair.gz

#Create header for output file
echo "Transcript Locus Arabidopsis_blast_hit Arabidopsis_GO_terms PFAM_hits PFAM_GO_terms Combined_Arabidopsis_PFAM_GO_terms Short_functional_description" | \
tr ' ' '\t' > ${output}-functional-annotations.tsv
#Loop over each gene and format data
cat new_proteins | while read line
do
	echo ${line}
	#Handle the Arabidopsis BLAST
	AtID=$(grep ${line} ${blast} | sort -r -k12 | head -1 | cut -f2)
	if [[ ! -z ${AtID} ]]
	then
		#Get the Arabidopsis Description
		AtDesc=$(grep ${AtID} TAIR10_short_functional_descriptions.txt | cut -f2)
		if [[ ! -z ${AtDesc} ]]
		then
			AtDesc=$(echo "Arabidopsis BLAST: ${AtDesc}" | cut -f2 | tr ' ' ';')
		else
			AtDesc=NA
		fi
		#Get Arabidopsis GO terms
		AtGO=$(grep ${AtID} gene_association.tair | cut -f5 | tr '|' '\n' | sort | uniq | tr '\n' '|' | \
			sed s/\|$//)
		if [ -z ${AtGO} ]
		then
			AtGO=NA
		fi
	else
		AtID=NA
		AtDesc=NA
		AtGO=NA
	fi
	#Get the Pfam domains
	grep ${line} ${output}.iprscan > tmp
	#If tmp is not empty
	if [ -s tmp ]
	then
		#List the PfamIDs
		PfamID=$(cut -f5 tmp | sort | uniq | tr '\n' '|' | sed s/\|$//)
		#Get the Pfam Descriptions
		PfamDesc=$(echo "PFAM: $(cut -f6 tmp | tr '\n' ',')" | tr ' ' ';' | sed s/\,$//)
		if [ -z ${PfamDesc} ]
		then
			PfamDesc=NA
		fi
		#List the PfamGO terms
		PfamGO=$(cut -f14 tmp | tr '|' '\n' | sort | uniq | grep -v "-" | tr '\n' '|' | sed s/\|$//)
		if [ -z ${PfamGO} ]
		then
			PfamGO=NA
		fi
	else
		PfamID=NA
		PfamDesc=NA
		PfamGO=NA
	fi
	if [[ ${AtDesc} != "NA" ]]
	then
		FuncDesc=${AtDesc}
	else
		if [[ ${PfamDesc} != "NA" ]]
		then
			FuncDesc=${PfamDesc}
		else
			if [[ ! -z $(grep ${line} old_proteins | cut -d ' ' -f5 | awk -v FS="|" '$3 > 0') ]]
			then
				FuncDesc="Expressed;gene;of;unknown;function"
			else
				FuncDesc="Hypothetical;gene;of;unknown;function"
			fi
		fi
	fi
	#Combine the GO term sets
	combinedGO=$(echo "${AtGO}|${PfamGO}" | tr '|' '\n' | grep -v NA | sort | uniq | tr '\n' '|' | sed s/\|$//)
	#Output the results
	echo "${line} ${line/\.*/} ${AtID} ${AtGO} ${PfamID} ${PfamGO} ${combinedGO} ${FuncDesc}" | \
	tr ' ' '\t' | tr ';' ' ' >> ${output}-functional-annotations.tsv
	#Remove tmp file
	rm tmp
done
