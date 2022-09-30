#!/bin/bash --login
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=50GB
#SBATCH --job-name=new_gene_functions
#SBATCH --output=../job_reports/%x-%j.SLURMout

#Set this variable to the path to wherever you have conda installed
conda="${HOME}/miniconda3"

#Set variables
theads=50
genotype="S1"
output="S1-v1.2"
blast="../../../comparative/diamond_blastp/Athaliana_Athaliana/Mguttatus_S1-Athaliana_Athaliana_orthogroup_filtered.m8"

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
path3="new_gene_functions"

#Make & cd to directory
if [ -d ${path3} ]
then
	cd ${path3}
else
	mkdir ${path3}
	cd ${path3}
fi

#Copy over the proteins
cp ../../final/pseudomolecule/annotations/v1/${genotype}-v1-proteins.fa ./
cp ../../final/pseudomolecule/annotations/v1.2/${genotype}-v1.2-proteins.fa ./

#map new ids onto the v1 roteins
${conda}/envs/maker/bin/map_fasta_ids ../liftoff/rename.map ${genotype}-v1-proteins.fa

#
grep \> ${genotype}-v1-proteins.fa | sed s/\>// > old_proteins
grep \> ${genotype}-v1.2-proteins.fa | sed s/\>// > new_proteins 

#Run interproscan
echo "Running interproscan"
${path2}/annotation/interproscan/interproscan.sh \
	-cpu ${threads} \
	-appl pfam \
	-goterms \
	-pa \
	-dp \
	-iprlookup \
	-t p \
	-f TSV \
	-i new.fa \
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
echo "Transcript Locus Arabidopsis_blast_hit Arabidopsis_GO_terms PFAM_hits PFAM_GO_terms Short_functional_description" | \
tr ' ' '\t' > ${output}-functional-annotations.tsv
#Loop over each gene and format data
cat new_proteins | while read line
do
	#Handle the Arabidopsis BLAST
	AtID=$(grep ${line} ${blast} | sort -r -k12 | head -1 | cut -f2)
	if [[ ! -z ${AtID} ]]
	then
		#Get the Arabidopsis Description
		AtDesc=$(echo "Arabidopsis BLAST: $(grep ${AtID} TAIR10_short_functional_descriptions.txt)" | \
			tr ' ' ';')
		#Get Arabidopsis GO terms
		AtGO=$(grep ${AtID} gene_association.tair | cut -f5 | tr '|' '\n' | sort | uniq | tr '\n' '|' | \
			sed s/\|$//)
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
		PfamDesc=$("PFAM: $(cut -f6 | tr '\n' ',')" | tr ' ' ';')
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
		PfamGO=NA
		PfamID=NA
	fi
	if [[ ${AtDesc} != "NA" ]]
	then
		FuncDesc=${AtDesc}
	else
		if [[ ${PfamDesc} != "NA" ]]
		then
			FuncDesc=${AtDesc}
		else
			if [[ $(grep ${line} old_proteins | cut -d ' ' -f5 | tr '|' '\t' | cut -f3) -gt 0 ]]
			then
				FuncDesc="Expressed gene of unknown function"
			else
				FuncDesc="Hypothetical gene of unknown function"
			fi
		fi
	fi
	#Output the results
	echo "${line} ${line/\.*/} ${AtID} ${AtGO} ${PfamID} ${PfamGO} ${FuncDesc}" | \
	tr ' ' '\t' | tr ';' ' ' >> ${output}-functional-annotations.tsv
	#Remove tmp file
	rm tmp
done
