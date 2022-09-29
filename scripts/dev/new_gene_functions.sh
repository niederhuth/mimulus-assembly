genotype="S1"
output="S1-v1.2"

mkdir v1
mv * v1
mkdir v1.2
cd v1.2

#Copy over the proteins
cp ../../../final/pseudomolecule/annotations/v1/${genotype}-v1-proteins.fa ./
cp ../../../final/pseudomolecule/annotations/v1.2/${genotype}-v1.2-proteins.fa ./

#map new ids onto the v1 roteins
map_fasta_ids ../../liftoff/rename.map ${genotype}-v1-proteins.fa

#
grep \> ${genotype}-v1-proteins.fa | sed s/\>// > old_proteins
grep \> ${genotype}-v1.2-proteins.fa | sed s/\>// > new_proteins 

cp ${genotype}-v1.2-proteins.fa new.fa
cat new_proteins | while read line
do
	new=$(grep ${line} old_proteins)
	if [[ ! -z ${new} ]]
	then
		sed -i s/"${line}"/"${new}"/ new.fa 
	else
		new="${line} AED:0.93 eAED:1.00 QI:200|-1|0|1|-1|1|1|0|97"
		sed -i s/"${line}"/"${new}"/ new.fa 
	fi
done


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

blast="../../../../comparative/diamond_blastp/Athaliana_Athaliana/Mguttatus_S1-Athaliana_Athaliana_orthogroup_filtered.m8"

cat new_proteins | while read line
do
	#Handle the Arabidopsis BLAST
	AtID=$(grep ${line} ${blast} | sort -r -k12 | head -1 | cut -f2)
	if [[ ! -z ${AtID} ]]
	then
		#Get the Arabidopsis Description
		AtDesc="Arabidopsis BLAST: $(grep ${AtID} TAIR10_short_functional_descriptions.txt)"
		#Get Arabidopsis GO terms
		AtGO=$(grep ${AtID} gene_association.tair | cut -f5 | tr '|' '\n' | sort | uniq | tr '\n' '|' | sed s/\|$//)
	else
		AtID=NA
		AtDesc=NA
		AtGO=NA
	fi
	#Get the Pfam domains
	grep ${line} ../v1/*.iprscan > tmp
	#If tmp is not empty
	if [ -s tmp ]
	then
		#List the PfamIDs
		PfamID=$(cut -f5 tmp | sort | uniq | tr '\n' '|' | sed s/\|$//)
		#Get the Pfam Descriptions
		PfamDesc="PFAM: $(cut -f6 | tr '\n' ',')"
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
	tr ' ' '\t' | tr ';' ' ' >> functional-annotations.tsv
	#Remove tmp file
	rm tmp
done
