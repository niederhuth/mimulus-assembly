
gene_list="L1genes"
species="Mguttatus_L1"
ref_species="Athaliana"

path2=$(pwd)

#set path to the SpeciesID.txt file
if [ -z ${speciesIDs} ]
then
	speciesIDs=$(ls ${path2}/orthofinder/*/WorkingDirectory/SpeciesIDs.txt)
fi
#Set path to N0.tsv orthofinder output
if [ -z ${orthogroups} ]
then
	orthogroups=$(ls ${path2}/orthofinder/*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv)
fi
column=$(expr 4 + $(grep ${ref_species} ${speciesIDs} | cut -d ' ' -f1 | sed s/\://))

cat ${gene_list} | while read line
do
grep ${line} ${orthogroups} > tmp
if [[ $(wc -l tmp) = 0 ]]
then
echo "${line} NA NA NA" | tr ' ' '\t' >> ${species}_${ref_species}_orthologs.tsv
else
hog=$(cut -f1 tmp) 
og=$(cut -f2 tmp)
ref_gene=$(cut -f${column} tmp | sed s/\ //g)
if [[ $(echo ${ref_gene} | wc -c) = 1 ]]
then
ref_gene="NA"
fi
echo "${line} ${hog} ${og} ${ref_gene}" | tr ' ' '\t' >> ${species}_${ref_species}_orthologs.tsv
fi 
rm tmp
done
