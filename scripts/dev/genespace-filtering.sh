

#Genes no BLAST hit
awk -v FS="\t" -v OFS="\t" '$1==""' L1_pangenomeDB.txt

#Possible tandem duplicates?
awk -v FS="\t" -v OFS="\t" '$10=="TRUE" && $12="FALSE" && $7=="L1"' L1_pangenomeDB.txt