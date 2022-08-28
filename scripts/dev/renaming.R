#
S1new <- read.table("../liftoff/S1/combined_order")

colnames(S1new) <- "gene"

S1new$genome <- gsub("_.*","",gsub("Mg","",S1new$gene))
S1new$chr <- gsub("g.*","",gsub("Mg.._","",S1new$gene))
S1new$pos <- gsub("MgL1.*","0",gsub("MgS1.*g","",S1new$gene))

for(i in 1:nrow(S1new)){
	if(S1new[i]$pos==0){
		S1new[i]$pos
	}
}