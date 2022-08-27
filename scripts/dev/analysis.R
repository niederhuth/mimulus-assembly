#Set some variables
tandem_cutoff=10 #Must be within this number of genes of a gene with same arrayID to be classified as "tandem"

#Load libraries
library(GENESPACE)
library(ggplot2)
library(dplyr)

#Read in the pangenome database
pgdb <- fread("results/L1_pangenomeDB.txt.gz",na.strings = c("", "NA"))
pgff <- fread("results/gffWithOgs.txt.gz",na.strings = c("", "NA"))

#Create two new tables for each of the genomes, removing genes not in the pangenome
L1 <- subset(pgdb, genome=="L1" & !is.na(pgChr))
S1 <- subset(pgdb, genome=="S1" & !is.na(pgChr))

#Create a separate list of unique genes to each genome, not in pangenome
L1u <- subset(pgdb, genome=="L1" & is.na(pgChr))
S1u <- subset(pgdb, genome=="S1" & is.na(pgChr))

#Combine the files based on the pangenome ID
LS <- full_join(L1,S1,by=("pgID"))

#Calculate putative CNV & PAV from pgID
cnv <- dcast(subset(pgdb, !is.na(pgChr)), pgChr + pgOrd + pgID ~ genome, value.var="id", fun.aggregate=length)

#If isDirectSyn & isArrayRep is TRUE for both paired genes in both genomes, it is a syntenic pair
#Subset syntenic genes and write out lists of these genes
synt <- subset(LS, isArrayRep.x & isDirectSyn.x & isArrayRep.y & isDirectSyn.y)
write.table(unique(synt$id.x),file="results/L1_syntelogs.txt",
	quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(unique(synt$id.y),file="results/S1_syntelogs.txt",
	quote=FALSE,col.names=FALSE,row.names=FALSE)
#Get High-confidence 1x1 syntelogs and write out lists of these genes
synt1x1 <- subset(synt, pgID %in% subset(cnv, L1==1 & S1==1)$pgID)
write.table(synt1x1$id.x,file="results/L1_1x1_syntelogs.txt",
	quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(synt1x1$id.y,file="results/S1_1x1_syntelogs.txt",
	quote=FALSE,col.names=FALSE,row.names=FALSE)
#Get sets of non-syntenic genes and write out lists of these genes
nonsynt <- setdiff(LS, synt)
L1nonsynt <- subset(nonsynt, !(id.x %in% synt$id.x))
S1nonsynt <- subset(nonsynt, !(id.y %in% synt$id.y))
write.table(c(unique(L1nonsynt$id.x),L1u$id),file="results/L1_nonsyntenic.txt",
	quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(c(unique(S1nonsynt$id.y),S1u$id),file="results/S1_nonsyntenic.txt",
	quote=FALSE,col.names=FALSE,row.names=FALSE)

#Get genes in arraypairs arrays
#Add the arrayID to pgdb
pge <- setorder(merge(pgdb,pgff[,c(2,15)],by="id"), pgChr, pgOrd, na.last = T)
#Count count arrays with more than one gene in them
arraycount <- as.data.frame(table(unique(pge[,c(1,15)])$arrayID))
#subset out those genes with an array count of more than 1
arraypairs <- setorder(subset(pge, arrayID %in% subset(arraycount, Freq > 1)$Var1), arrayID)
#Create output dataframe
columns=c("genome","arrayID","id","isTandem","geneCount","arrayGenes")
tandem <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(tandem) <- columns
#Now for a complicated loop to identify tandem arrays
for(i in unique(arraypairs$arrayID)){
	#some arraypairs arrays may belong to multiple orthogroups, collapse these so each gene is represented once
	x <- setorder(unique(subset(arraypairs, arrayID==i)[,c(1,8,9,10,15)]), ord)
	genome <- unique(x$genome)
	#Check to make sure all of these are on the same chromosome!
	if(uniqueN(x$chr)==1){
		#Most array pairs are simple arraypairs duplicates and faster to process
		#If the number of genes is 2 and the distance is <= our cutoff classify as arraypairs duplicates
		if(nrow(x)==2 & abs(x[1]$ord - x[2]$ord) <= tandem_cutoff){
			for(gene in x$id){
				tandem <- rbind(tandem,data.frame(genome=genome,arrayID=i,id=gene,isTandem=TRUE,geneCount=2,
								arrayGenes=paste(x$id,collapse=", ")))
			} 
		#If the number of genes is 2 and distance > cutoff, classify as a dispersed
		} else if(nrow(x)==2 & abs(x[1]$ord - x[2]$ord) > tandem_cutoff) {
			tandem <- rbind(tandem,data.frame(genome=genome,arrayID=i,id=gene,isTandem=FALSE,geneCount=NA,
							arrayGenes=NA))
		#If the arraypairs have more than one gene, it requires more work
		} else {
			#Create an empty dataframe to calculate pairwise distances between all genes
			x2 <- data.frame(matrix(ncol=nrow(x),nrow=nrow(x)))
			#Set row & column names for the table
			row.names(x2) <- colnames(x2) <- x$id
			#Calculate absolute distances between all genes
			for(gene in 1:nrow(x)){
				x2[,gene] <- abs(x$ord - x[gene]$ord)
			}
			for(gene in x$id){
				#tandem if a gene's closest array pair is greater than the cutoff
				#If so, classify as FALSE
				if(min(x2[x2[,gene] != 0,gene]) > tandem_cutoff){
					tandem <- rbind(tandem,data.frame(genome=genome,arrayID=i,id=gene,isTandem=FALSE,geneCount=NA,
									arrayGenes=NA))
				#Otherwise, the gene is part of an arraypairs array
				} else {
					#Get a list of all those genes within the cutoff distance of the arraypairs array
					x3 <- row.names(x2[x2[gene,] <= tandem_cutoff,])
					#Now we are going to loop back over those genes in that initial list and search
					#for other genes which may be within range of the arraypairs array, but not the initial gene
					#these will get added to the list
					newgenes <- 1 #set this to initialize the for loop
					#Keep the loop going until no new genes are added
					while(newgenes != 0){
						#Create x4 list same as x3
						x4 <- x3
						#Loop over genes in x3
						for(gene2 in x3){
							#Extend out the array and add the new genes to x4
							x4 <- unique(c(x4,row.names(x2[x2[gene2,] <= tandem_cutoff,])))
						}
						#Set newgenes to the difference between x4 & x3
						#The loop will terminate when no new genes are added
						newgenes <- length(x4) - length(x3)
						#Set x3 to now be equal to x4
						x3 <- x4
					}
					#Now we wll report it
					tandem <- rbind(tandem,data.frame(genome=genome,arrayID=i,id=gene,isTandem=TRUE,
								geneCount=length(x3),arrayGenes=paste(x3,collapse=", ")))
				}
			}
		}	
	#If the genes are not on the same chromosome, they can't be a arraypairs array. Throw a warning.
	} else {
		print(paste("Warning! Array ID",unique(x$arrayID),"genes on different chromosomes. Skipping",sep=" "))
	}
}
#Cleanup
rm(x,x2,x3,x4,gene,gene2,i,columns)
#Output tables of the tandem arrays
write.table(subset(tandem,genome=="L1" & isTandem),file="results/L1_tandem_arrays.tsv",sep="\t",quote=FALSE,row.names=FALSE)
write.table(subset(tandem,genome=="S1" & isTandem),file="results/S1_tandem_arrays.tsv",sep="\t",quote=FALSE,row.names=FALSE)

#








#Subset out genes that are nonsyntenic orthologs to those in the syntenic pairs
L1_NSog <- subset(nonsynt, id.y %in% synt$id.y)
S1_NSog <- subset(nonsynt, id.x %in% synt$id.x)

#
NS <- setdiff(setdiff(nonsynt,L1_NSog),S1_NSog)




pc <- pgdb
#pgff <- fread("results/gffWithOgs.txt.gz",na.strings = c("", "NA"))

#Count occurrence of a gene in pc
count <- as.data.frame(table(pc$id))
pc$genecount <- count$Freq[match(pc$id, count$Var1)]
#Count occurances of pgID in pc
count <- as.data.frame(table(pc$pgID))
pc$pgIDcount <- count$Freq[match(pc$pgID, count$Var1)]
#Count occurrence of an pgID for L1 in pc
count <- as.data.frame(table(pc[pc$genome=="L1",]$pgID))
pc$L1pgIDcount <- count$Freq[match(pc$pgID, count$Var1)]
pc$L1pgIDcount <- ifelse(is.na(pc$L1pgIDcount),0,pc$L1pgIDcount)
#Count occurrence of an pgID for S1 in pc
count <- as.data.frame(table(pc[pc$genome=="S1",]$pgID))
pc$S1pgIDcount <- count$Freq[match(pc$pgID, count$Var1)]
pc$S1pgIDcount <- ifelse(is.na(pc$S1pgIDcount),0,pc$S1pgIDcount)
#Count occurrence of an orthogroup in pc
count <- as.data.frame(table(pc$og))
pc$ogcount <- count$Freq[match(pc$og, count$Var1)]
#Count occurrence of an orthogroup for L1 in pc
count <- as.data.frame(table(pc[pc$genome=="L1",]$og))
pc$L1ogcount <- count$Freq[match(pc$og, count$Var1)]
pc$L1ogcount <- ifelse(is.na(pc$L1ogcount),0,pc$L1ogcount)
#Count occurrence of an orthogroup for S1 in pc
count <- as.data.frame(table(pc[pc$genome=="S1",]$og))
pc$S1ogcount <- count$Freq[match(pc$og, count$Var1)]
pc$S1ogcount <- ifelse(is.na(pc$S1ogcount),0,pc$S1ogcount)
#cleanup a bit
rm(count)
