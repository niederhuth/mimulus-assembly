#Read in the pangenome database
pgdb <- fread("results/L1_pangenomeDB.txt.gz",na.strings = c("", "NA"))
pgff <- fread("results/gffWithOgs.txt.gz",na.strings = c("", "NA"))

#Add the arrayID
pge <- setorder(merge(pgdb,pgff[,c(2,15)],by="id"), pgChr, pgOrd, na.last = T)

arraycount <- as.data.frame(table(pge$arrayID))
tandem <- setorder(pge[pge$arrayID %in% arraycount[arraycount$Freq > 1,]$Var1,], arrayID)

#Create two new tables for each of the genomes, removing genes not in the pangenome
L1 <- subset(pgdb, genome=="L1" & !is.na(pgChr))
S1 <- subset(pgdb, genome=="S1" & !is.na(pgChr))
#Combine the files based on the pangenome ID
df1 <- full_join(L1,S1,by=("pgID"))

#If isDirectSyn & isArrayRep is TRUE for both paired genes in both genomes, it is a syntenic pair
synt <- subset(df1, isArrayRep.x & isDirectSyn.x & isArrayRep.y & isDirectSyn.y)
#Get non-syntenic gene pairs
nonsynt <- setdiff(df1, synt)

#Subset out genes that are nonsyntenic orthologs to those in the syntenic pairs
L1_NSog <- subset(nonsynt, id.y %in% synt$id.y)
S1_NSog <- subset(nonsynt, id.x %in% synt$id.x)

#
tmp <- setdiff(setdiff(nonsynt,L1_NSog),S1_NSog)




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

L1 <- pgdb[pgdb$genome=="L1",]
S1 <- pgdb[pgdb$genome=="S1",]

library(dplyr)

tmp <- full_join(L1,S1,,by=("pgID"))
tmp <- subset(ppgdb, isArrayRep & isDirectSyn)

#
tmp <- subset(ppgdb, isArrayRep & isDirectSyn)






test

pc$class <- NA

if()

	pc[pc$L1ogcount==1 & pc$S1ogcount==1 & pc$genecount==1 & pc$isDirectSyn,]






L1 <- pc[pc$L1ogcount==1 & pc$S1ogcount==1 & pdgb$genecount==1,]$genome
S1 <- pc[pc$L1ogcount==1 & pc$S1ogcount==1,]$genome