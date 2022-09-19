#
df <- read.table("final_valid_gene_order.tsv")

colnames(df) <- c("chr","source","start","stop","gene")

df$up_pos <- NA
new_count=10
for(i in 1:nrow(df)){
	if(df[i,]$source == "maker"){
		maker_pos <- df[i,]$gene
	} else if(i == 1 & df[i,]$source == "Liftoff"){
		if(nrow(head(df[df$chr == df[i,]$chr & df$source == "maker",],1)) == 1){
			pos_num=as.numeric(gsub(".*g","",head(df[df$chr == df[i,]$chr & df$source == "maker",],1)$gene))-10
			pos_num=paste(paste(rep(0,5-nchar(pos_num)),collapse=""),pos_num,sep="")
			df[i,]$up_pos <- maker_pos <- paste("Mg",gsub("Mg","",gsub("_.*","",head(df[df$source=="maker",]$gene,1))),
				"_",gsub("contig_.*","UN",gsub("01","1",gsub("chr","0",df[i,]$chr))),"g",pos_num,sep="")
		} else if(nrow(head(df[df$chr == df[i,]$chr & df$source == "maker",],1)) == 0){
			pos_num=as.numeric(gsub(".*g","",tail(df[df$source=="maker",],1)$gene))+new_count
			pos_num=paste(paste(rep(0,5-nchar(pos_num)),collapse=""),pos_num,sep="")
			df[i,]$up_pos <- maker_pos <- paste("Mg",gsub("Mg","",gsub("_.*","",head(df[df$source=="maker",]$gene,1))),
				"_",gsub("contig_.*","UN",gsub("01","1",gsub("chr","0",df[i,]$chr))),"g",pos_num,sep="")
			new_count=new_count+10
		}
	} else if(df[i,]$chr != df[i-1,]$chr & df[i,]$source == "Liftoff"){
		if(nrow(head(df[df$chr == df[i,]$chr & df$source == "maker",],1)) == 1){
			pos_num=as.numeric(gsub(".*g","",head(df[df$chr == df[i,]$chr & df$source == "maker",],1)$gene))-10
			pos_num=paste(paste(rep(0,5-nchar(pos_num)),collapse=""),pos_num,sep="")
			df[i,]$up_pos <- maker_pos <- paste("Mg",gsub("Mg","",gsub("_.*","",head(df[df$source=="maker",]$gene,1))),
				"_",gsub("contig_.*","UN",gsub("01","1",gsub("chr","0",df[i,]$chr))),"g",pos_num,sep="")
		} else if(nrow(head(df[df$chr == df[i,]$chr & df$source == "maker",],1)) == 0){
			pos_num=as.numeric(gsub(".*g","",tail(df[df$source=="maker",],1)$gene))+new_count
			pos_num=paste(paste(rep(0,5-nchar(pos_num)),collapse=""),pos_num,sep="")
			df[i,]$up_pos <- maker_pos <- paste("Mg",gsub("Mg","",gsub("_.*","",head(df[df$source=="maker",]$gene,1))),
				"_",gsub("contig_.*","UN",gsub("01","1",gsub("chr","0",df[i,]$chr))),"g",pos_num,sep="")
			new_count=new_count+10
		}
	} else if(df[i,]$source == "Liftoff") {
		df[i,]$up_pos <- maker_pos
	}
}
#
df$count <- NA
for(i in 1:nrow(df)){
	if(df[i,]$source == "Liftoff" ){
		df[i,]$count <- nrow(subset(df, chr == df[i,]$chr & up_pos == df[i,]$up_pos))
	}
}
#
df$overlap_down_gene <- df$overlap_down <- df$overlap_up_gene <- df$overlap_up <- NA
for(i in 1:nrow(df)){
	if(i != 1){
		if(df[i,]$start <= df[i-1,]$stop & df[i-1,]$source == "maker" & df[i,]$chr == df[i-1,]$chr){
			df[i,]$overlap_up <- TRUE
			df[i,]$overlap_up_gene <- df[i-1,]$gene
		} else {
			df[i,]$overlap_up <- FALSE
		}
	}
	if(i != nrow(df)){ 
		if(df[i,]$stop <= df[i+1,]$start & df[i+1,]$source == "maker" & df[i,]$chr == df[i+1,]$chr){
			df[i,]$overlap_down <- TRUE
			df[i,]$overlap_down_gene <- df[i+1,]$gene
		} else {
			df[i,]$overlap_down <- FALSE
		}
	}
}
#Rename genes
df2 <- subset(df, source == "Liftoff")
df2$new_pos <- df2$center_gene <- df2$region_center <- NA
df3 <- unique(df2[,c("chr","up_pos")])
for(i in 1:nrow(df3)){
	tmp1 <- subset(df2, chr == df3[i,1] & up_pos == df3[i,2])
	#If only one new gene is added between two existing genes
	if(nrow(tmp1) == 1){
		#if the gene overlaps with both the upstream & downstream gene, it is assigned a 5
		if(tmp1[1,]$overlap_up == TRUE & tmp1[1,]$overlap_down == TRUE & tmp1[1,]$overlap_up_gene != tmp1[1,]$overlap_down_gene){
			df2[df2$gene == tmp1[1,]$gene,]$new_pos <- gsub("0$","5",tmp1[1,]$up_pos)
		#if the gene overlaps with the upstream gene and not the downstream gene, it is assigned a 2
		} else if(tmp1[1,]$overlap_up == TRUE & tmp1[1,]$overlap_down == FALSE) { 
			df2[df2$gene == tmp1[1,]$gene,]$new_pos <- gsub("0$","2",tmp1[1,]$up_pos)
		#if the gene overlaps with the downstream gene and not the upstream gene, it is assigned an 8
		} else if(tmp1[1,]$overlap_up == FALSE & tmp1[1,]$overlap_down == TRUE) { 
			df2[df2$gene == tmp1[1,]$gene,]$new_pos <- gsub("0$","8",tmp1[1,]$up_pos)
		#if the gene doesn't overlap with either, then it is assigned a 5
		} else if(tmp1[1,]$overlap_up == FALSE & tmp1[1,]$overlap_down == FALSE) {
			df2[df2$gene == tmp1[1,]$gene,]$new_pos <- gsub("0$","5",tmp1[1,]$up_pos)
		}
	#If there is more than one gene in the region
	} else if(nrow(tmp1) > 1 & nrow(tmp1) <= 9){
		#Get the center point for each new gene
		tmp1$gene_center <- NA
		for(x in 1:nrow(tmp1)){
			tmp1[x,]$gene_center <- mean(c(tmp1[x,]$start,tmp1[x,]$stop))
		}
		#Calculate the center position & range for the region
		region_center <- mean(c(df[as.numeric(row.names(tmp1[1,]))-1,]$stop,df[as.numeric(row.names(tmp1[nrow(tmp1),]))+1,]$start))
		region_range <- region_center - df[as.numeric(row.names(tmp1[1,]))-1,]$stop
		#Identify the central most gene
		tmp1$center_gene <- NA
		if(nrow(subset(tmp1,start < region_center & stop > region_center)) == 1){
			tmp1$center_gene <- subset(tmp1,start < region_center & stop > region_center)$gene
		} else if(nrow(subset(tmp1,start < region_center & stop > region_center)) >1){
			closest <- which.min(abs(tmp1$gene_center - region_center))
			tmp1$center_gene <- tmp1[x,]$gene
		}
		if(TRUE %in% is.na(tmp1$center_gene)){
			closest <- which.min(abs(tmp1$gene_center - region_center))
			tmp1$center_gene <- tmp1[closest,]$gene
		}
		#Get the number of genes upstream & downstream of the center gene
		up_genes <- tmp1[tmp1$start < tmp1[tmp1$gene == tmp1$center_gene,]$start,]
		down_genes <- tmp1[tmp1$start > tmp1[tmp1$gene == tmp1$center_gene,]$start,]
		range <- abs(tmp1[tmp1$gene == tmp1$center_gene,]$gene_center - region_center) > region_range/4
		#If there are more than 4 genes upstream and less than 4 downstream
		if(nrow(up_genes) > 4 & nrow(down_genes) < 4){
			#Adjust the position of the center gene accordingly
			center_count <- nrow(up_genes) + 1
			df2[df2$gene == tmp1[tmp1$gene==tmp1$center_gene,]$gene,]$new_pos <- gsub("0$",center_count,tmp1[tmp1$gene==tmp1$center_gene,]$up_pos)
			df2[df2$gene == tmp1[tmp1$gene==tmp1$center_gene,]$gene,]$center_gene <- tmp1[1,]$center_gene
		#If there are more than 4 genes downstream and less than 4 upstream
		} else if(nrow(up_genes) < 4 & nrow(down_genes) > 4) {
			#Adjust the position of the center gene accordingly
			center_count <- 10 - nrow(down_genes) - 1
			df2[df2$gene == tmp1[tmp1$gene==tmp1$center_gene,]$gene,]$new_pos <- gsub("0$",center_count,tmp1[tmp1$gene==tmp1$center_gene,]$up_pos)
			df2[df2$gene == tmp1[tmp1$gene==tmp1$center_gene,]$gene,]$center_gene <- tmp1[1,]$center_gene
		} else {
			#Check how close to the center the center gene really is, if ts greater than a quarter then 1/4th the range then we adjust
			if(abs(tmp1[tmp1$gene == tmp1$center_gene,]$gene_center - region_center) > region_range/4){
				#If the difference is greater than zero, then its downstream of the center, so change the count to 6
				if(tmp1[tmp1$gene == tmp1$center_gene,]$gene_center - region_center > 0){
					center_count <- 6
				#Otherwise its upstream, so change the count to 4
				} else {
					center_count <- 4
				}
			#Otherwise it stays as 5
			} else {
				center_count <- 5
			}
			#Set the name for the center gene
			df2[df2$gene == tmp1[tmp1$gene==tmp1$center_gene,]$gene,]$new_pos <- gsub("0$",center_count,tmp1[tmp1$gene==tmp1$center_gene,]$up_pos)
			df2[df2$gene == tmp1[x,]$gene,]$center_gene <- tmp1[1,]$center_gene
			df2[df2$gene == tmp1[x,]$gene,]$region_center <- region_center
		}
		#Now we deal with the other genes
		if((center_count - 1) == nrow(up_genes) || (center_count - 1) == ((floor((center_count - 1)/2)*nrow(up_genes)))){
			count=0
		} else {
			count=1
		}
		#If the number of upstream genes is more than 1
		if(nrow(up_genes) > 1){
			for(x in 1:nrow(up_genes)){
				count=count+floor((center_count - 1)/nrow(up_genes))
				df2[df2$gene == tmp1[x,]$gene,]$new_pos <- gsub("0$",count,tmp1[x,]$up_pos)
				df2[df2$gene == tmp1[x,]$gene,]$center_gene <- tmp1[1,]$center_gene
				df2[df2$gene == tmp1[x,]$gene,]$center_gene <- tmp1[1,]$center_gene
				df2[df2$gene == tmp1[x,]$gene,]$region_center <- region_center
			}
		} else if(nrow(up_genes) == 1){
			count=count+floor((center_count - 1)/2)
			df2[df2$gene == up_genes[1,]$gene,]$new_pos <- gsub("0$",count,up_genes[1,]$up_pos)
			df2[df2$gene == up_genes[1,]$gene,]$center_gene <- up_genes[1,]$center_gene
			df2[df2$gene == up_genes[1,]$gene,]$region_center <- region_center 
		}
		count=center_count
		if(nrow(down_genes) > 1){
			for(x in (nrow(tmp1)-nrow(down_genes)+1):nrow(tmp1)){
				count=count+floor((9 - center_count)/nrow(down_genes))
				df2[df2$gene == tmp1[x,]$gene,]$new_pos <- gsub("0$",count,tmp1[x,]$up_pos)
				df2[df2$gene == tmp1[x,]$gene,]$center_gene <- tmp1[1,]$center_gene
				df2[df2$gene == tmp1[x,]$gene,]$region_center <- region_center
			} 
		} else if(nrow(down_genes) == 1){
			count=count+floor((center_count - 1)/2)
			df2[df2$gene == down_genes[1,]$gene,]$new_pos <- gsub("0$",count,down_genes[1,]$up_pos)
			df2[df2$gene == down_genes[1,]$gene,]$center_gene <- down_genes[1,]$center_gene
			df2[df2$gene == down_genes[1,]$gene,]$region_center <- region_center 
		}
	} else {
		print("Warning greater than 9 genes in region")
		for(x in 1:nrow(tmp1)){
			count=paste("x",x,sep="")
			df2[df2$gene == tmp1[x,]$gene,]$new_pos <- gsub("0$",count,tmp1[x,]$up_pos)
			df2[df2$gene == tmp1[x,]$gene,]$center_gene <- tmp1[1,]$center_gene
			df2[df2$gene == tmp1[x,]$gene,]$region_center <- region_center
		}
	}
}
#
write.table(df2[,c(5,14)],file="rename.map",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)







