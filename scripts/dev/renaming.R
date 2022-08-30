#
df <- read.table("combined_order")

colnames(df) <- c("chr","source","start","stop","gene")

df$up_pos <- NA
for(i in 1:nrow(df)){
	if(df[i,]$source == "maker"){
		maker_pos <- df[i,]$gene
	} else if(i == 1 & df[i,]$source == "Liftoff"){

		df[i,]$up_pos <- maker_pos <- paste("Mg",gsub("Mg","",gsub("_.*","",head(df[df$source=="maker",]$gene,1))),
			"_",gsub("01","1",gsub("chr","0",df[i,]$chr)),"g","00000",sep="")
	} else if(df[i,]$chr != df[i-1,]$chr & df[i,]$source == "Liftoff"){
		df[i,]$up_pos <- maker_pos <- paste("Mg",gsub("Mg","",gsub("_.*","",head(df[df$source=="maker",]$gene,1))),
			"_",gsub("01","1",gsub("chr","0",df[i,]$chr)),"g","00000",sep="")
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
	if(nrow(tmp1) == 1){
		if(tmp1[1,]$overlap_up == TRUE & tmp1[1,]$overlap_down == TRUE & tmp1[1,]$overlap_up_gene != tmp1[1,]$overlap_down_gene){
			df2[df2$gene == tmp1[1,]$gene,]$new_pos <- gsub("0$","5",tmp1[1,]$up_pos)
		} else if(tmp1[1,]$overlap_up == TRUE & tmp1[1,]$overlap_down == FALSE) { 
			df2[df2$gene == tmp1[1,]$gene,]$new_pos <- gsub("0$","2",tmp1[1,]$up_pos)
		} else if(tmp1[1,]$overlap_up == FALSE & tmp1[1,]$overlap_down == TRUE) { 
			df2[df2$gene == tmp1[1,]$gene,]$new_pos <- gsub("0$","8",tmp1[1,]$up_pos)
		} else if(tmp1[1,]$overlap_up == FALSE & tmp1[1,]$overlap_down == FALSE) {
			df2[df2$gene == tmp1[1,]$gene,]$new_pos <- gsub("0$","5",tmp1[1,]$up_pos)
		}
	} else if(nrow(tmp1) > 1){
		tmp1$gene_center <- NA
		for(x in 1:nrow(tmp1)){
			tmp1[x,]$gene_center <- mean(c(tmp1[x,]$start,tmp1[x,]$stop))
		}
		region_center <- mean(c(df[as.numeric(row.names(tmp1[1,]))-1,]$stop,df[as.numeric(row.names(tmp1[nrow(tmp1),]))+1,]$start))
		region_range <- region_center - df[as.numeric(row.names(tmp1[1,]))-1,]$stop
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
		up_genes <- tmp1[tmp1$start < tmp1[tmp1$gene == tmp1$center_gene,]$start,]
		down_genes <- tmp1[tmp1$start > tmp1[tmp1$gene == tmp1$center_gene,]$start,]
		range <- abs(tmp1[tmp1$gene == tmp1$center_gene,]$gene_center - region_center) > region_range/4
		if(nrow(up_genes) > 4 & nrow(down_genes) < 4){
			center_count <- nrow(up_genes) + 1
			df2[df2$gene == tmp1[tmp1$gene==tmp1$center_gene,]$gene,]$new_pos <- gsub("0$",center_count,tmp1[tmp1$gene==tmp1$center_gene,]$up_pos)
			df2[df2$gene == tmp1[tmp1$gene==tmp1$center_gene,]$gene,]$center_gene <- tmp1[1,]$center_gene
		} else if(nrow(up_genes) < 4 & nrow(down_genes) > 4) {
			center_count <- 10 - nrow(down_genes) - 1
			df2[df2$gene == tmp1[tmp1$gene==tmp1$center_gene,]$gene,]$new_pos <- gsub("0$",center_count,tmp1[tmp1$gene==tmp1$center_gene,]$up_pos)
			df2[df2$gene == tmp1[tmp1$gene==tmp1$center_gene,]$gene,]$center_gene <- tmp1[1,]$center_gene
		} else {
			if(abs(tmp1[tmp1$gene == tmp1$center_gene,]$gene_center - region_center) > region_range/4){
				if(tmp1[tmp1$gene == tmp1$center_gene,]$gene_center - region_center > 0){
					center_count <- 6
				} else {
					center_count <- 4
				}
			} else {
				center_count <- 5
			}
			df2[df2$gene == tmp1[tmp1$gene==tmp1$center_gene,]$gene,]$new_pos <- gsub("0$",center_count,tmp1[tmp1$gene==tmp1$center_gene,]$up_pos)
			df2[df2$gene == tmp1[x,]$gene,]$center_gene <- tmp1[1,]$center_gene
				df2[df2$gene == tmp1[x,]$gene,]$region_center <- region_center
		}
		if((center_count - 1) > 3 & nrow(up_genes) == 3){
			count=1
		} else {
			count=1
		}
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
	}
}
#
write.table(df2[,c(5,14)],file="rename.map",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)







