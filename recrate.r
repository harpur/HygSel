#This is a QD analysis of recombination rate of genes within 
#selected windows. This is to satisfy reviewer 1 that the 
#candidate genes are not in areas of low recombination


#load library -----------------------------------------------------------------
library(hygsel)

#Load data sets ---------------------------------------------------------------
data(gff)
data(associatedGenes)
recom <- read.table("Chrom10kbrec.txt",header <- T) #from Liu et al. (2015; Genome Biology) 
recom$chrom <- gsub("Group","", recom$Group)

#extract genes from GFF -------------------------------------------------------
associated.genes <- gff[which(gff$GB %in% hi.genes$GB),]
associated.genes <- associated.genes[which(associated.genes$type=='gene'),]
associated.genes <- gff[which(gff$type == 'gene'),]
associated.genes <- associated.genes[-(grep("^17[.]",associated.genes$chrom)),]

#Extract recombination rate per gene ------------------------------------------
	#this is a silly loop, but it works.

rec.genes <- c()
chrom <- intersect(recom$chrom, associated.genes$chrom)
for(i in chrom){
	deg.temp <- associated.genes[associated.genes$chrom==i,]
	win.temp <- recom[which(as.character(recom$chrom)==as.character(i)),]
	blah <- outer(as.numeric(deg.temp$start), as.numeric(as.character(win.temp$End)), "<=") 
	blah1 <- outer(as.numeric(deg.temp$end), as.numeric(as.character(win.temp$Start)), ">=") 
	blah <- (which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp <- deg.temp[blah[,1],]
	temp <- cbind(temp, win.temp[blah[,2],])
	rec.genes   <- rbind(temp,rec.genes )
	print(i)
}


#Is recomb rate different for associated vs not? ------------------------------
rec.genes <- rec.genes[which(rec.genes$Recrate>0),]
rec.model <- aov(log10(rec.genes$Recrate)~rec.genes$associated.genes)
summary(rec.model)







