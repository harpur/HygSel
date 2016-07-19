#Compare FST among lineages 



#Load Fst and Pi and GFF

cm.fst = read.table(file="M_vs_C.weir.fst",header=T, colClasses = c("character", "numeric",  "numeric")); names(cm.fst)[3] = "cm.fst"
am.fst = read.table(file="M_vs_A.weir.fst",header=T, colClasses = c("character", "numeric",  "numeric")); names(am.fst)[3] = "am.fst"
ca.fst = read.table(file="A_vs_C.weir.fst",header=T, colClasses = c("character", "numeric",  "numeric")); names(ca.fst)[3] = "ca.fst"


# Mung Data frames ------------------------------------------
#create SNP ID
cm.fst$SNP = paste(cm.fst$CHROM, cm.fst$POS, sep="_")
ca.fst$SNP = paste(ca.fst$CHROM, ca.fst$POS, sep="_"); ca.fst = ca.fst[c(3,4)]
am.fst$SNP = paste(am.fst$CHROM, am.fst$POS, sep="_"); am.fst = am.fst[c(3,4)]


#merge
amc.fst = merge(cm.fst, ca.fst, by = "SNP")
amc.fst = merge(amc.fst, am.fst, by = "SNP")
amc.fst = amc.fst[complete.cases(amc.fst),]
amc.fst$CHROM = paste("Group", amc.fst$CHROM, sep="")


# FST ----------------------------------------
high.Fst = c()
chrom = intersect(regions$GRP, amc.fst$CHROM)
for(i in chrom){
	win.temp = regions[regions$GRP==i,]
	deg.temp = amc.fst[which(as.character(amc.fst$CHROM)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$GPOS)), ">=") 
	blah1=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$end)), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	high.Fst = rbind(temp,high.Fst)
	print(i)
}


high.Fst.genic = c()
chrom = intersect(gff$chrom, amc.fst$CHROM)
for(i in chrom){
	win.temp = gff[gff$chrom==i,]
	deg.temp = amc.fst[which(as.character(amc.fst$CHROM)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$start)), ">=") 
	blah1=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$end)), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	high.Fst.genic  = rbind(temp,high.Fst.genic )
	print(i)
}


 
 
boxplot(high.Fst.genic$cm.fst, high.Fst$cm.fst, high.Fst.genic$am.fst, high.Fst$am.fst,high.Fst.genic$ca.fst, high.Fst$ca.fst )
 
 

