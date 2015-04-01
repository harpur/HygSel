######
# Drone TD and Pi.
#####



#pop1 contains all FAS selected bees. pop2, control, pop3 MAS
#Pi
	vcftools --vcf DroneSelection.vcf  --keep pop1.txt --window-pi 1000 --out Dronepop1 --max-missing 0.25
	vcftools --vcf DroneSelection.vcf  --keep pop2.txt --window-pi 1000 --out Dronepop2 --max-missing 0.25
	vcftools --vcf DroneSelection.vcf  --keep pop3.txt --window-pi 1000 --out Dronepop3 --max-missing 0.25

#TD
	vcftools --vcf DroneSelection.vcf  --keep pop1.txt --TajimaD 1000 --out Dronepop1 --max-missing 0.25
	vcftools --vcf DroneSelection.vcf   --keep pop2.txt --TajimaD 1000 --out Dronepop2 --max-missing 0.25
	vcftools --vcf DroneSelection.vcf   --keep pop3.txt --TajimaD 1000 --out Dronepop3 --max-missing 0.25





#pop1 contains all MAS selected bees. pop2, control, pop3 FAS
R
fileTD="Dronepop3.Tajima.D"

#Identify  low TD and Pi windows from VCFtools outputs (.pi and .TajimaD)
#low defined as <-2 SD
td=read.table(file=fileTD,header=T,colClass=c("character",rep("numeric",3)))
#td=td[-grep("17[.]",td$CHROM),]
#td=td[-grep("18[.]",td$CHROM),]
#td=td[which(td$N_SNPS>15),]
td$norm=(td$TajimaD-mean(td$TajimaD))/sd(td$TajimaD)
td$BIN_END=td$BIN_START+999
#chroms=aggregate(td$CHROM,by=list(td$CHROM),length);chroms=chroms$Group.1[chroms$x>100]
#td=td[td$CHROM %in% chroms,]
#td$st=paste(td$CHROM,"_",td$BIN_START,sep="")
tdp3=td

filePi="Dronepop3.windowed.pi"
pi=read.table(file=filePi,header=T,colClass=c("character",rep("numeric",3)))
#x=grep("17[.]",pi$CHROM)
#pi=pi[-x,]
#x=grep("18[.]",pi$CHROM)
#pi=pi[-x,]
#pi=pi[which(pi$N_VARIANTS>15),]
#chroms=aggregate(pi$CHROM,by=list(pi$CHROM),length);chroms=chroms$Group.1[chroms$x>100]
#pi=pi[pi$CHROM %in% chroms,]
pi$norm=(as.numeric(pi$PI)-mean(as.numeric(pi$PI)))/sd(as.numeric(pi$PI))
pip3=pi





#tdpX and pipX store the population-specific extimates of Pi and tajima's D
#Recall: pop1 contains all MAS selected bees. pop2, control, pop3 FAS
	#I want SNPs selected in FAS (POP3!!!!!)





tdcand=tdp3[which(tdp3$norm<(quantile(tdp3$norm,0.1))),];tdcand$st=NULL
picand=pip3[which(pip3$norm<(quantile(pip3$norm,0.1))),]




chroms=intersect(picand$CHROM, tdcand$CHROM)
overlapping=c()
for(i in chroms){
	blah1=outer(unlist(tdcand$BIN_START[tdcand$CHROM==i]),unlist(picand$BIN_END[picand$CHROM==i]), "<=") 
	blah=outer(unlist(tdcand$BIN_END[tdcand$CHROM==i]), unlist(picand$BIN_END[picand$CHROM==i]), ">=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #row is df1, col is df2
	overlap=tdcand[(blah[,1]),];overlap$Pi=picand[blah[,2],5]
	print(i)
	overlapping=rbind(overlapping,overlap)	
	} 
##Above is incorrect

snps=c()
for(i in 1:nrow(overlapping)){
	test=cbind(rep(overlapping[i,1],1000),c(as.numeric(overlapping[i,2]):as.numeric(overlapping[i,6])))
	print(i)
	snps=rbind(test,snps)
} 
 
write.table(snps, file="LowPiTDCandidatesP3.snp", col.names=F,row.names=F,quote=F)
#Did the same thing for MAS selected.
#Using these SNPs as candidates first. If necessary, I'll come back and re-thin the set. 
	#There are too many, I'll filter the list with FST Data too.
	
	


























#tdpX and pipX store the population-specific extimates of Pi and tajima's D
#Recall: pop1 contains all MAS selected bees. pop2, control, pop3 FAS
#What I care about are windows of low diversity in p1,p3 relative to p2
	#p1 relative to p2
		#note, some windows not calcuble in some pops....I'll continue with this in mind
	tdp1$win=paste(tdp1$CHROM, tdp1$BIN_START, sep="_");tdp2$win=paste(tdp2$CHROM, tdp2$BIN_START, sep="_");
	test=merge(tdp1,tdp2,by="win")
	#95/108 were the same
	with(test,plot(norm.x,norm.y)) #x is p1, y is p2, so low in p1 but not low in p2
	#ok, so none of them are "outliers" in the p2 
	#relative to other windows???


#p3 relative to p2
		#note, some windows not calcuble in some pops....I'll continue with this in mind
	tdp3$win=paste(tdp3$CHROM, tdp3$BIN_START, sep="_");tdp2$win=paste(tdp2$CHROM, tdp2$BIN_START, sep="_");
	test=merge(tdp3,tdp2,by="win")
	#95/108 were the same
	with(test,plot(norm.x,norm.y)) #x is p1, y is p2, so low in p1 but not low in p2
	#ok, so none of them are "outliers" in the p2 
	#relative to other windows???


#p3 relative to p1
		#note, some windows not calcuble in some pops....I'll continue with this in mind
	tdp3$win=paste(tdp3$CHROM, tdp3$BIN_START, sep="_");tdp1$win=paste(tdp1$CHROM, tdp1$BIN_START, sep="_");
	test=merge(tdp3,tdp1,by="win")
	#95/108 were the same
	with(test,plot(norm.x,norm.y)) #x is p1, y is p2, so low in p1 but not low in p2
	#ok, so none of them are "outliers" in the p2 
	#relative to other windows???

	
#tdpX and tdpX store the population-specific extimates of Pi and tajima's D
#Recall: pop1 contains all MAS selected bees. pop2, control, pop3 FAS
#What I care about are windows of low diversity in p1,p3 relative to p2
	#p1 relative to p2
		#note, some windows not calcuble in some pops....I'll continue with this in mind
	tdp1$win=paste(tdp1$CHROM, tdp1$BIN_START, sep="_");tdp2$win=paste(tdp2$CHROM, tdp2$BIN_START, sep="_");
	test=merge(tdp1,tdp2,by="win")
	#95/108 were the same
	with(test,plot(norm.x,norm.y)) #x is p1, y is p2, so low in p1 but not low in p2
	#ok, so none of them are "outliers" in the p2 
	#relative to other windows???


#p3 relative to p2
		#note, some windows not calcuble in some pops....I'll continue with this in mind
	tdp3$win=paste(tdp3$CHROM, tdp3$BIN_START, sep="_");tdp2$win=paste(tdp2$CHROM, tdp2$BIN_START, sep="_");
	test=merge(tdp3,tdp2,by="win")
	#95/108 were the same
	with(test,plot(norm.x,norm.y)) #x is p1, y is p2, so low in p1 but not low in p2
	#ok, so none of them are "outliers" in the p2 
	#relative to other windows???


#p3 relative to p1
		#note, some windows not calcuble in some pops....I'll continue with this in mind
	tdp3$win=paste(tdp3$CHROM, tdp3$BIN_START, sep="_");tdp1$win=paste(tdp1$CHROM, tdp1$BIN_START, sep="_");
	test=merge(tdp3,tdp1,by="win")
	#95/108 were the same
	with(test,plot(norm.x,norm.y)) #x is p1, y is p2, so low in p1 but not low in p2
	#ok, so none of them are "outliers" in the p2 
	#relative to other windows???

	
	






	
##
	
	
as=read.table(file="plink.assoc.linear",header=T)
as=as[which(as$TEST=="ADD" & !is.na(as$P)),]
	
as$CHROM=gsub(":.*","",as$SNP)
	
#use tdp1
chroms=intersect(as$CHROM, tdp1$CHROM)
overlapped=c()
for(i in chroms){
	df1=as[which(as$CHROM==i),]
	df2=tdp1[which(tdp1$CHROM==i),]
	blah1=outer(unlist(df1$BP),unlist(df2$BIN_END), "<=") 
	blah=outer(unlist(df1$BP),unlist(df2$BIN_START), ">=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #row is df1, col is df2
	overlap=df1[(blah[,1]),];overlap$TD=df2[blah[,2],4]
	overlapped=rbind(overlap,overlapped)
	print(i)
}













	
	
	

	
	