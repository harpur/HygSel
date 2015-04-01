 ###
 #Fst- TD and Pi-based approaches for Drone Selection: 

 
 #BASH Pre-amble
 #########################################################################
 
cd /media/data1/forty3/drone/vcf_drone

#plink --file DroneSelection --maf 0.1 --recode --make-bed --noweb --mind 0.2 --out CleandDrone 
 
#pop1 contains all MAS selected bees. pop2, control, pop3 FAS
#Run VCFtools to estimate point-estimates of FST
vcftools --maf 0.05 --vcf DroneSelectionFinal.recode.vcf --weir-fst-pop pop1.txt --weir-fst-pop pop2.txt --out pop1_vs_pop2 
vcftools --maf 0.05 --vcf DroneSelectionFinal.recode.vcf --weir-fst-pop pop3.txt --weir-fst-pop pop2.txt --out pop3_vs_pop2 
vcftools --maf 0.05 --vcf DroneSelectionFinal.recode.vcf --weir-fst-pop pop3.txt --weir-fst-pop pop1.txt --out pop1_vs_pop3 
 
 
 
vcftools --maf 0.05 --vcf DroneSelection.vcf --weir-fst-pop highHB.txt --weir-fst-pop lowHB.txt --out high_vs_low --max-missing 0.25 




#One more, all selected VS control (2 vs 1+3=Selpop.txt)
vcftools --maf 0.05 --vcf DroneSelectionFinal.recode.vcf --weir-fst-pop pop2.txt --weir-fst-pop Selpop.txt --out pop2_vs_sel 


 
 
 
 
 #1Kb window
 vcftools --maf 0.1 --vcf DroneSelection.vcf --weir-fst-pop pop3.txt  --fst-window-size 1000 --weir-fst-pop pop2.txt --out pop3_vs_pop2 --max-missing 0.25
 sed '/-/d' pop3_vs_pop2.windowed.weir.fst > pop3_vs_pop2.window.fst

 #########################################################################
 
 
 
 
 
 
 
 
 
 
#Functions
MultiWayOverlapper = function(win.start,win.end,gene.start,gene.end,gene.list) {
  #this is a monster, but basically, looks within each row for genes  overlapping with whatever you want
  blah=outer(as.numeric(unlist(win.start)), as.numeric(unlist(gene.end)), "<=") 
  blah1=outer(as.numeric(unlist(win.end)), as.numeric(unlist(gene.start)), ">=") 
  blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
  if(!is.null(nrow(blah))){return(as.character(unlist(gene.list)[blah[,2]]))}	
  }

lenunique<-function(x){length(unique(x))}



#########################################################################
#Fst in Pop2 vs 3 (control vs FAS)
fst23=read.table(file="/media/data1/forty3/drone/FST/pop3_vs_pop2.weir.fst",header=T, colClasses=c("character","numeric","numeric"))
names(fst23)[3]="FST"
fst23=fst23[which(fst23$FST>=0 & !is.na(fst23$FST)),] 
pdf(file="/media/data1/forty3/drone/FST/HistogramofFASvsControlfst23.pdf")
hist(fst23$FST,col="lightgreen",main=NULL,xlab="fst23")
dev.off()
fst23$SNP=paste(fst23$CHROM,fst23$POS, sep=":")
fst23=fst23[,-c(1,2)]





map=read.table(file="DroneSelection.map");names(map)=c("CHROM", "SNP","NAL","POS");map$NAL=NULL
fst23=merge(fst23,map,by="SNP")
fst23=fst23[order(fst23$CHROM, fst23$POS), ]


source("/media/data1/forty3/brock/scripts/movingavg.r")
Fst=with(fst23, aggregate(FST,by=list(CHROM),function(x) movingAverage(x,n=1000,centered=F)))
FstPositions=with(fst23, aggregate(POS,by=list(CHROM),function(x) movingAverage(x,n=1000,centered=F)))
Fst$x=apply(Fst, 1, function(x) as.numeric(unlist(x[2]))[-c(1:100)])
FstPositions$x=apply(FstPositions, 1, function(x) as.numeric(unlist(x[2]))[-c(1:100)])
Fst$Pos=FstPositions$x;rm(FstPositions)
Fst=Fst[with(Fst, order(as.numeric(Group.1))), ]



#PLOTTING
x11();plot(as.numeric(unlist(Fst$Pos[Fst$Group.1=="16"])),as.numeric(unlist(Fst$x[Fst$Group.1=="16"])), pch=19)
#abline(v=c(HighFstPositions[2,1]),col="red")
#x11();plot(as.numeric(unlist(Fst$x[Fst$Group.1=="1"])), pch=19)
highFstWindow=quantile(unlist(Fst$x),0.98)
HighFstPositions=as.matrix(apply(Fst, 1, function(x) unlist(as.numeric(unlist(x[3]))[which((as.numeric(unlist(x[2]))>=highFstWindow))])))
		#Fst$Group.1 for order of Chromosomes.

##
#To get gene soverlapping these regions, I'll load the GFF file.
#Genes in each region of pi.high and TD.high will overlap with values in the "starts" and "ends"
gff=read.table(file="/media/data1/forty3/brock/scripts/NCBIGFF.txt",col.names=c("chrom","GB","st","en","type"))




genes=gff[which(gff$type=="gene"),]
sts=aggregate(genes$st,by=list(genes$chrom),function(x) x)
ends=aggregate(genes$en,by=list(genes$chrom),function(x) x)
scaff.genes=aggregate(as.character(genes$GB),by=list(genes$chrom),function(x) as.character(x))
#Merge everything together, and then find overlap of metrics (TD and PI) and genes.
#output,for each scaffold, the indices to get out genes.





Find.High.Regions=cbind(gene=scaff.genes,gstart=sts,gend=ends,HighFstPositions)
	#function(win.start,win.end,gene.start,gene.end,gene.list)
Fst.genes=apply(Find.High.Regions, 1, function(x) MultiWayOverlapper(x[7],x[7],x[4],x[6],x[2]) )
Fst.genes=lapply(Fst.genes,unique)


GenicSNPs=genes[genes$GB %in% unlist(Fst.genes),]
#write.table(GenicSNPs[,c(1,3,4,2)],file="GenicSNPsHighFst",col.names=F,row.names=F,quote=F)



snps=c()
for(i in 1:nrow(GenicSNPs)){
	test=cbind(rep(GenicSNPs[i,1],length(c(as.numeric(GenicSNPs[i,3]):as.numeric(GenicSNPs[i,4]))) ),c(as.numeric(GenicSNPs[i,3]):as.numeric(GenicSNPs[i,4])))
	print(i)
	snps=rbind(test,snps)
} 
 
write.table(snps, file="GenicFstP3.snp", col.names=F,row.names=F,quote=F)




#I'd like to get out the highest FST SNPs from these genes.

chroms=intersect(fst23$CHROM, GenicSNPs$chrom)
GenicFST=c()
for(i in chroms){
	#i=3
	fst=fst23[fst23$CHROM==i,];genic=GenicSNPs[GenicSNPs$chrom==i,]
	blah1=outer(unlist(fst$POS),unlist(genic$st), ">=") 
	blah=outer(unlist(fst$POS), unlist(genic$en), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #row is df1, col is df2
	overlap=fst[(blah[,1]),];overlap$GB=genic[blah[,2],2]
	print(i)
	GenicFST=rbind(GenicFST,overlap)	
	} 
rm(fst,genic,blah,blah1,i)



HighestGenicFST=with(GenicFST, aggregate(FST, by=list(GB), max))
names(HighestGenicFST)[1]="GB"
test=merge(HighestGenicFST, GenicFST, by="GB")
test=test[test$x==test$FST,]
test$pos2=gsub(".*_","", test$SNP)
test$pos1=gsub("_.*","", test$SNP)
HighestGenicFST=test
	#188 genes
pdf(file="/media/data1/forty3/drone/R/SelvsConFstPlot.pdf"); hist(fst23$FST,col="lightgreen",xlab="Fst",main=NULL);dev.off()

	
#took highest FST SNPs (>0.99) in that region and looked for associations
test=test[which(test$FST>quantile(fst23$FST,0.98)),]
table(test$CHROM)
write.table(test[,c(7)], file="GenicHighestFSTSNPs.snp", col.names=F,row.names=F, quote=F)


#Randomly chose one of the highest SNPs (if there are many)
x=test[unlist(tapply(1:nrow(test),test$GB,function(x) {if (length(x)>1) {sample(x,min(length(x),1))}else {x}})),]
x=x[complete.cases(x),]
write.table(test[,c(7,8)], file="GenicHighestFSTSNPs.snp", col.names=F,row.names=F, quote=F)
write.table(x[,c(7,8)], file="GenicHighestFSTSNPs.snp", col.names=F,row.names=F, quote=F)






save.image(file="/media/data1/forty3/drone/R/SelectionScansPop3.RData")








#FST FIGURE:
par(mfrow=c(8, 2),mar = rep(2, 4))
for(i in 1:16){
plot(as.numeric(unlist(Fst$Pos[Fst$Group.1==i])),as.numeric(unlist(Fst$x[Fst$Group.1==i])), 
	type="l",
	axes=FALSE,
	xlab=NA,
	ylab=NA,
	ylim=c(0,0.4)
	)
	abline(h=0.17,col="red",lty=2)
}


#FST FIGURE:
pdf(file="FstOutlierExampleChromosome6.pdf")
par(mfrow=c(2, 1),mar = rep(4, 4))
for(i in c(1,6)){
plot(as.numeric(unlist(Fst$Pos[Fst$Group.1==i])),as.numeric(unlist(Fst$x[Fst$Group.1==i])), 
	type="l",
	xlab="Genomic Position",
	ylab="Fst",
	ylim=c(0,0.4)
	)
	abline(h=c(0.06,0.13),col="darkgrey",lty=2)
	abline(h=0.17,col="red",lty=2)
	#abline(v=c(1331658, 1337245),lty=3,col="blue")
}

#quantile(unlist(Fst[1,2]))

#> median(unlist(Fst[1,2]))+1.5*IQR(unlist(Fst[1,2]))
#[1] 0.12626
#> median(unlist(Fst[1,2]))-1.5*IQR(unlist(Fst[1,2]))
#[1] 0.05564069
dev.off()




#Pi and TD don't get calculated when Fst=1 BECAUSE THERE ARE NO SNPS....ugh. 
 #########################################################################
 
 
 
 
 
 
 #########################################################################
#Fst in Pop2 vs 1 (control vs FAS)
fst21=read.table(file="/media/data1/forty3/drone/FST/pop1_vs_pop2.weir.fst",header=T, colClasses=c("character","numeric","numeric"))
names(fst21)[3]="FST"
fst21=fst21[which(fst21$FST>=0 & !is.na(fst21$Fst)),] 
pdf(file="/media/data1/forty3/drone/FST/HistogramofFASvsControlfst21.pdf")
hist(fst21$FST,col="lightgreen",main=NULL,xlab="FST")
dev.off()
fst21$SNP=paste(fst21$CHROM,fst21$POS, sep=":")
fst21=fst21[,-c(1,2)]





map=read.table(file="DroneSelection.map");names(map)=c("CHROM", "SNP","NAL","POS");map$NAL=NULL
fst21=merge(fst21,map,by="SNP")
fst21=fst21[order(fst21$CHROM, fst21$POS), ]


source("/media/data1/forty3/brock/scripts/movingavg.r")
Fst=with(fst21, aggregate(FST,by=list(CHROM),function(x) movingAverage(x,n=1000,centered=F)))
FstPositions=with(fst21, aggregate(POS,by=list(CHROM),function(x) movingAverage(x,n=1000,centered=F)))
Fst$x=apply(Fst, 1, function(x) as.numeric(unlist(x[2]))[-c(1:100)])
FstPositions$x=apply(FstPositions, 1, function(x) as.numeric(unlist(x[2]))[-c(1:100)])
Fst$Pos=FstPositions$x;rm(FstPositions)
Fst=Fst[with(Fst, order(as.numeric(Group.1))), ]



#PLOTTING
#x11();plot(as.numeric(unlist(Fst$Pos[Fst$Group.1=="16"])),as.numeric(unlist(Fst$x[Fst$Group.1=="16"])), pch=19)
#abline(v=c(HighFstPositions[2,1]),col="red")
#x11();plot(as.numeric(unlist(Fst$x[Fst$Group.1=="1"])), pch=19)
highFstWindow=quantile(unlist(Fst$x),0.98)
HighFstPositions=as.matrix(apply(Fst, 1, function(x) unlist(as.numeric(unlist(x[3]))[which((as.numeric(unlist(x[2]))>=highFstWindow))])))
		#Fst$Group.1 for order of Chromosomes.

##
#To get gene soverlapping these regions, I'll load the GFF file.
#Genes in each region of pi.high and TD.high will overlap with values in the "starts" and "ends"
gff=read.table(file="/media/data1/forty3/brock/scripts/NCBIGFF.txt",col.names=c("chrom","GB","st","en","type"))




genes=gff[which(gff$type=="gene"),]
sts=aggregate(genes$st,by=list(genes$chrom),function(x) x)
ends=aggregate(genes$en,by=list(genes$chrom),function(x) x)
scaff.genes=aggregate(as.character(genes$GB),by=list(genes$chrom),function(x) as.character(x))
#Merge everything together, and then find overlap of metrics (TD and PI) and genes.
#output,for each scaffold, the indices to get out genes.





Find.High.Regions=cbind(gene=scaff.genes,gstart=sts,gend=ends,HighFstPositions)
	#function(win.start,win.end,gene.start,gene.end,gene.list)
Fst.genes=apply(Find.High.Regions, 1, function(x) MultiWayOverlapper(x[7],x[7],x[4],x[6],x[2]) )
Fst.genes=lapply(Fst.genes,unique)


GenicSNPs=genes[genes$GB %in% unlist(Fst.genes),]
#write.table(GenicSNPs[,c(1,3,4,2)],file="GenicSNPsHighFst",col.names=F,row.names=F,quote=F)



#snps=c()
#or(i in 1:nrow(GenicSNPs)){
#	test=cbind(rep(GenicSNPs[i,1],length(c(as.numeric(GenicSNPs[i,3]):as.numeric(GenicSNPs[i,4]))) ),c(as.numeric(GenicSNPs[i,3]):as.numeric(GenicSNPs[i,4])))
#	print(i)
#	snps=rbind(test,snps)
#} 
#write.table(snps, file="GenicFstP3.snp", col.names=F,row.names=F,quote=F)




#I'd like to get out the highest FST SNPs from these genes.

chroms=intersect(fst21$CHROM, GenicSNPs$chrom)
GenicFST=c()
for(i in chroms){
	#i=3
	fst=fst21[fst21$CHROM==i,];genic=GenicSNPs[GenicSNPs$chrom==i,]
	blah1=outer(unlist(fst$POS),unlist(genic$st), ">=") 
	blah=outer(unlist(fst$POS), unlist(genic$en), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #row is df1, col is df2
	overlap=fst[(blah[,1]),];overlap$GB=genic[blah[,2],2]
	print(i)
	GenicFST=rbind(GenicFST,overlap)	
	} 
rm(fst,genic,blah,blah1,i)



HighestGenicFST=with(GenicFST, aggregate(FST, by=list(GB), max))
names(HighestGenicFST)[1]="GB"
test=merge(HighestGenicFST, GenicFST, by="GB")
test=test[test$x==test$FST,]
test$pos2=gsub(".*_","", test$SNP)
test$pos1=gsub("_.*","", test$SNP)
HighestGenicFST=test
	#188 genes
pdf(file="/media/data1/forty3/drone/R/SelvsConFstPlot21.pdf"); hist(fst21$FST,col="lightgreen",xlab="Fst",main=NULL);dev.off()

	
#took highest FST SNPs (>0.99) in that region and looked for associations
test=test[which(test$FST>quantile(fst21$FST,0.98)),]
table(test$CHROM)
write.table(test[,c(7)], file="GenicHighestFSTSNPs21.snp", col.names=F,row.names=F, quote=F)


#Randomly chose one of the highest SNPs (if there are many)
x=test[unlist(tapply(1:nrow(test),test$GB,function(x) {if (length(x)>1) {sample(x,min(length(x),1))}else {x}})),]
x=x[complete.cases(x),]
write.table(test[,c(7,8)], file="GenicHighestFSTSNPs.snp", col.names=F,row.names=F, quote=F)
write.table(x[,c(7,8)], file="GenicHighestFSTSNPs.snp", col.names=F,row.names=F, quote=F)






save.image(file="/media/data1/forty3/drone/R/SelectionScansPop12.RData")








#FST FIGURE:
par(mfrow=c(8, 2),mar = rep(2, 4))
for(i in 1:16){
plot(as.numeric(unlist(Fst$Pos[Fst$Group.1==i])),as.numeric(unlist(Fst$x[Fst$Group.1==i])), 
	type="l",
	axes=FALSE,
	xlab=NA,
	ylab=NA,
	ylim=c(0,0.4)
	)
	abline(h=0.17,col="red",lty=2)
}


#FST FIGURE:
pdf(file="FstOutlierExampleChromosome6.pdf")
par(mfrow=c(2, 1),mar = rep(4, 4))
for(i in c(1,6)){
plot(as.numeric(unlist(Fst$Pos[Fst$Group.1==i])),as.numeric(unlist(Fst$x[Fst$Group.1==i])), 
	type="l",
	xlab="Genomic Position",
	ylab="Fst",
	ylim=c(0,0.4)
	)
	abline(h=c(0.06,0.13),col="darkgrey",lty=2)
	abline(h=0.17,col="red",lty=2)
}

#quantile(unlist(Fst[1,2]))

#> median(unlist(Fst[1,2]))+1.5*IQR(unlist(Fst[1,2]))
#[1] 0.12626
#> median(unlist(Fst[1,2]))-1.5*IQR(unlist(Fst[1,2]))
#[1] 0.05564069
dev.off()




#Pi and TD don't get calculated when Fst=1 BECAUSE THERE ARE NO SNPS....ugh. 
 #########################################################################
 
 
 
 
 
 
 
 assoc=read.table(file="HLCandidate.assoc.linear.adjusted",header=T)
 assoc=assoc[assoc$FDR_BY<0.05,]
 
 
 
 
 
 
 
 
 
 
 
