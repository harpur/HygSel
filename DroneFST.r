#!/usr/bin/Rscript
###
# Drone FST functions in R


#Requires 1 PLINK-formated FST file (e.g. *.weir.fst)
#Outputs all SNPs within high FST windows





###
#Arguments
args = commandArgs(trailingOnly = TRUE)
FSTFILE = agrs[1]
winSize = agrs[2] #I've been using 10000

###
#Functions
source("/media/data1/forty3/brock/scripts/GOgetter.r")
source("/media/data1/forty3/brock/scripts/VarFunct.r")
source("/media/data1/forty3/brock/scripts/movingavg.r")





###
# Read in FST file, merge with chromosome positions 
fst23 = read.table(file=FSTFILE,header=T,colClasses=rep("character",3));names(fst23)[3]="FST"
fst23$SNP = paste(fst23$CHROM,fst23$POS, sep="_")
fst23$CHROM=paste("Group",fst23$CHROM,sep="")
#write.list(fst23[c(4,1,2)],file="/media/data1/forty3/drone/ScaffConvert/AHBFST")
#run  perl scaffold_to_chr.pl I've removed this because it's already been done.
posits = read.table(file="/media/data1/forty3/drone/ScaffConvert/AHBFST_onChr.txt",header=F)
#aggregate(posits$V3, by=list(posits$V2),range) #data checks
fst.raw=cbind(posits,fst23)
fst23 = fst.raw[c(2,3,6)];names(fst23)=c("CHROM","POS","FST")
fst23$CHROM=gsub("chr","",fst23$CHROM)
fst23 = fst23[as.numeric(fst23$FST)>0,]
fst23=fst23[with(fst23, order(as.numeric(POS))), ]


###
# Moving average window
Fst=with(fst23, aggregate(as.numeric(FST),by=list(CHROM),function(x) movingAverage(x,n=winSize,centered=F)))
FstPositions=with(fst23, aggregate(as.numeric(POS),by=list(CHROM),function(x) movingAverage(x,n=winSize,centered=F)))
Fst$x=apply(Fst, 1, function(x) as.numeric(unlist(x[2]))[-c(1:1000)])
FstPositions$x=apply(FstPositions, 1, function(x) as.numeric(unlist(x[2]))[-c(1:1000)])
Fst$Pos=FstPositions$x;rm(FstPositions)
Fst=Fst[with(Fst, order(as.numeric(Group.1))), ]


###
#Save
save.name = gsub(".weir.fst", "", FSTFILE)
save.image(file=paste(save.name, ".RDATA", sep="") )




highFstWindow98=quantile(unlist(Fst$x),0.98)
highFstWindow95=quantile(unlist(Fst$x),0.95)







load(file="/media/data1/forty3/drone/FST/SelvsCon/SelectedvsControlFST.RData") 




















HighFstPositions98=as.matrix(apply(Fst, 1, function(x) unlist(as.numeric(unlist(x[3]))[which((as.numeric(unlist(x[2]))>=highFstWindow98))])))
#HighFstPositions95=as.matrix(apply(Fst, 1, function(x) unlist(as.numeric(unlist(x[3]))[which((as.numeric(unlist(x[2]))>=highFstWindow95))])))
#HighFstPositions90=as.matrix(apply(Fst, 1, function(x) unlist(as.numeric(unlist(x[3]))[which((as.numeric(unlist(x[2]))>=highFstWindow90))])))


dataset=HighFstPositions98
HighSNPs98=c()
chr = apply(dataset, 1,function(x) length(as.numeric(unlist(x))))
chr = names(chr[chr>0])
for (i in chr){
	test=fst23[fst23$CHROM==i,]
	low=as.numeric(unlist(dataset[i,1]))-1000
	high=as.numeric(unlist(dataset[i,1]))+1000
	blah=outer(as.numeric(unlist(low)), as.numeric(unlist(test$POS)), "<=") 
	blah1=outer(as.numeric(unlist(high)), as.numeric(unlist(test$POS)), ">=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	blah=(test[blah[,2],])
	blah=blah[!duplicated(blah),]
	HighSNPs98=rbind(HighSNPs98,blah)
	print(i)
	#x11();plot(as.numeric(unlist(Fst$Pos[Fst$Group.1==i])),as.numeric(unlist(Fst$x[Fst$Group.1==i])), pch=19)
	#abline(v=blah$PO,col="red")
}



#With these SNPs, create "regions list" for association work:
dataset = HighSNPs98
scaff = gsub(":.*","",dataset$SNP)
Unscaff = unique(scaff) # 19 unique
for (i in Unscaff){
	print(i)
	snps = dataset$SNP[which(scaff==i)]
	num=which(unique(scaff)==i)
	sink(file=paste("Candidates98_", "ALL", ".set",sep=""),append=T)
	cat("\n")
	cat(paste("Set",num,sep=""))
	cat("\n")
	write.table(snps,col.names=F,quote=F,row.names=F)
	cat("END")
	cat("\n")
	sink()
}


#Made 2 files Candiates98.set, and Candiates95.set, formatted as below, for PLINK SET tests

#Set1
#1.32:131761
#1.32:131040
#1.32:131037
#1.32:131008
#...
#END

#Also wrote out candiadates for selective association
write.table(HighSNPs98$SNP,file="CandidateSNPs98.snps",col.names=F,quote=F,row.names=F)


save.image(file="DroneSelection.RData")