#!/usr/bin/Rscript
###
# Moving AVG from PLINK FST


#Requires 1 PLINK-formated FST file (e.g. *.weir.fst)
#Outputs all SNPs within high FST windows





###
#Arguments
args = commandArgs(trailingOnly = TRUE)
FSTFILE = args[1]
winSize = args[2] #I've been using 10000

###
#Functions
source("/media/data1/forty3/brock/scripts/GOgetter.r")
source("/media/data1/forty3/brock/scripts/VarFunct.r")
source("/media/data1/forty3/brock/scripts/movingavg.r")
print("functions loaded")




###
# Read in FST file, merge with chromosome positions 
fst23 = read.table(file=FSTFILE,header=T,colClasses=rep("character",3));names(fst23)[3]="FST"
fst23$SNP = paste(fst23$CHROM,fst23$POS, sep="_")
fst23$CHROM=paste("Group",fst23$CHROM,sep="")
write.list(fst23[c(4,1,2)],file="/media/data1/forty3/drone/ScaffConvert/AHBFST")

#Convert scaffold positions to CHR positions with scaffold_to_chr_BAH.pl
system("perl /media/data1/forty3/drone/ScaffConvert/scaffold_to_chr_BAH.pl /media/data1/forty3/drone/ScaffConvert/scaffolds_on_chr.txt /media/data1/forty3/drone/ScaffConvert/AHBFST")
posits = read.table(file="/media/data1/forty3/drone/ScaffConvert/AHBFST_onChr.txt",header=F)
#aggregate(posits$V3, by=list(posits$V2),range) #data checks
fst.raw=cbind(posits,fst23)
fst23 = fst.raw[c(2,3,6)];names(fst23)=c("CHROM","POS","FST")
fst23$CHROM=gsub("chr","",fst23$CHROM)
fst23 = fst23[as.numeric(fst23$FST)>0,]
fst23=fst23[with(fst23, order(as.numeric(POS))), ]
print("Fst File loaded")

###
# Moving average window
Fst=with(fst23, aggregate(as.numeric(FST),by=list(CHROM),function(x) movingAverage(x,n=as.numeric(winSize),centered=F)))
FstPositions=with(fst23, aggregate(as.numeric(POS),by=list(CHROM),function(x) movingAverage(x,n=as.numeric(winSize),centered=F)))
Fst$x=apply(Fst, 1, function(x) as.numeric(unlist(x[2]))[-c(1:1000)])
FstPositions$x=apply(FstPositions, 1, function(x) as.numeric(unlist(x[2]))[-c(1:1000)])
Fst$Pos=FstPositions$x;rm(FstPositions)
Fst=Fst[with(Fst, order(as.numeric(Group.1))), ]
print("Windows complete")



###
#Save
save.name = gsub(".weir.fst", "", FSTFILE)
save.image(file=paste(save.name, ".RDATA", sep="") )





