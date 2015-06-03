#!/usr/bin/Rscript




source("/media/data1/forty3/brock/scripts/GOgetter.r")
source("/media/data1/forty3/brock/scripts/VarFunct.r")
source("/media/data1/forty3/brock/scripts/movingavg.r")




args = commandArgs(trailingOnly = TRUE)
hapsfile = args[1]

#hapsfile="DroneSamps_16.phased.haps.out"

haps=read.table(file=hapsfile)
hapsfile=gsub(".haps.out",".sample",hapsfile)
samp=read.table(file=hapsfile,header=T)
pheno=read.table(file="/media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt");names(pheno)[1]="ID_1";pheno$V2=NULL
samp=samp[-1,];samp=samp[-3]
samp=merge(samp,pheno,by="ID_1");samp$plink_pheno=NULL

#samp + haps = PED file


#samp=rbind(samp,samp);samp=samp[order(samp[,1]),]
haps$zer=rep("0", nrow(haps))
maps=haps[c(1,2,ncol(haps), 3)]
haps=haps[-c(1,2,ncol(haps), 3,4,5)]
haps1=as.matrix(haps)
haps1[haps1=="0"]=23;haps1[haps1=="1"]=2;haps1[haps1=="23"]=1;
thaps=cbind(maps, haps1)
write.list(thaps, file=paste(hapsfile, ".tped",sep=""))
write.list(samp, file=paste(hapsfile, ".tfam",sep=""))
write.list(maps, file=paste(hapsfile, ".map",sep=""))


