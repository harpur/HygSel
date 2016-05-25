###
# Scaff Map to Chr Map



#converts a .map file of scaffolds to a .map file (same name) of chromosomes
#usage: R CMD BATCH 


#!/usr/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
source("/media/data1/forty3/brock/scripts/VarFunct.r")
fil = args[1]
out = args[2]
print(fil)
map = read.table(file=fil,header=F)
map$CHROM = gsub(":.*","",map$V2)
map$CHROM = paste ("Group", map$CHROM,sep="")


write.list(map[c(2,5,4)],file="/media/data1/forty3/drone/ScaffConvert/CandidateMap")
system("perl /media/data1/forty3/drone/ScaffConvert/scaffold_to_chr_BAH.pl /media/data1/forty3/drone/ScaffConvert/scaffolds_on_chr.txt /media/data1/forty3/drone/ScaffConvert/CandidateMap")

#map = read.table(file=fil,header=T)
conv = read.table(file="/media/data1/forty3/drone/ScaffConvert/CandidateMap_onChr.txt",header=F)
conv$V2 = gsub("chr","",conv$V2)
conv = cbind(map,conv)
write.list(conv[c(7,2,3,8)],file=paste(out, ".chrom.map",sep=""))

#LocusID Scaffold Pos1 Post2*optional