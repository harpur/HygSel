###
# LAMP Analyses
###


#follows LAMP.sh, takes output, analyzes it relative to anal.r ranges



#!/usr/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
chrs = args[1] #16

#load data sets and packages -------------------------------
source("/media/data1/forty3/brock/GenomeR/VarFunct.r")
rang = read.table(file="/media/data1/forty3/drone/working/hapFLK_RANGES", header=F)
rang5 = rang[rang$V1 == chrs,]

#compile LAMP outputs -----------------
lamp = read.table(paste("chr",chrs,".out",sep=""),header=F)
lamp = as.matrix(lamp, nr = 82)
#lamp[lamp==0] = "yellow"
#lamp[lamp==1] = "black"
#lamp[lamp==2] = "red"
pos = read.table(paste("pos_",chrs,sep=""),header=F)


#extract SNPs within range ---------
	#lamp.rang = c()

blah=outer(as.numeric(rang5$V3), as.numeric(as.character(pos$V1)), "<=") 
blah1=outer(as.numeric(rang5$V4), as.numeric(as.character(pos$V1)), ">=") 
blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) 
	#temp = rang5[blah[,1],]
	#temp = cbind(temp, pos[blah[,2],])
	#lamp.rang = rbind(temp,lamp.rang)
	#names(lamp.rang) = c("chr", "grp", "start", "end", "pos")

#blah[,2] now contains the indices of all columns in lamp that are within my selected regions.
	#regions = lamp[,blah[,2]]
sel.SNP = rep(0, nrow(pos))
sel.SNP[blah[,2]] = "1"


pop = 0 #0,1,2 = C,M,A
sel.C = apply(lamp[c(1:58),], 2, function(x) length(x[x==pop])/length(x) )
con.C = apply(lamp[c(59:82),], 2, function(x) length(x[x==pop])/length(x) )

pop = 1 #0,1,2 = C,M,A
sel.M = apply(lamp[c(1:58),], 2, function(x) length(x[x==pop])/length(x) )
con.M = apply(lamp[c(59:82),], 2, function(x) length(x[x==pop])/length(x) )

pop = 2 #0,1,2 = C,M,A
sel.A = apply(lamp[c(1:58),], 2, function(x) length(x[x==pop])/length(x) )
con.A = apply(lamp[c(59:82),], 2, function(x) length(x[x==pop])/length(x) )

#sel.C.r = apply(regions[c(1:58),], 2, function(x) length(x[x==pop])/length(x) )
#con.C.r = apply(regions[c(59:82),], 2, function(x) length(x[x==pop])/length(x) )

pos$chr = rep(chrs, nrow(pos))
pos$sel.C = sel.C
pos$con.C = con.C
pos$sel.M = sel.M
pos$con.M = con.M
pos$sel.A = sel.A
pos$con.A = con.A
pos$sel.SNP = sel.SNP


write.list(pos)
#boxplot(sel.C,sel.C.r, con.C,con.C.r,notch=T)
#plot(pos$V1, sel.C-con.C,pch=19)

names(lamp) = c("grp","pos","chr", "sel.C", "con.C", "sel.M", "con.M", "sel.A", "con.A", "sel.SNP")











