###
# Plot 6.2 
###





source("/media/data1/forty3/drone/git/GenomeR/VCFFunctions.r")





i = "out.recode.vcf" 
vcf = Read.VCF()
samps = read.table(file="/media/data1/forty3/drone/git/DroneSamps.txt",header=F)


pos=sapply(vcf,function(x) return((x[2])))
pos = pos[-1]
allele  = sapply(vcf,function(x) substr(unlist(x)[10:134],1 ,1 )) 
allele = allele[,-1]


#identify major allele in control pop [c(93:125),]
maj = apply(allele[c(93:125),], 2, function(x) names(sort(summary(as.factor(x)), decreasing=T)[1]))
maj = matrix(maj, nrow(allele), nc = ncol(allele), byrow = T)
z = maj == allele
allele = z
allele[allele=="TRUE"] = "grey" #white is control maj allele
allele[allele!="grey"] = "black" #black is control min allele



plot(pos,y=rep(1,length(pos)), col=allele[1,],pch=15, ylim = c(1,92))
for(i in 2:92){
	points(pos,y=rep(i,length(pos)), col=allele[i,],pch=15)

}



x11();plot(pos,y=rep(1,length(pos)), col=allele[93,],pch=15, ylim = c(92,126))
for(i in 94:125){
	points(pos,y=rep(i,length(pos)), col=allele[i,],pch=15)

}


#3173, 3507, 3521
#3822, 3824, 3780

x11();plot(pos,y=rep(1,length(pos)), col=allele[38,],pch=15, ylim = c(1,126))
for(i in c(39,40,55,54,53,77:80,102:104, 105:107,93:95)){
	points(pos,y=rep(i,length(pos)), col=allele[i,],pch=15)

}






 [1]"3145-1" "3145-2" "3145-3" "3152-1" "3152-2" "3152-3" "3154-1"
 [8] "3154-2" "3154-3" "3155-1" "3155-2" "3155-3" "3156-1" "3156-2" "3156-3"
 [16] "3160-1" "3160-2" "3160-3" "3162-1" "3162-2" "3162-3" "3163-1" "3163-2"
 [24] "3163-3" "3164-1" "3164-2" "3164-3" "3165-1" "3165-2" "3165-3" "3169-1"
 [32] "3169-2" "3169-3" "3169-4" "3172-1" "3172-2" "3172-3" "3173-1" "3173-2"
 [40] "3173-3" "3175-1" "3175-2" "3175-3" "3498-1" "3498-2" "3498-3" "3504-1"
 [48] "3504-2" "3504-3" "3505-1" "3505-2" "3505-3" "3507-1" "3507-2" "3507-3"
 [56] "3508-1" "3508-2" "3508-3" "3510-1" "3510-2" "3510-3" "3516-1" "3516-2"
 [64] "3516-3" "3517-1" "3517-2" "3517-3" "3518-1" "3518-2" "3518-3" "3519-1"
 [72] "3519-2" "3519-3" "3520-1" "3520-2" "3520-3" "3521-1" "3521-2" "3521-3"
 [80] "3521-4" "3522-1" "3522-2" "3522-3" "3523-1" "3523-2" "3523-3" "3525-1"
 [88] "3525-2" "3525-3" "3774-1" "3774-2" "3774-3" 
 
[93]"3780-1" "3780-2" "3780-3"
[96] "3800-1" "3800-2" "3800-3" "3803-1" "3803-2" "3803-3" "3822-1" "3822-2"
[104] "3822-3" "3824-1" "3824-2" "3824-3" "3828-1" "3828-2" "3828-3" "3841-1"
[112] "3841-2" "3841-3" "3858-1" "3858-3" "3862-1" "3862-2" "3862-3" "3862-4"
[120] "3867-1" "3867-2" "3867-3" "3870-1" "3870-2" "3870-3"

[c(93:125),]






