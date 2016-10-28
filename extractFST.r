###
# Extract high FST sites
###





#!/usr/bin/Rscript


ma = read.table(file="MA.weir.fst",header=T,colClasses = c("character", "character", "numeric"))
cm = read.table(file="CM.weir.fst",header=T,colClasses = c("character", "character", "numeric")) 
ac = read.table(file="CA.weir.fst",header=T,colClasses = c("character", "character", "numeric")) 

names(ma)[3]= "ma"
names(ac)[3]= "ac"
names(cm)[3]= "cm"

#ma = ma[which(ma$fst>0.95),]
#ac = ac[which(ac$fst>0.95),]
#cm = cm[which(cm$fst>0.95),]

fst = ma 
fst$ac = ac$ac
fst$cm = cm$cm

fst = fst[which(fst$ma>=0),]
fst = fst[which(fst$cm>=0),]
fst = fst[which(fst$ac>=0),]
#fst = fst[fst$ma>0.9 | fst$ac>0.9,]

#ones  =rowSums(fst[3:5]=="1")
#fst = fst[ones>0,]
#aims = rbind(ma, cm, ac)
#aims = aims[!duplicated(aims$POS),]


write.table(fst[c(1,2)], file="aims", col.names=F, row.names=F, quote=F)
