#!/usr/bin/Rscript
###
# Concatenate and assess significance of haplotype scanned GWAS
###
	

#assumes file names are consistent with previous output hap_<CHR#>.qassoc.hap.snp (sorry)	
	

#required packages --------------
require(MASS)
	
#take in arguments ---------------------------
args <- commandArgs(trailingOnly = TRUE)
chrs = args[1] #16
directs = args[2] #/media/data1/forty3/drone/hapstats/
Pvalue = args[3] #0.05
win.size = arg[4] #1041


#run through each file, concatenate ------------------
hap.df = c()
for(i in c(1:chrs)){
	#1=5
	hap = read.table(file=paste(directs,"hapflk_",i,".hapflk",sep=""),header=T)
	hap.df = rbind(hap, hap.df)
	hap.df = hap.df[order(hap.df$chr),]
}



#calculate Pvalues ------------------------------
mod=rlm(hap.df$hapflk~1)
mu=mod$coefficients[1]
ss=mod$s
pvalue=1-pnorm(hap.df$hapflk,mean=mu,sd=ss)
hap.df$p = -log10(pvalue)

# clean and save --------------------------------
rm(list=setdiff(ls(), "hap.df"))
save.image(file="hapFLK.RData")

























###########################################################
# i = 6
#all.df = hap.df
hap.df = all.df[all.df$chr==i,]
hap.df[hap.df$p>10,] #tweak
P = 0.05
haps = read.table(file=paste("/media/data1/forty3/drone/hapstats/hap_",i,".qassoc.hap.snp",sep=""), header=F, skip = 1)
haps$V8=NULL
haps.df = aggregate(haps$V7, by = list(haps$V1, haps$V9), function(x) length(x[which(as.numeric(x)<P)])/length(x))
haps.df$len = aggregate(haps$V7, by = list(haps$V1, haps$V9), function(x) length(x[!is.na(as.numeric(x))]))$x
haps.df$sig.len = aggregate(haps$V7, by = list(haps$V1, haps$V9), function(x) length(x[which(as.numeric(x)<P)]))$x
haps.df$sig.mean = aggregate(haps$V7, by = list(haps$V1, haps$V9), function(x) mean(x))$x
haps.df$WIN = gsub("WIN","",haps.df$Group.1)
haps.df = haps.df[order(as.numeric(haps.df$WIN)),]

names(haps.df)[2] = "rs"

P = 0.05
haps.df$runmed = runmed((haps.df$x),k=1041) #tweak
plot(haps.df$runmed)


test = merge(haps.df, hap.df, by = "rs")
test = test[order(as.numeric(test$pos)),]



plot(test$pos, 3+test$p, ylim = c(0,18))
points(test$pos[test$runmed>0.1], (test$runmed[test$runmed>0.1]))
abline(h=10, col="red")
test[test$p>10,]

	#i = 5
	
	#SNP = apply(hap[8], 1, function(x) unlist(strsplit(as.character(x), "[|]"))[5])
	#hap$SNP = SNP
	hap.df = aggregate(hap$V7, by = list(hap$V1, hap$V9), function(x) length(x[which(as.numeric(x)<P)])/length(x))
	hap.df$len = aggregate(hap$V7, by = list(hap$V1, hap$V9), function(x) length(x[!is.na(as.numeric(x))]))$x
	hap.df$sig.len = aggregate(hap$V7, by = list(hap$V1, hap$V9), function(x) length(x[which(as.numeric(x)<P)]))$x
	#hap.df$mid.SNP = aggregate(hap$SNP, by = list(hap$LOCUS), function(x) x[1])
	hap.df$WIN = gsub("WIN","",hap.df$Group.1)
	hap.df = hap.df[order(as.numeric(hap.df$WIN)),]
	write.table(hap.df, file="haplen001.summary", append=T, col.names=F, row.names=F, quote=F)
	print(i)







haps = read.table(file="/media/data1/forty3/drone/hapstats/haplen001.summary", header=F)
haps$chr = gsub("[.].*","",haps$V2) 
names(haps)[2] = "rs"



test = merge(hap.df, haps, by ="rs")
test = test[order(test$pos),]
test = test[order(test$chr.x),]

plot(3+test$p, col=test$chr.x,ylim=c(0,19))
points(x=rep(1:nrow(test)), y=test$runmed, col=test$chr.x)



p.rank = rank(test$p, na.last ="keep")
f.rank = rank(test$V3, na.last ="keep")


p.rank01 = p.rank/ (max(p.rank,    na.rm = T) + 1)
f.rank01 = f.rank/ (max(f.rank,    na.rm = T) + 1)

p.z = qnorm(p.rank01) #I THINK?
f.z = qnorm(f.rank01)


test.mz <- (p.z + f.z)/2
p.test.mz <- pnorm(test.mz, 0, 1/sqrt(2), lower.tail = F) #  zbar ~ N(0, 1/sqrt(3))

css = -log10(p.test.mz)

test$css = css
test = test[order(test$POS),]
test = test[order(c(test$CHR, test$POS)),]






names(hapsC)[2] = "rs"
test = merge(test, hapsC, by ="rs")
	


