##
# Analyze Hap Data
###



#load libraries ----------
library(IRanges)
source("/media/data1/forty3/brock/GenomeR/Q_correct.r")
source("/media/data1/forty3/brock/GenomeR/VarFunct.r")

#load dataframes -----------------------
	#cd /media/data1/forty3/drone/working
load(file="HapFLKGWAS.RData")
	#file contents:	
	#load(file="/media/data1/forty3/drone/PLINK/hapflk-1.3.0/chr/hapFLK.RData")
	#load(file="/media/data1/forty3/drone/hapstats/hapGWAS.RData") #win.size = 101 
	#load(file="/media/data1/forty3/drone/hapstats/hapGWAS_C.RData")
	#source("/media/data1/forty3/brock/GenomeR/Q_correct.r")

#Load Pi and GFF
pi = read.table(file="/media/data1/forty3/drone/vcf_drone/sel.Tajima.D",header=T, colClasses = c("character", "numeric","numeric","numeric"))
gff = read.table("/media/data1/forty3/R/GFF3_2",header=T, colClasses = c("character", "character", "numeric","numeric","character","character","character","character"))

# Mung Data frames ------------------------------------------
#create SNP ID
pi$SNP = paste(pi$CHROM, pi$BIN_START, sep="_")

#create scaff ID
pi$GRP = paste("Group", pi$CHROM, sep="")


#convert scaffold to chromosome
write.list(pi[c(5,6,2)],file="/media/data1/forty3/drone/ScaffConvert/CandidateMap")
system("perl /media/data1/forty3/drone/ScaffConvert/scaffold_to_chr_BAH.pl /media/data1/forty3/drone/ScaffConvert/scaffolds_on_chr.txt /media/data1/forty3/drone/ScaffConvert/CandidateMap")
conv = read.table(file="/media/data1/forty3/drone/ScaffConvert/CandidateMap_onChr.txt",header=F)
conv$V2 = gsub("chr","",conv$V2)
pi$CHROM = conv$V2
pi$start = conv$V3

	
	
	
	

#Calculate hapFLK significance cutoff via Storey -----------------
names(haps)[2] = "rs"
names(haps.C)[2] = "rs"
all.df = merge(haps, hap.df, by = "rs")
all.df = merge(all.df, haps.C, by = "rs")

x=(10^(-1*(all.df$p)))
q = qvalue1(x)$q
-log10(range(x[q<0.01]))
Q.cutoff = 5 #for Q<<0.01
rm(q,x)


	
#extract high regions and identify windows -----------
hi.df = all.df[all.df$p > Q.cutoff,]
hi.df$end = hi.df$pos + 10000
	#6914 SNPs?
hi.df$group2 = rep("NA", nrow(hi.df))

for(k in unique(hi.df$chr.x)){
	ir <- IRanges(hi.df$pos[which(hi.df$chr.x==k)], hi.df$end[which(hi.df$chr.x==k)]) 
	hi.df$group2[which(hi.df$chr.x==k)] <- subjectHits(findOverlaps(ir, reduce(ir)))
}

rang=aggregate(pos~chr.x + group2,data = hi.df, min)
rang.1=aggregate(pos~chr.x + group2, data = hi.df, function(x) max(x))
rang$end = rang.1$pos

rang = rang[order(rang$pos),]
rang = rang[order(rang$chr.x),]
write.list(rang, file="hapFLK_RANGES")

#Which have significant GWAS? -----------
gwas.rang = c()
chr.x = intersect(all.df$chr.x, rang$chr.x)
for(i in chr.x){
	win.temp = all.df[all.df$chr.x==i,]
	deg.temp = rang[which(as.character(rang$chr.x)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$pos), as.numeric(as.character(win.temp$pos)), "<=") 
	blah1=outer(as.numeric(deg.temp$end), as.numeric(as.character(win.temp$pos)), ">=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	gwas.rang = rbind(temp,gwas.rang)
	print(i)
}

nrow(gwas.rang[gwas.rang$V3.x>0 & gwas.rang$V3.y>0,])/nrow(gwas.rang)
head(gwas.rang[gwas.rang$V3.x>0 & gwas.rang$V3.y>0,])


GWAS.SNPs = (gwas.rang[gwas.rang$V3.x>0 & gwas.rang$V3.y>0,])
GWAS.SNPs = GWAS.SNPs[c(1,2,3,4,14)]; names(GWAS.SNPs)[c(3,5)] = c("start","pos")


#plot HapFLK ---------------------------------
source("hapFLKPlot.r")
	#this plots the Pvalue of hapflk for each chromosome and also plots the QTL regions that are known.
	#"pplotHygiene.tiff

	
	
#GWAS Plots -----------------------------------	
	



	
	
#get out genes -------------
all.df$chrom = gsub(":.*","",all.df$rs)
all.df$gpos = gsub(".*:","",all.df$rs)


gff = gff[which(gff$type=="gene"),]
#gff$chrom = paste("Group", gff$chrom , sep="") 

hi.genes = c()
chrom = intersect(all.df$chrom, gff$chrom)
for(i in chrom){
	win.temp = all.df[all.df$chrom==i,]
	deg.temp = gff[which(as.character(gff$chrom)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$start), as.numeric(as.character(win.temp$gpos)), "<=") 
	blah1=outer(as.numeric(deg.temp$end), as.numeric(as.character(win.temp$gpos)), ">=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	hi.genes = rbind(temp,hi.genes)
	print(i)
}

high.df = hi.genes[hi.genes$p>Q.cutoff,]


#GO Analysis --------------------
	#used DAVID with fly terms - Still all neurological, developmental
	


# TD and Fst for Selected regions/SNPs -------------

Pi.Fst = c()
pi$GRP = gsub("Group","",pi$GRP)
chrom = intersect(all.df$chrom, pi$GRP)
for(i in chrom){
	win.temp = all.df[all.df$chrom==i,]
	deg.temp = pi[which(as.character(pi$GRP)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$BIN_START)+1000, as.numeric(as.character(win.temp$gpos)), ">=") 
	blah1=outer(as.numeric(deg.temp$BIN_START), as.numeric(as.character(win.temp$gpos)), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	Pi.Fst = rbind(temp,Pi.Fst)
	print(i)
}

#significant TD and gros more significant as you increase TD
 t.test(Pi.Fst$TajimaD, Pi.Fst$TajimaD[Pi.Fst$p>5])

        Welch Two Sample t-test

data:  Pi.Fst$TajimaD and Pi.Fst$TajimaD[Pi.Fst$p > 5]
t = 3.1678, df = 6927.1, p-value = 0.001543
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.01471635 0.06249996
sample estimates:
mean of x mean of y
0.3771302 0.3385221


 t.test(Pi.Fst$TajimaD, Pi.Fst$TajimaD[Pi.Fst$p>9])

        Welch Two Sample t-test

data:  Pi.Fst$TajimaD and Pi.Fst$TajimaD[Pi.Fst$p > 9]
t = 11.188, df = 930.66, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.3250766 0.4633763
sample estimates:
  mean of x   mean of y
 0.37713024 -0.01709623

 
 
 
#Differential Admixture ------------------------------
	#obtained from LAMP.sh, etc. 
lamp = read.table(file="LAMP.summary")
names(lamp) = c("pos", "chr", "sel.c", "con.c", "sel.m", "con.m", "sel.a", "con.a", "sel.snp")

 
 
 
t.test(lamp$sel.c,lamp$sel.c[lamp$sel.snp=="1"])

        Welch Two Sample t-test

data:  lamp$sel.c and lamp$sel.c[lamp$sel.snp == "1"]
t = -5.9499, df = 688.11, p-value = 4.271e-09
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.05383107 -0.02711856
sample estimates:
mean of x mean of y
0.6445355 0.6850103

#sel bees have more "C" in hyg loci than rest of genome; con bees do not
#and more c lineage than con bees at same sites.
 
 
 
 
 
 
 
 
 
 
 







#
Pi.genic = c()
chrom = intersect(gff$chrom, pi$GRP)
for(i in chrom){
	win.temp = gff[gff$chrom==i,]
	deg.temp = pi[which(as.character(pi$GRP)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$BIN_START)+1000, as.numeric(as.character(win.temp$start)), ">=") 
	blah1=outer(as.numeric(deg.temp$BIN_START), as.numeric(as.character(win.temp$end)), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	Pi.genic = rbind(temp,Pi.genic)
	print(i)
}


t.test(Pi.Fst$TajimaD, Pi.genic$TajimaD)
 t.test(Pi.Fst$TajimaD, Pi.genic$TajimaD)

        Welch Two Sample t-test

data:  Pi.Fst$TajimaD and Pi.genic$TajimaD
t = -3.1492, df = 646.57, p-value = 0.001713
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.26918210 -0.06241697
sample estimates:
 mean of x  mean of y
0.08766019 0.25345972




# Fst in Selected in high regions ----------------------------------------
fst$FSTSEL = as.numeric(fst$FSTSEL)
FSTSEL  = fst$FSTSEL
FSTSEL[FSTSEL < 0] = 0
fst$FSTSEL = FSTSEL

wilcox.test(fst$FSTSEL,  fst$FSTSEL[fst$FSTP<0.05])


hi.Fst = c()
chrom = intersect(regions$GRP, fst$GRP)
for(i in chrom){
	win.temp = regions[regions$GRP==i,]
	deg.temp = fst[which(as.character(fst$GRP)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$GPOS)), ">=") 
	blah1=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$end)), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	hi.Fst = rbind(temp,hi.Fst)
	print(i)
}

boxplot(-log10(hi.Fst$FSTSEL), -log10(fst$FSTSEL[!(fst$SNP %in% hi.Fst$SNP)]), notch=T)
t.test(-log10(hi.Fst$FSTSEL), -log10(fst$FSTSEL[!(fst$SNP %in% hi.Fst$SNP)]))
wilcox.test(-log10(hi.Fst$FSTSEL), -log10(fst$FSTSEL[!(fst$SNP %in% hi.Fst$SNP)]))

































#differential Admixture ------------------------
	#LAMP.sh for thsi analysis, outputting through LAMPanalysis.r


	

	
	
	
	
	
	
	

#Overlap with DEGs? ----------------------------------------
degs = read.table(file="/media/data1/forty3/drone/git/data/boutinDEGs.txt",header=T)
degs[which(degs$GB %in% hi.genes$GB),]
# (3 or 0)


#Overlap with DEPs? ----------------------------------------
deps = read.table(file="/media/data1/forty3/drone/git/data/FosterDeps.txt",header=T)
deps[which(deps$GB %in% hi.genes$GB),]
#none.

	
	
	
	
	
	
	
	
	
	
	
	











###########################################################
# i = 5
#all.df = hap.df
#hap.df = all.df[all.df$chr==i,]
#hap.df[hap.df$p>10,] #tweak

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
haps.df$runmed = runmed((haps.df$x),k=101) #tweak

test = merge(all.df, haps.df, by = "rs")
test = test[order(test$pos),]


plot(test$pos, 3+test$p, ylim = c(0,19))
points(test$pos[test$runmed.y>0.1], test$runmed.y[test$runmed.y>0.1])
points(test$pos[test$runmed.y>0.2], 1+test$runmed.y[test$runmed.y>0.2])











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
	


