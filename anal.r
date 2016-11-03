##
# Analyze Hap Data
###



#load libraries ----------
library(IRanges)
source("/media/data1/forty3/brock/GenomeR/Q_correct.r")
source("/media/data1/forty3/brock/GenomeR/VarFunct.r")

#load dataframes -----------------------
	#cd /media/data1/forty3/drone/working
load(file="working.RData") #load(file="HapFLKGWAS.RData")

	#file contents:	
	#load(file="/media/data1/forty3/drone/PLINK/hapflk-1.3.0/chr/hapFLK.RData")
	#load(file="/media/data1/forty3/drone/hapstats/hapGWAS.RData") #win.size = 101 
	#load(file="/media/data1/forty3/drone/hapstats/hapGWAS_C.RData")
	#source("/media/data1/forty3/brock/GenomeR/Q_correct.r")

#Load Pi and GFF
pi = read.table(file="/media/data1/forty3/drone/vcf_drone/sel.Tajima.D",header=T, colClasses = c("character", "numeric","numeric","numeric"))
gff = read.table("/media/data1/forty3/R/GFF3_2",header=T, colClasses = c("character", "character", "numeric","numeric","character","character","character","character"))


#load in immune gene list 
imm = read.table("/media/data1/forty3/brock/immune/EvansNewImm.txt", header =T)


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

###
#Selection analysis
########################################################	
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
hi.df$end = hi.df$pos + 500
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
write.list(hi.df$rs, file="hapFLK_SNPS")

#plot HapFLK ---------------------------------
source("hapFLKPlot.r")
	#this plots the Pvalue of hapflk for each chromosome and also plots the QTL regions that are known.
	#"pplotHygiene.tiff

	
#MAF and HWE within selected regions------------------
frq = read.table(file="/media/data1/forty3/drone/vcf_drone/S.frq",header=F,colClasses = c("character","numeric", "numeric", "character","character"),skip = 1)

MAF = as.numeric(gsub("*.[:]","",frq$V6))
MAF[MAF >= 0.5] = (1 - MAF[MAF >= 0.5])
frq$MAF = MAF

frq$rs = paste(frq$V1, frq$V2, sep=":")
test = merge(frq, all.df, by = "rs")	
	
	
	
hwe = read.table(file="/media/data1/forty3/drone/vcf_drone/S.hwe",header=F,colClasses = c("character","numeric",  "character","character","character","numeric"),skip = 1)
hwe$V6 = -log10(hwe$V6 + 1e-06 )

hwe$rs = paste(hwe$V1, hwe$V2, sep=":")
test = merge(hwe, all.df, by = "rs")	
boxplot(test$V6, test$V6[test$p>5])		

#significantly more homozygotes in selected windows. 	
	# t.test(test$V6, test$V6[test$p>5])
	#
	#       Welch Two Sample t-test
	#data:  test$V6 and test$V6[test$p > 5]
	#t = -4.2485, df = 5342.6, p-value = 2.188e-05


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
#
#        Welch Two Sample t-test
#
#data:  Pi.Fst$TajimaD and Pi.Fst$TajimaD[Pi.Fst$p > 5]
#t = 3.1678, df = 6927.1, p-value = 0.001543
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
# 0.01471635 0.06249996
#sample estimates:
#mean of x mean of y
#0.3771302 0.3385221
#
#
t.test(Pi.Fst$TajimaD, Pi.Fst$TajimaD[Pi.Fst$p>9])
#
#        Welch Two Sample t-test
#
#data:  Pi.Fst$TajimaD and Pi.Fst$TajimaD[Pi.Fst$p > 9]
#t = 11.188, df = 930.66, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
# 0.3250766 0.4633763
#sample estimates:
#  mean of x   mean of y
# 0.37713024 -0.01709623
#
#



###
#SNP Characterization
########################################################	

hap.snps = read.table("CAND-up_dwn.snpeff.eff", header=F, skip = 3, sep="\t", colClasses = c(rep("character", 24)))
hap.snps$rs = paste(hap.snps$V1, hap.snps$V2, sep=":")
hap.snps$V16 = gsub("[:].*", "", hap.snps$V16)
#snps = snps[which(snps$V16!="INTERGENIC"),]
#snps = snps[which(snps$V16!="INTRAGENIC"),]
#snps = snps[which(snps$V16!="SYNONYMOUS_STOP"),]
write.list(aggregate(as.factor(V16)~rs, data=hap.snps, table),file="hapSNPSummary")
hap.summ = read.table(file="hapSNPSummary",header=F)
names(hap.summ) = c("GB",names(table(hap.snps$V16)))

#all SNPs
snps = read.table("HYG-up_dwn.snpeff.eff", header=F, skip = 3, sep="\t", colClasses = c(rep("character", 24)))
snps$rs = paste(snps$V1, snps$V2, sep=":")
snps$V16 = gsub("[:].*", "", snps$V16)
write.list(aggregate(as.factor(V16)~rs, data=snps, table),file="hapSNPSummary")
summ = read.table(file="hapSNPSummary",header=F)
names(summ) = c("GB",names(table(snps$V16)))

















write.table(table(test$V16[test$V10 %in% imm$GB]))
write.table(table(test$V16[test$V10 %in% hap.summ$GB]))


write.list(aggregate(as.factor(V16[test$p>5])~V10[test$p>5], data=test, table), file = "NSYNSYMMARY")
hap.summ = read.table(file="NSYNSYMMARY",header=F)
names(hap.summ) = c("GB",names(table(test$V16[test$p>5])))

write.list(aggregate(as.factor(V16[test$p<5])~V10[test$p<5], data=test, table), file = "ALLNSYNSYMMARY")
summ = read.table(file="ALLNSYNSYMMARY",header=F)
names(summ) = c("GB",names(table(test$V16[test$p<5])))


hap.snps = snps[which(snps$V11 %in% hap.summ$GB),]













###
#Gene Characterization
########################################################	

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

high.genes = unique(high.df$GB)




#Overlap with DEGs? ----------------------------------------
degs = read.table(file="/media/data1/forty3/drone/git/data/boutinDEGs.txt",header=T)
degs[which(degs$GB %in% high.genes),]
#        Gene Chromosome      GB Reg
#11 LOC552229          1 GB46514 DWN


#Overlap with DEPs? ----------------------------------------
deps = read.table(file="/media/data1/forty3/drone/git/data/FosterDeps.txt",header=T)
deps[which(deps$GB %in% high.genes),]
#none.


#overlap with Beye? -------------------------------------------
beye = read.table("clipboard", header=F)
beye[which(beye$V1 %in% high.genes),]
	#GB43119	FBgn0038035	lig3 (11.20)
	#GB44556	FBgn0027609	morgue (5.14)
	#GB53743	0	LOC408740 (3.9)
	#GB54885	FBgn0004859	Ci (4.4)








# Developmental biased genes? ------------------------
	#I'm pulling in the supplemental table S5 from Pires et al, 2016)
	#http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0146447#sec022
	#Using this, I want to see if my 112 genes are expressed during specific LH stages

#C1 = only detected in mature oocytes or transcripts expressed in mature oocytes that had decreased expression until 18–24 h.
#C2 = Class II included mRNAs and miRNAs that were undetected in mature oocytes but had significant expression in the 0–2 h, 0–6 h and 18–24 h 
#C3 =  expressed in mature oocytes with decreased (or absent) expression during cleavage (0–2 h and 0–6 h) and increased expression in the 18–24 h embryos.	
	 #The transcripts of class III mRNAs in the haploid and diploid embryos were involved in the process of “cell division.”
exp = read.table(file="/media/data1/forty3/drone/git/data/PiresEtalS5.txt",header=T)	
	
table(exp$pattern[which(exp$GB %in% high.genes)])
#C1 C2 C3
# 9 11  2

table(exp$pattern)
#  C1   C2   C3
# 656 2680  263




# Old versus new genes? ------------------------		
trg = read.table(file="/media/data1/forty3/drone/git/data/PNAS_TRG.txt",header=T) #this is straigth from my PNAS paper (SD5)
test = trg 
qwd = rep("N", nrow(test)); 
qwd[test$GB %in% high.genes] = "sig" 
test$sig = qwd
aggregate(test$sig, by = list(test$Taxa), table)	
#      Group.1  x.N x.sig
#1        Apis   88     0
#2     Apoidea  214     1
#3 Hymenoptera 1306    15
#4     Insecta 8620    66
	
# Expression? ------------------------		
	#Use Johnson's Data (Johnson_TableS3.xlsx)
exp = read.table(file="/media/data1/forty3/drone/git/data/johnson_expression.txt",header=T)
test = exp 
qwd = rep("N", nrow(test)); 
qwd[test$GB %in% high.genes] = "sig" 
test$sig = qwd	
aggregate(test$sig, by = list(test$Taxonomy), table)		
#      Group.1  x.N x.sig
#1   Arthropod  211     3
#2         Bee   24     0
#3   Conserved 5377    37
#4 Hymenoptera  389     3
#5      Insect  561     2
#6      Orphan  191     0
#7      O-TRGs  800     7




#load datasets (Pi, TD, gamma)
gamm = read.table(file = "/media/data1/forty3/brock/immune/Gene_gamma.txt", header=T)
kaks = read.table(file="/media/data1/forty3/brock/immune/KAKS.out")
names(kaks) =  c("GB", "c.ks", "m.ks", "y.ks", "a.ks", "c.ka", "m.ka", "y.ka", "a.ka" )
kaks$kaks = kaks$c.ka/kaks$c.ks
kaks = kaks[kaks$c.ks>0,]
Pi = read.table(file="/media/data1/forty3/brock/immune/C.snpeff1.eff.Pi")
names(Pi) = c("GB", "NS","NN", "NumS", "NumN", "PiS", "PiN")


test = merge(Pi, kaks, by = "GB")



	
#Gamma
test = gamm ; qwd = rep("N", nrow(test)); 
qwd[test$GB %in%  high.genes] = "HYG" 
test$qwd = qwd



#Gamma vs immune 
imm = read.table("/media/data1/forty3/brock/immune/EvansNewImm.txt", header =T)
#imm.gen = imm[imm$Class %in% c("Recognition", "Signalling","Effector"),]

test = gamm ; qwd = rep("N", nrow(test)); 
qwd[test$GB %in%  high.genes] = "HYG" 

qwd[test$GB %in%  imm$GB] = "IMM" 
test$qwd = qwd

boxplot(test$gamma~test$qwd)

pairwise.wilcox.test(test$gamma,test$qwd)



 
#Selection in C and M lineages ------------
all.df$chrom = gsub(":.*","",all.df$rs)
all.df$gpos = gsub(".*:","",all.df$rs)

# TD 
ctd = read.table(file="/media/data1/forty3/drone/vcf_drone/C.Tajima.D",header=T, colClasses = c("character", "numeric", "numeric","numeric"))
ntd = read.table(file="/media/data1/forty3/drone/vcf_drone/N.Tajima.D",header=T, colClasses = c("character", "numeric", "numeric","numeric"))
mtd = read.table(file="/media/data1/forty3/drone/vcf_drone/M.Tajima.D",header=T, colClasses = c("character", "numeric", "numeric","numeric"))



# TD and Fst for Selected regions/SNPs -------------
mtd.gene = c()
chrom = intersect(gff$chrom, mtd$CHROM)
for(i in chrom){
	win.temp = gff[gff$chrom==i,]
	deg.temp = mtd[which(as.character(mtd$CHROM)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$BIN_START)+1000, as.numeric(as.character(win.temp$end)), "<=") 
	blah1=outer(as.numeric(deg.temp$BIN_START), as.numeric(as.character(win.temp$start)), ">=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	mtd.gene  = rbind(temp,mtd.gene  )
	print(i)
}



ctd.gene = c()
chrom = intersect(gff$chrom, ctd$CHROM)
for(i in chrom){
	win.temp = gff[gff$chrom==i,]
	deg.temp = ctd[which(as.character(ctd$CHROM)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$BIN_START)+1000, as.numeric(as.character(win.temp$end)), "<=") 
	blah1=outer(as.numeric(deg.temp$BIN_START), as.numeric(as.character(win.temp$start)), ">=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	ctd.gene  = rbind(temp,ctd.gene  )
	print(i)
}


ntd.gene = c()
chrom = intersect(gff$chrom, ntd$CHROM)
for(i in chrom){
	win.temp = gff[gff$chrom==i,]
	deg.temp = ntd[which(as.character(ntd$CHROM)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$BIN_START)+1000, as.numeric(as.character(win.temp$end)), "<=") 
	blah1=outer(as.numeric(deg.temp$BIN_START), as.numeric(as.character(win.temp$start)), ">=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	ntd.gene  = rbind(temp,ntd.gene  )
	print(i)
}



ctd = aggregate(ctd.gene$TajimaD, by = list(ctd.gene$GB ), mean)
mtd = aggregate(mtd.gene$TajimaD, by = list(mtd.gene$GB ), mean)
ntd = aggregate(ntd.gene$TajimaD, by = list(ntd.gene$GB ), mean)



qwd = rep("N", nrow(ctd)); 
qwd[ctd$Group.1 %in% high.genes] = "sig" 
ctd$sig = qwd

qwd = rep("N", nrow(mtd)); 
qwd[mtd$Group.1 %in% high.genes] = "sig" 
mtd$sig = qwd

qwd = rep("N", nrow(ntd)); 
qwd[ntd$Group.1 %in% higher.genes] = "sig" 
ntd$sig = qwd


nrow(mtd[mtd$x<median(mtd$x) & mtd$sig == "sig" ,])/nrow(mtd[ mtd$sig == "sig" ,])
nrow(ctd[ctd$x<median(ctd$x) & ctd$sig == "sig" ,])/nrow(ctd[ ctd$sig == "sig" ,])
nrow(ntd[ntd$x<median(ntd$x) & ntd$sig == "sig" ,])/nrow(ntd[ ntd$sig == "sig" ,])

#test against immune genes
qwd = rep("N", nrow(ctd)); 
qwd[ctd$Group.1 %in% imm$GB] = "sig" 
ctd$sig = qwd

qwd = rep("N", nrow(mtd)); 
qwd[mtd$Group.1 %in% imm$GB] = "sig" 
mtd$sig = qwd


nrow(mtd[mtd$x<median(mtd$x) & mtd$sig == "sig" ,])/nrow(mtd[ mtd$sig == "sig" ,])
nrow(ctd[ctd$x<median(ctd$x) & ctd$sig == "sig" ,])/nrow(ctd[ ctd$sig == "sig" ,])

nrow(mtd[mtd$x<median(mtd$x),])/nrow(mtd)
nrow(ctd[ctd$x<median(ctd$x),])/nrow(ctd)
nrow(ntd[ntd$x<median(ntd$x),])/nrow(ntd)

#selection in CANADA
    #Group.1           x sig
#1311 GB42184  0.02824837 sig
#3572 GB45765 -0.55997950 sig
#3577 GB45774 -1.02116309 sig
#5436 GB49098 -0.04325323 sig


#Compare FST in regions among AMC ---------------

				amc.fst = All.Data
				for(fil in dir(pattern="/media/data1/forty3/drone/vcf_drone/*.weir.fst")){
					#fil = "A_vs_CONT.weir.fst"
					temp_fst = read.table(file=fil,header=T, colClasses = c("character", "numeric", "numeric"))
					temp_fst$SNP = paste(temp_fst$CHROM, "_", temp_fst$POS,sep="")
					temp_fst = temp_fst[-c(1,2)]
					fix.fst  = temp_fst$WEIR_AND_COCKERHAM_FST
					fix.fst[fix.fst<0]= 0
					temp_fst$WEIR_AND_COCKERHAM_FST = fix.fst
					fil = gsub(".weir.fst","",fil)
					names(temp_fst)=c(fil,"SNP")
					amc.fst = merge(amc.fst, temp_fst, by = "SNP")
					print(fil)
				}
				rm(fil, temp_fst)

				amc.fst$C.CON = amc.fst$C_vs_CONT/(amc.fst$C_vs_CONT + amc.fst$M_vs_CONT)
				amc.fst$C.SEL = amc.fst$C_vs_SEL/(amc.fst$C_vs_SEL + amc.fst$M_vs_SEL)



				AMC.high.fst = c()
				chrom = intersect(rang$chrom, amc.fst$CHROM)
				for(i in chrom){
					win.temp = rang[rang$chrom==i,]
					deg.temp = amc.fst[which(as.character(amc.fst$CHROM)==as.character(i)),]
					blah=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$start)), ">=") 
					blah1=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$end)), "<=") 
					blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
					temp = deg.temp[blah[,1],]
					temp = cbind(temp, win.temp[blah[,2],])
					AMC.high.fst  = rbind(temp,AMC.high.fst )
					print(i)
				}

				rm(temp)




				x11();boxplot(AMC.high.fst$M_vs_SEL, amc.fst$M_vs_SEL, AMC.high.fst$M_vs_CONT,amc.fst$M_vs_CONT, notch = T)
				x11();boxplot(AMC.high.fst$C_vs_SEL, amc.fst$C_vs_SEL, AMC.high.fst$C_vs_CONT,amc.fst$C_vs_CONT, notch = T)

				x11();boxplot(AMC.high.fst$M_vs_C, amc.fst$M_vs_C) #


				#Judging by this, it seems Selected bees at these regions have become "more C" and less "M" than control. 


				x11();boxplot(AMC.high.fst$C.CON, amc.fst$C.SEL, notch = T)














###
#GWAS analysis
########################################################	

#Which have significant GWAS? -----------
hi.df = all.df[all.df$p > Q.cutoff,]




gwas.rang = c()
chr.x = intersect(all.df$chr.x, rang$chr.x)
for(i in chr.x){
	win.temp = all.df[all.df$chr.x==i,]
	deg.temp = rang[which(as.character(rang$chr.x)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$pos)-1000, as.numeric(as.character(win.temp$pos)), "<=") 
	blah1=outer(as.numeric(deg.temp$end)+1000, as.numeric(as.character(win.temp$pos)), ">=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	gwas.rang = rbind(temp,gwas.rang)
	print(i)
}


nrow(gwas.rang[gwas.rang$V3.y>0,])/nrow(gwas.rang)
nrow(gwas.rang[gwas.rang$V3.y>0 & gwas.rang$p>5,])/nrow(gwas.rang)





nrow(gwas.rang[gwas.rang$V3.x>0 & gwas.rang$V3.y>0,])/nrow(gwas.rang)
head(gwas.rang[gwas.rang$V3.x>0 & gwas.rang$V3.y>0,])



GWAS.SNPs = (gwas.rang[gwas.rang$V3.x>0 & gwas.rang$V3.y>0,])
GWAS.SNPs = GWAS.SNPs[c(1,2,3,4,14)]; names(GWAS.SNPs)[c(3,5)] = c("start","pos")

		

	
#GWAS Plots -----------------------------------	

	#1) as the correlation between trait and genotype increases, the effect of differentiation following selection increases, fixing alleles in selected populations (N&P 1997). In the regions of the genome we've found acted on by selection, MAF drops significantly relative to the rest of the genome (test) not allowing us to fully identify candidate robustly associated in a GWAS within selected populations. However, we find that N regions are near to selected
	#2) within the unselected population, we are able to detect an effect of X

#Gwas power? -----------------------------------
	#gwas_power.r























	
	
	



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
	


 
 
 
#Differential Admixture ------------------------------
	#obtained from LAMP.sh, etc. 
lamp = read.table(file="LAMP_hygAMC.summary")
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
 
 
 
 
 
 
 



































############################################



 







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
	


