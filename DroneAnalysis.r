###
# FST analyses
###






#See README, I estimated Fst as Weir & Cockerham's Fst for two populations and pFst a likelihood ratio test for allele frequency differences between populations.


# Load functions and dataframes ------------------------------
library(IRanges)
library(Hmisc)
library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)
library(scales)
library(data.table)
source("/media/data1/forty3/drone/git/GenomeR/VarFunct.r")


#Load Fst and Pi and GFF
pfst = read.table(file="/media/data1/forty3/drone/FST/pFST/pFST.out",header=F, colClasses = c("character", "numeric", "numeric"))
pfst.fas = read.table(file="/media/data1/forty3/drone/FST/pFST/pFST.FAS.out",header=F, colClasses = c("character", "numeric", "numeric"))
pfst.mas = read.table(file="/media/data1/forty3/drone/FST/pFST/pFST.MAS.out",header=F, colClasses = c("character", "numeric", "numeric"))
fst = read.table(file="/media/data1/forty3/drone/FST/pFST/wcFST.out",header=F, colClasses = c("character", "numeric", "numeric", "numeric"))
pi = read.table(file="/media/data1/forty3/drone/vcf_drone/sel.Tajima.D",header=T, colClasses = c("character", "numeric","numeric","numeric"))
gff = read.table("/media/data1/forty3/R/GFF3_2",header=T, colClasses = c("character", "character", "numeric","numeric","character","character","character","character"))


# Mung Data frames ------------------------------------------
	
#create SNP ID
pfst$SNP = paste(pfst$V1, pfst$V2, sep="_")
pfst.mas$SNP = paste(pfst.mas$V1, pfst.mas$V2, sep="_")
pfst.fas$SNP = paste(pfst.fas$V1, pfst.fas$V2, sep="_")
fst$SNP = paste(fst$V1, fst$V2, sep="_")
pi$SNP = paste(pi$CHROM, pi$BIN_START, sep="_")

#create scaff ID
fst$V1 = paste("Group", fst$V1, sep="")
pfst$V1 = paste("Group", pfst$V1, sep="")
pfst.mas$V1 = paste("Group", pfst.mas$V1, sep="")
pfst.fas$V1 = paste("Group", pfst.fas$V1, sep="")
pi$GRP = paste("Group", pi$CHROM, sep="")

#merge
pfst.mas= pfst.mas[c(3,4)];names(pfst.mas)[1] = "pMAS"
pfst.fas = pfst.fas[c(3,4)];names(pfst.fas)[1] = "pFAS"
pfst = pfst[c(3,4)];names(pfst)[1] = "P"


fst = data.table(fst)
pfst.mas = data.table(pfst.mas)
pfst.fas = data.table(pfst.fas)
pfst = data.table(pfst)
fst = merge(fst, pfst.mas, by="SNP")
fst = merge(fst, pfst.fas, by="SNP")
fst = merge(fst, pfst, by="SNP")
fst = data.frame(fst)

fst = fst[c(1,2,3,6,7,8,9)]; names(fst) = c("SNP","GRP","GPOS","FST","pMAS", "pFAS", "FSTP")


#convert scaffold to chromosome
write.list(fst[c(1,2,3)],file="/media/data1/forty3/drone/ScaffConvert/CandidateMap")
system("perl /media/data1/forty3/drone/ScaffConvert/scaffold_to_chr_BAH.pl /media/data1/forty3/drone/ScaffConvert/scaffolds_on_chr.txt /media/data1/forty3/drone/ScaffConvert/CandidateMap")
conv = read.table(file="/media/data1/forty3/drone/ScaffConvert/CandidateMap_onChr.txt",header=F)
conv$V2 = gsub("chr","",conv$V2)
fst$CHROM = conv$V2
fst$POS = conv$V3


#convert scaffold to chromosome
write.list(pi[c(6,7,2,3)],file="/media/data1/forty3/drone/ScaffConvert/CandidateMap")
system("perl /media/data1/forty3/drone/ScaffConvert/scaffold_to_chr_BAH.pl /media/data1/forty3/drone/ScaffConvert/scaffolds_on_chr.txt /media/data1/forty3/drone/ScaffConvert/CandidateMap")
conv = read.table(file="/media/data1/forty3/drone/ScaffConvert/CandidateMap_onChr.txt",header=F)
conv$V2 = gsub("chr","",conv$V2)
pi$CHROM = conv$V2
pi$start = conv$V3
pi$end = conv$V4




# Analyses -----------------------------
#All.Data = fst
test = All.Data

#Identify High Sites with logP
test$pMAS = -1*log10(test$pMAS)
test$pFAS = -1*log10(test$pFAS)
test$FSTP = -1*log10(test$FSTP)

#Identify clusters of high FST regions
	#used creeping window approach
window.size = 1000 #size of window to be crept and gaps to skip, in BP 
bin.size = 5 # minimum number of SNPs in a window (bin size, see plot)
test = test[with(test, order(POS)),]

creep.all = c()
for(k in unique(test$CHROM)){
	
	#Gather SNPs and FST as numeric variables
	snp = as.numeric(test$POS[test$CHROM==k])
	fst = as.numeric(test$pMAS[test$CHROM==k]) 
	fst2 = as.numeric(test$pFAS[test$CHROM==k])
	fst3 = as.numeric(test$FSTP[test$CHROM==k]) 	
	#Pc = as.numeric(test$P[test$CHROM==k]) 
	#dp = as.numeric(test$MEAN_DEPTH[test$CHROM==k]) 
	
	#fst in Hvs L
	#estimate difference between SNPs	
	p = as.numeric(snp)
	p = p[-1]
	p = c(p, NA)
	diff = abs(as.numeric(snp)-p)
	diff = diff[!is.na(diff)]

	x=rep(0, (max(snp)+ window.size))
	x[snp]=1
	x=cumsum(x)
	endsnp=x[snp+ window.size]
	n=length(fst)
	vec=rep(0,n)
	vec2=rep(0,n)
	vec3=rep(0,n)
	vec4=rep(0,n)
		
	len=rep(0,n)
	end=rep(0,n)
	
	for (i in 1:(n-1)){
			if (diff[i] < window.size ){
				snp_len=length(fst[i:endsnp[i]])
				
				fst_mean=mean(fst[i:endsnp[i]])
				vec[i]=fst_mean				
				
				fst_mean=mean(fst2[i:endsnp[i]])
				vec2[i]=fst_mean				
				
				fst_mean=mean(fst3[i:endsnp[i]])
				vec3[i]=fst_mean					

				
				len[i]=snp_len
				end[i]=(snp[endsnp[i]])
				
								
			}else{
				vec[i]=NA
				len[i]=NA
				end[i]=(snp[endsnp[i]])
				
			}
	print(i/n)
}


	creeper = data.frame(cbind(vec, vec2, vec3, len, snp, end))
	creeper$CHROM = rep(k, nrow(creeper))
	names(creeper)=c("pMAS","pFAS","FSTP", "num_snps", "start", "end","chrom")
	creep.all = rbind(creep.all, creeper)


}
	
save.image(file="RAWOUT1KbtotalFST.RData")
creep.all = creep.all[which(complete.cases(creep.all$num_snps)),]


#Trim creeper windos based on their number of SNPs 
len = round(max(creep.all$num_snps)/10) 
creep.all$bin <- as.numeric(cut2(creep.all$num_snps, g=len))
#var.plot = aggregate(creep.all$pMAS, by = list(creep.all$bin), function(x) sd(x, na.rm=T))
#plot(var.plot,pch=16)
creeper = creep.all
creeper=creeper[(creeper$bin > 3),]





high = creeper[which(creeper$FSTP > 5  ),] #or signi/1000
high$group2 = rep("NA", nrow(high))

for(k in unique(high$chrom)){
	ir <- IRanges(high$start[which(high$chrom==k)], high$end[which(high$chrom==k)]) 
	high$group2[which(high$chrom==k)] <- subjectHits(findOverlaps(ir, reduce(ir)))
}


rang=aggregate(start~chrom + group2,data = high, min)
rang.1=aggregate(end~chrom + group2, data = high, function(x) max(x))
rang$end = rang.1$end




#made rang
write.list(rang, file="ClusteredHighSNPsCreeper1kTotal")





# Genes withing high FST regions ----------------------------------------
#I'll pull geneswithin "rang"

hi.Fst = c()
chrom = intersect(rang$chrom, All.Data$CHROM)
for(i in chrom){
	win.temp = rang[rang$chrom==i,]
	deg.temp = All.Data[which(as.character(All.Data$CHROM)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$start)), ">=") 
	blah1=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$end)), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	hi.Fst = rbind(temp,hi.Fst)
	print(i)
}

regions = aggregate(GPOS~GRP+group2, data=hi.Fst, function(x) min(x))
regions$end = aggregate(GPOS~GRP+group2, data=hi.Fst, function(x) max(x))$GPOS



# Genes? ------------------------------------------------
gff = gff[which(gff$type=="gene"),]
gff$chrom = paste("Group", gff$chrom , sep="") 


hi.genes = c()
chrom = intersect(regions$GRP, gff$chrom)
for(i in chrom){
	win.temp = regions[regions$GRP==i,]
	deg.temp = gff[which(as.character(gff$chrom)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$start), as.numeric(as.character(win.temp$end)), "<=") 
	blah1=outer(as.numeric(deg.temp$end), as.numeric(as.character(win.temp$GPOS)), ">=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	hi.genes = rbind(temp,hi.genes)
	print(i)
}

high.genes = unique(hi.genes$GB)

out = hi.genes
out$chrom = gsub("Group", "", out$chrom)
out$start = out$start - 250
out$end = 250 + out$end
write.list(out[c(1,3,4)], file="high.bed")

#write.list(regions[c(1,3,4)], file="VERYhigh.bed")

#Characterizing SNPs -------------------------------------


charac = read.table(file="/media/data1/forty3/drone/git/data/SelectedSitesNSYNSYN.txt",header=T) #output from SNPEFF
write.list(aggregate(Effect~GB, data=charac, table), file="NSYNSYMMARY")
summ = read.table(file="NSYNSYMMARY",header=F)
names(summ) = c("GB","NSYN","SSA","SSR","STOPG","STOPL","SYN")

#					indels = read.table(file="/media/data1/forty3/drone/vcf_drone/Indels.txt",header=T)
#					write.list(aggregate(Effect~GB, data=indels, table), file="indelsSYMMARY")
#					in.summ = read.table(file="indelsSYMMARY",header=F)
#					names(in.summ) = c("GB","CCDEL","CCINS","DEL","INS","FS","SSA","SSD","SSR","STARTL","STOPG","STOPL")
#					in.sum[in.sum$GB %in% summ$GB,]

test = merge(All.Data, charac, by = "SNP")

test = (test[which(-log10(test$FSTP)>2.5),])
write.list(aggregate(Effect~GB, data=test, table), file="ONLYNSYNSYMMARY")
summ.nsyn = read.table(file="ONLYNSYNSYMMARY",header=F)
names(summ.nsyn) = c("GB","NSYN","SSA","SSR","STOPG","STOPL","SYN")

#more "function" in high windows? ----------------------------------
	#are there more genic SNPs in high windows?
#compare HYG.snpeff.eff to HYG.high.snpeff.eff

system("grep -c INTERGENIC  HYG.snpeff.eff") #/ system("wc -l HYG.snpeff.eff")
system("grep -c INTERGENIC  HYG.high.snpeff.eff") #/ system("wc -l HYG.high.snpeff.eff")




#Overlap with DEGs? ----------------------------------------
degs = read.table(file="/media/data1/forty3/drone/git/data/boutinDEGs.txt",header=T)
degs[which(degs$GB %in% hi.genes$GB),]
# (3 or 0)


#Overlap with DEPs? ----------------------------------------
deps = read.table(file="/media/data1/forty3/drone/git/data/FosterDeps.txt",header=T)
deps[which(deps$GB %in% hi.genes$GB),]
#none.



# Pi vs FST ----------------------------------------
Pi.Fst = c()
chrom = intersect(regions$GRP, pi$GRP)
for(i in chrom){
	win.temp = regions[regions$GRP==i,]
	deg.temp = pi[which(as.character(pi$GRP)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$BIN_START)+1000, as.numeric(as.character(win.temp$GPOS)), ">=") 
	blah1=outer(as.numeric(deg.temp$BIN_START), as.numeric(as.character(win.temp$end)), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	Pi.Fst = rbind(temp,Pi.Fst)
	print(i)
}


 t.test(Pi.Fst$TajimaD, pi$TajimaD)

        Welch Two Sample t-test

data:  Pi.Fst$TajimaD and pi$TajimaD
t = -3.5112, df = 645.18, p-value = 0.0004771
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.28809005 -0.08143491
sample estimates:
 mean of x  mean of y
0.08766019 0.27242267


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








# TD vs FST ----------------------------------------
	#I USED REGIONS FOR FTP>5 HERE!
ctd = read.table(file="/media/data1/forty3/drone/vcf_drone/C.Tajima.D",header=T, colClasses = c("character", "numeric", "numeric","numeric"))
#ctd = read.table(file="/media/data1/forty3/drone/vcf_drone/S.Tajima.D",header=T, colClasses = c("character", "numeric", "numeric","numeric"))
mtd = read.table(file="/media/data1/forty3/drone/vcf_drone/M.Tajima.D",header=T, colClasses = c("character", "numeric", "numeric","numeric"))
mtd$CHROM = paste("Group", mtd$CHROM, sep=""); mtd$SNP = paste(mtd$CHROM, mtd$BIN_START,sep=":"); names(mtd)[4] = "TDM"
ctd$CHROM = paste("Group", ctd$CHROM, sep=""); ctd$SNP = paste(ctd$CHROM, ctd$BIN_START,sep=":"); names(ctd)[4] = "TDC"
#ctd = ctd[c(4,5)]
#td = merge(mtd, ctd, by = "SNP"); 
td = mtd
td = td[td$N_SNPS>1,]

td.Fst = c()
chrom = intersect(regions$GRP, td$CHROM)
for(i in chrom){
	win.temp = regions[regions$GRP==i,]
	deg.temp = td[which(as.character(td$CHROM)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$BIN_START)+1000, as.numeric(as.character(win.temp$GPOS)), ">=") 
	blah1=outer(as.numeric(deg.temp$BIN_START), as.numeric(as.character(win.temp$GPOS)), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	td.Fst = rbind(temp,td.Fst)
	print(i)
}

head(td.Fst[td.Fst$CHROM=="Group6.2",])

boxplot(td.Fst$TDC, td$TDC)#,notch=T)

wilcox.test(td.Fst$TDC, td$TDC)








td.Fst$TDC, td$TDC

#For C lineage!
boxplot(td.Fst$TDC[td.Fst$CHROM=="Group6.2"], td$TDC)
> wilcox.test(td.Fst$TDC[td.Fst$CHROM=="Group6.2"], td$TDC)

        Wilcoxon rank sum test with continuity correction

data:  td.Fst$TDC[td.Fst$CHROM == "Group6.2"] and td$TDC
W = 574450, p-value = 0.001399
alternative hypothesis: true location shift is not equal to 0



#####
td$grp = rep("all", nrow(td))
td.Fst$grp = rep("high", nrow(td.Fst));td.Fst = td.Fst[c()]



#
#td.NA.res = rbind(td[c(4,6)], td.Fst[c(4,10)])
#






#Compare FST in regions among AMC ---------------

amc.fst = All.Data
for(fil in dir(pattern="*.weir.fst")){
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




#A vs M vs C ---------------









#	
	
# QWD biased genes? ------------------------			
	#high.genes.perm in QWD List

load(file="/media/data1/forty3/brock/expression_data/DWQ_expression.RData")
#for exclusive genes:
Q<-as.character(Grozinger2007QueenGenes1$GB)
W<-as.character(Grozinger2007WorkerGenes1$GB)
g<-as.character(ZayedDroneWorker$GB); D<-g[ZayedDroneWorker$W.D.1=="D" & (ZayedDroneWorker$casteF1_fdr<.05)]
D1<-setdiff(D,W); D1<-setdiff(D1,Q); W1<-setdiff(W,D); Q1<-setdiff(Q,D)	

test = gff ; qwd = rep("N", nrow(test)); 
qwd[test$GB %in% Q1] = "Q" 
qwd[test$GB %in% W1] = "W" 
qwd[test$GB %in% D1] = "D"   	
test$qwd = qwd

qwd = rep("N", nrow(test)); 
qwd[test$GB %in% high.genes] = "sig" 
test$sig = qwd
aggregate(test$sig, by = list(test$qwd), table)
#  Group.1   x.N x.sig
#1       D   930     3
#2       N 13580    42
#3       Q   403     0
#4       W   355     1



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
 #7 14  0	

#Expressed mostly in mature oocytes and decrease expression until 24hr and/or increase expression
	
	
	
	

#Expression in BeeSpace experiments? ------------------------


	

####
#Pre-amble, are significant genes longer than others?
load(file="/media/data1/forty3/R/frames/BeeSpaceQ.RData")
GRlist=BeeSpaceQ

p=matrix(nr=27, nc=2)
comps=c(3:29)
for (i in 3:29){
	#i=3
	d=which(comps==i)
	total=GRlist[c(1,i)] #change middle for all comparisons (from 3-29)
	names(total)=c("GB","comp")
	total=total[!is.na(total$comp),]
		
	hi.list=c(nrow(total[which(total$GB %in% high.genes & total$comp<0.05),]), nrow(total[which(total$GB %in% high.genes & total$comp>0.05),]) )
	else.list=c(nrow(total[which(!(total$GB %in% high.genes & total$comp<0.05)),]), nrow(total[which(!(total$GB %in% high.genes & total$comp>0.05)),]) )
	test = matrix(nrow=2, ncol=2)
	test[1,]=hi.list
	test[2,]=else.list
	test=fisher.test(test)
	
	print(names(GRlist[i]))
	print(test$p.value)
	
}
 x=FDRcontrol(pvec=(p[,1]))
 #None significant 



	
	
	
#GO analysis ------------------------------------



FBGN=read.table(file="/media/data1/forty3/brock/scripts/FBGN.txt",header=T)

GO.getter <- function(df){
	output = FBGN[which(FBGN$GB %in% df$GB),]
	return(output)
	}
	
#RAN DAVID BP_ALL 
#
#Category	Term	Count	%	PValue	Genes	List Total	Pop Hits	Pop Total	Fold Enrichment	Bonferroni	Benjamini	FDR
#GOTERM_BP_ALL	GO:0043279~response to alkaloid	3	5.769230769230769	0.00611611017132975	FBGN0013984, FBGN0025631, FBGN0000473	40	13	4210	24.288461538461537	0.9815707003971634	0.9815707003971634	8.7547212991953
#GOTERM_BP_ALL	GO:0014070~response to organic cyclic substance	3	5.769230769230769	0.00611611017132975	FBGN0013984, FBGN0025631, FBGN0000473	40	13	4210	24.288461538461537	0.9815707003971634	0.9815707003971634	8.7547212991953
#GOTERM_BP_ALL	GO:0008354~germ cell migration	3	5.769230769230769	0.018715216751410337	FBGN0003391, FBGN0005659, FBGN0003507	40	23	4210	13.728260869565217	0.9999954439669752	0.997865513404851	24.5834186972517
#GOTERM_BP_ALL	GO:0007186~G-protein coupled receptor protein signaling pathway	5	9.615384615384617	0.023932321362161214	FBGN0025631, FBGN0037976, FBGN0036742, FBGN0004435, FBGN0004852	40	120	4210	4.385416666666666	0.9999998582785352	0.9947863099299389	30.354594655877488
#GOTERM_BP_ALL	GO:0060249~anatomical structure homeostasis	3	5.769230769230769	0.034826422463473485	FBGN0013984, FBGN0003391, FBGN0003507	40	32	4210	9.8671875	0.9999999999049108	0.9968772818288929	41.102591260415444
#GOTERM_BP_ALL	GO:0000902~cell morphogenesis	7	13.461538461538462	0.06233031424500982	FBGN0013984, FBGN0003391, FBGN0013759, FBGN0003721, FBGN0035101, FBGN0004435, FBGN0011817	40	309	4210	2.3843042071197407	1.0	0.9997704412630757	61.753593653667885
#GOTERM_BP_ALL	GO:0060429~epithelium development	5	9.615384615384617	0.06341072666175532	FBGN0013984, FBGN0003391, FBGN0025631, FBGN0003507, FBGN0011817	40	164	4210	3.2088414634146343	1.0	0.9991813731582373	62.406462194186204
#GOTERM_BP_ALL	GO:0050896~response to stimulus	10	19.230769230769234	0.07583987461759972	FBGN0013984, FBGN0013759, FBGN0025631, FBGN0037976, FBGN0032147, FBGN0004435, FBGN0000473, FBGN0003507, FBGN0011817, FBGN0025808	40	576	4210	1.8272569444444446	1.0	0.9993476329297016	69.20598147031278
#GOTERM_BP_ALL	GO:0042592~homeostatic process	4	7.6923076923076925	0.07879295919100246	FBGN0013984, FBGN0003391, FBGN0025631, FBGN0003507	40	109	4210	3.8623853211009176	1.0	0.998742317241724	70.64322578097165
#GOTERM_BP_ALL	GO:0065008~regulation of biological quality	8	15.384615384615385	0.08499859577840965	FBGN0013984, FBGN0003391, FBGN0013759, FBGN0025631, FBGN0035101, FBGN0036742, FBGN0003507, FBGN0016984	40	417	4210	2.019184652278178	1.0	0.9983800285994078	73.46191860075264
#GOTERM_BP_ALL	GO:0048610~reproductive cellular process	6	11.538461538461538	0.08838085728518927	FBGN0013984, FBGN0003391, FBGN0003721, FBGN0005659, FBGN0003507, FBGN0016984	40	259	4210	2.4382239382239383	1.0	0.9975797808387583	74.88977071965044
#GOTERM_BP_ALL	GO:0032989~cellular component morphogenesis	7	13.461538461538462	0.0971302788484947	FBGN0013984, FBGN0003391, FBGN0013759, FBGN0003721, FBGN0035101, FBGN0004435, FBGN0011817	40	347	4210	2.123198847262248	1.0	0.9976351038968776	78.2578961701996
#GOTERM_BP_ALL	GO:0030866~cortical actin cytoskeleton organization	2	3.8461538461538463	0.0974196072989818	FBGN0025631, FBGN0016984	40	11	4210	19.136363636363637	1.0	0.9961531097580288	78.3617151383581
#GOTERM_BP_ALL	GO:0030865~cortical cytoskeleton organization	2	3.8461538461538463	0.0974196072989818	FBGN0025631, FBGN0016984	40	11	4210	19.136363636363637	1.0	0.9961531097580288	78.3617151383581
#GOTERM_BP_ALL	GO:0048812~neuron projection morphogenesis	5	9.615384615384617	0.09835031017337892	FBGN0013984, FBGN0003391, FBGN0003721, FBGN0035101, FBGN0004435	40	191	4210	2.755235602094241	1.0	0.994396837041965	78.69254833756551
#GOTERM_BP_ALL	GO:0006928~cell motion	5	9.615384615384617	0.09979233627714601	FBGN0013984, FBGN0003391, FBGN0005659, FBGN0004435, FBGN0003507	40	192	4210	2.7408854166666665	1.0	0.9924675575210851	79.19583164418135
#GOTERM_BP_ALL	GO:0031175~neuron projection development	5	9.615384615384617	0.09979233627714601	FBGN0013984, FBGN0003391, FBGN0003721, FBGN0035101, FBGN0004435	40	192	4210	2.7408854166666665	1.0	0.9924675575210851	79.19583164418135
#	
	
	
	
		
	
# Do they make up any known networks? Are they TFs? Are they central? ------------------------			
	#high.genes.perm vs TRN 
trn = read.table(file="/media/data1/forty3/drone/git/data/trn.txt",header=T)
test = trn 
qwd = rep("N", nrow(test)); 
qwd[test$GB %in% high.genes] = "sig" 
test$sig = qwd
aggregate(test$sig, by = list(test$gamma), table)





	
	
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
#3 Hymenoptera 1304    17
#4     Insecta 8614    72

	
	#Most of them are very old indeed
	
	
	
	
	
# Expression? ------------------------		
	#Use Johnson's Data (Johnson_TableS3.xlsx)
exp = read.table(file="/media/data1/forty3/drone/git/data/johnson_expression.txt",header=T)
test = exp 
qwd = rep("N", nrow(test)); 
qwd[test$GB %in% high.genes] = "sig" 
test$sig = qwd	
aggregate(test$sig, by = list(test$Taxonomy), table)		
#1   Arthropod  213     1
#2         Bee   24     0
#3   Conserved 5400    14
#4 Hymenoptera  391     1
#5      Insect  563     0
#6      Orphan  191     0
#7      O-TRGs  807     0




#Not highly expressed (relative to all expression) in all tissues.
boxplot(log10(1+test$brn[test$brn>25])~test$sig[test$brn>25])




	#check insect and hymenoptera and O-TRG
test2 = test[which(test$brn>25),] 
boxplot(log10(1+test$brn[test$Taxonomy=="O-TRGs"])~test$sig[test$Taxonomy=="O-TRGs"])
#They are not over-expressed in any single Nurse tissue (Brain included).


test = test[-c(1,2)]
aggregate(.~sig, data = test, function(x) c(mean(x), sd(x)/sqrt(length(x))))

boxplot(log10(1+test$brn[test$brn>25])~test$sig[test$brn>25])



	
# Are genes enriched in clusters? -----------------------
clust = read.table(file="/media/data1/forty3/drone/git/data/gt_networks.txt",header=T) 
qwd = rep("N", nrow(clust)); 
qwd[clust$GB %in% high.genes] = "sig" 
clust$sig = qwd
aggregate(clust$sig, by = list(clust$cluster), table)

	
qwd = rep("N", nrow(clust)); 
qwd[clust$GB %in% degs$GB] = "sig" 
clust$sig = qwd
aggregate(clust$sig, by = list(clust$cluster), table)
	
	
qwd = rep("N", nrow(clust)); 
qwd[clust$GB %in% deps$GB] = "sig" 
clust$sig = qwd
aggregate(clust$sig, by = list(clust$cluster), table)



hyg = c(high.genes, as.character(deps$GB), as.character(degs$GB))
qwd = rep("N", nrow(clust)); 
qwd[clust$GB %in% hyg] = "sig" 
clust$sig = qwd
aggregate(clust$sig, by = list(clust$cluster), table)




#not enriched.



#Finally, long-term evolution of these genes (aka, gamma) -------------------

	
# Longer-term selection ------------------------			
	#Gamma in hygiene-associated genes? TD? Pi?

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

boxplot(test$gamma~test$qwd)
	#HYG is higher.....
perm.test(test$gamma, length(high.genes), test$gamma[(test$qwd=="HYG")], 10000)		
		# P (hyg gene Gamma) is 0.0000000295818, so it's significant. (with hyg genes)
#So these genes are under selection voer long term.	And likely more genes under purifying selection in rest of genome
				

nrow(test[which(test$qwd=="HYG" & test$gamma>1),])/nrow(test[which(test$qwd=="HYG"),])
	#18% of them under selection.

gamma.res = test




































#under selection in C?

vcftools --vcf HYG.vcf --keep /media/data1/forty3/drone/git/data/s.txt --TajimaD 1000 --maf 0.05  --out S

Cpi = read.table(file="/media/data1/forty3/drone/vcf_drone/S.Tajima.D",header=T, colClasses = c("character", "numeric","numeric","numeric"))
Cpi$GRP = paste("Group", Cpi$CHROM, sep="")



# Pi vs FST ----------------------------------------
CPi.Fst = c()
chrom = intersect(regions$GRP, Cpi$GRP)
for(i in chrom){
	win.temp = regions[regions$GRP==i,]
	deg.temp = pi[which(as.character(pi$GRP)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$BIN_START)+1000, as.numeric(as.character(win.temp$GPOS)), ">=") 
	blah1=outer(as.numeric(deg.temp$BIN_START), as.numeric(as.character(win.temp$end)), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	CPi.Fst = rbind(temp,CPi.Fst)
	print(i)
}


 t.test(CPi.Fst$TajimaD, Cpi$TajimaD)







test = kaks ; qwd = rep("N", nrow(test)); 
qwd[test$GB %in%  high.genes] = "HYG" 
test$qwd = qwd

boxplot(test$kaks~test$qwd)
aggregate(kaks~qwd, data=test, mean )



test = test ; qwd = rep("N", nrow(test)); 
qwd = rep("N", nrow(test)); 
qwd[test$GB %in%  high.genes] = "HYG" 
test$qwd = qwd

NS = test$PiN/test$PiS
NS[is.na(NS)] = 0
NS[!is.finite(NS)] = 0
test$NvS = NS
#test = test[which(test$NvS>0),]
test = test[which(test$NumS + test$NumN > 5),]

aggregate(NvS~qwd, data=test, mean )
summary(aov(test$NvS~test$qwd))

#quick TD



td = read.table("/media/data1/forty3/brock/balancedSNPs/data/C.recode.vcf.ETD",header=F)
td = aggregate(V4~V1, median, data=td)


test = td ; qwd = rep("N", nrow(test)); 




qwd = rep("N", nrow(test)); 
qwd[test$V1 %in%  high.genes] = "HYG" 
test$qwd = qwd
 
aggregate(V4~qwd, data=test, mean )
summary(aov(test$V4~test$qwd))




#
gamm = read.table(file = "/media/data1/forty3/brock/immune/Gene_gamma.txt", header=T)
kaks = read.table(file="/media/data1/forty3/brock/immune/KAKS.out")
names(kaks) =  c("GB", "c.ks", "m.ks", "y.ks", "a.ks", "c.ka", "m.ka", "y.ka", "a.ka" )
kaks$kaks = kaks$m.ka/kaks$m.ks
kaks = kaks[kaks$m.ks>0,]
Pi = read.table(file="/media/data1/forty3/brock/immune/C.snpeff1.eff.Pi")
names(Pi) = c("GB", "NS","NN", "NumS", "NumN", "PiS", "PiN")


test = merge(Pi, kaks, by = "GB")


test = test ; qwd = rep("N", nrow(test)); 
qwd = rep("N", nrow(test)); 
qwd[test$GB %in%  high.genes] = "HYG" 
test$qwd = qwd







test = test[which(test$PiN>0 & test$PiS>0),]

 x = lm(log10(test$kaks)~ log10(test$PiN))
  plot(log10(test$kaks)~ log10(test$PiN), col = as.factor(qwd),pch=19)

test$res = residuals(x)
aggregate(res~qwd, data=test, mean )
summary(aov(test$res~test$qwd))







