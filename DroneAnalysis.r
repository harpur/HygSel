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
pi = read.table(file="/media/data1/forty3/drone/vcf_drone/sel.windowed.pi",header=T, colClasses = c("character", "numeric", "numeric","numeric","numeric"))
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
All.Data = fst
test = fst

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
	fst3 = as.numeric(test$FST[test$CHROM==k]) 	
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
creeper=creeper[(creeper$bin >3),]





high = creeper[which(creeper$pMAS >2.5 & creeper$pFAS >2.5  ),] #or signi/1000
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
	blah=outer(as.numeric(deg.temp$start-500), as.numeric(as.character(win.temp$end)), "<=") 
	blah1=outer(as.numeric(deg.temp$end+500), as.numeric(as.character(win.temp$GPOS)), ">=") 
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
chrom = intersect(rang$chrom, pi$CHROM)
for(i in chrom){
	win.temp = rang[rang$chrom==i,]
	deg.temp = pi[which(as.character(pi$CHROM)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$end), as.numeric(as.character(win.temp$start)), ">=") 
	blah1=outer(as.numeric(deg.temp$start), as.numeric(as.character(win.temp$end)), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	Pi.Fst = rbind(temp,Pi.Fst)
	print(i)
}



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


#Judging by this, it seems Selected bees at these regions have become "more C" and less "M" than control. 


x11();boxplot(AMC.high.fst$C.CON, amc.fst$C.SEL, notch = T)





#Quick plot to show purging effect:

x11();hist(AMC.high.fst$C.CON)#, amc.fst$C.SEL, notch = T)
















	
	
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
#2     Apoidea  215     0
#3 Hymenoptera 1313     8
#4     Insecta 8657    29
	
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




#not enriched...:(


