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
pfst.fas = read.table(file="/media/data1/forty3/drone/FST/pFST/pFST.FAS.out",header=F, colClasses = c("character", "numeric", "numeric"))
pfst.mas = read.table(file="/media/data1/forty3/drone/FST/pFST/pFST.MAS.out",header=F, colClasses = c("character", "numeric", "numeric"))
fst = read.table(file="/media/data1/forty3/drone/FST/pFST/wcFST.out",header=F, colClasses = c("character", "numeric", "numeric", "numeric"))
pi = read.table(file="/media/data1/forty3/drone/vcf_drone/sel.windowed.pi",header=T, colClasses = c("character", "numeric", "numeric","numeric","numeric"))
gff = read.table("/media/data1/forty3/R/GFF3_2",header=T, colClasses = c("character", "character", "numeric","numeric","character","character","character","character"))


# Mung Data frames ------------------------------------------
	
#create SNP ID
pfst.mas$SNP = paste(pfst.mas$V1, pfst.mas$V2, sep="_")
pfst.fas$SNP = paste(pfst.fas$V1, pfst.fas$V2, sep="_")
fst$SNP = paste(fst$V1, fst$V2, sep="_")
pi$SNP = paste(pi$CHROM, pi$BIN_START, sep="_")

#create scaff ID
fst$V1 = paste("Group", fst$V1, sep="")
pfst.mas$V1 = paste("Group", pfst.mas$V1, sep="")
pfst.fas$V1 = paste("Group", pfst.fas$V1, sep="")
pi$GRP = paste("Group", pi$CHROM, sep="")

#merge
pfst.mas= pfst.mas[c(3,4)];names(pfst.mas)[1] = "pMAS"
pfst.fas = pfst.fas[c(3,4)];names(pfst.fas)[1] = "pFAS"


fst = data.table(fst)
pfst.mas = data.table(pfst.mas)
pfst.fas = data.table(pfst.fas)
fst = merge(fst, pfst.mas, by="SNP")
fst = merge(fst, pfst.fas, by="SNP")
fst = data.frame(fst)

fst = fst[c(1,2,3,6,7,8)]; names(fst) = c("SNP","GRP","GPOS","FST","pMAS", "pFAS")


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

#identify high sites 
test$pMAS = -1*log10(test$pMAS)
test$pFAS = -1*log10(test$pFAS)

#Identify clusters of high FST regions
	#used creeping window approach


window.size = 5000 #size of window to be crept and gaps to skip, in BP 
bin.size = 5 # minimum number of SNPs in a window (bin size, see plot)
test = test[with(test, order(POS)),]

creep.all = c()
for(k in unique(test$CHROM)){
	
	#Gather SNPs and FST as numeric variables
	snp = as.numeric(test$POS[test$CHROM==k])
	fst = as.numeric(test$pMAS[test$CHROM==k]) 
	fst2 = as.numeric(test$pFAS[test$CHROM==k]) 
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
						
				len[i]=snp_len
				end[i]=(snp[endsnp[i]])
				
								
			}else{
				vec[i]=NA
				len[i]=NA
				end[i]=(snp[endsnp[i]])
				
			}
	print(i/n)
}


	creeper = data.frame(cbind(vec, vec2, len, snp, end))
	creeper$CHROM = rep(k, nrow(creeper))
	names(creeper)=c("pMAS","pFAS","num_snps", "start", "end","chrom")
	creep.all = rbind(creep.all, creeper)


}
	
save.image(file="RAWOUT.RData")
creep.all = creep.all[which(complete.cases(creep.all$num_snps)),]




len = round(max(creep.all$num_snps)/10) 
creep.all$bin <- as.numeric(cut2(creep.all$num_snps, g=len))
#var.plot = aggregate(creep.all$Sp, by = list(creep.all$bin), function(x) sd(x, na.rm=T))
#plot(var.plot,pch=16)
creeper = creep.all
creeper=creeper[(creeper$bin >5),]





high = creeper[which(creeper$pMAS >2 & creeper$pFAS >2  ),] #or signi/1000
high$group2 = rep("NA", nrow(high))

for(k in unique(high$chrom)){
	ir <- IRanges(high$start[which(high$chrom==k)], high$end[which(high$chrom==k)]) 
	high$group2[which(high$chrom==k)] <- subjectHits(findOverlaps(ir, reduce(ir)))
}


rang=aggregate(start~chrom + group2,data = high, min)
rang.1=aggregate(end~chrom + group2, data = high, function(x) max(x))
rang$end = rang.1$end




#made rang
write.list(rang, file="ClusteredHighSNPsCreeper5kb")





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










