
	
# For high SNPs -------------------------
library("FactoMineR")

source("/media/data1/forty3/brock/scripts/VCFFunctions.r")
source("/media/data1/forty3/drone/git/GenomeR/")


pheno = read.table(file="/media/data1/forty3/drone/git/data/DroneSamps.txt")
i = "highPruned.recode.vcf" 
vcf = Read.VCF()	
 #change later
i = "ShighPruned.recode.vcf" 
vcf.S = Read.VCF()	


#Identify major allele in Selected pop ----------------------
s.allele  = sapply(vcf.S,function(x)  substr(unlist(x)[10:39],1 ,3 )) 
s.allele = s.allele[,-1]
alt.freq = apply(s.allele, 2, function(x) (2*length(grep("1/1",x)) + length(grep("0/1",x)))/(2*length(x))) #freq of alt allele
alt.freq[alt.freq>=0.5] = 1
alt.freq[alt.freq<0.5] = 0 #now highlights the S major
#alt.freq = paste(alt.freq, alt.freq, sep = "/")


#Extract Alleles fro each Control Bee ---------------------
	#this is, for now, hard-coded
samps = unlist(vcf[1])[10:20]
pos = sapply(vcf,function(x) return((x[2])))
pos = pos[-1]
chr = sapply(vcf,function(x) return((x[1])))[-1]
snp = paste(chr, pos, sep="_")
allele  = sapply(vcf,function(x)  substr(unlist(x)[10:20],1 ,3 )) 
allele = allele[,-1]

allele1 = substr(allele,1,1)
allele2 = substr(allele,3,3)
majeS  = matrix(alt.freq , nrow(allele), nc = ncol(allele), byrow = T)
allele1 = allele1 == majeS
allele2 = allele2 == majeS

add.allele = allele1 + allele2 #table of additive alleleic effects
dom.allele = add.allele
dom.allele[dom.allele=="2"] = 1


# Dominant Effect Analysis ---------------------------
allele = data.frame(dom.allele)
names(allele) = paste("SNP", snp, sep="_")
allele$V2 = samps
allele = merge(allele, pheno, by = "V2")
num = apply(allele, 2, function(x) length(unique(x)))
allele = allele[,which(num>1)]


ps = c()
for(i in 2:(ncol(allele)-2)){
	#afit <- lm(allele$V4~as.factor(unlist(allele[i])))
	ps=c(ps, summary(aov(allele$V4~as.factor(unlist(allele[i]))))[[1]][["Pr(>F)"]][[1]])
}




x = data.frame(SNP  = names(allele)[2:(ncol(allele)-2)])
x$p = ps
x$SNP = gsub("SNP_","", x$SNP)
x$grp = paste("Group",gsub("_.*","",x$SNP),sep="")
x$POS = gsub(".*_","",x$SNP)


x[x$p<0.05,]
#             SNP           p        grp     POS
#42 13.12_1978613 0.020662138 Group13.12 1978613
#59   16.6_343197 0.020662138  Group16.6  343197
#82    4.5_239941 0.028479563   Group4.5  239941
#92    6.2_522941 0.006333605   Group6.2  522941





















#####################################################################################


#Identify major allele in Selected pop ----------------------
s.allele  = sapply(vcf.S,function(x)  substr(unlist(x)[10:39],1 ,3 )) 
s.allele = s.allele[,-1]
alt.freq = apply(s.allele, 2, function(x) (2*length(grep("1/1",x)) + length(grep("0/1",x)))/(2*length(x))) #freq of alt allele
alt.freq[alt.freq>=0.5] = 1
alt.freq[alt.freq<0.5] = 0 #now highlights the S major
#alt.freq = paste(alt.freq, alt.freq, sep = "/")


#Extract Alleles fro each Control Bee ---------------------
	#this is, for now, hard-coded
samps = unlist(vcf.S[1])[10:39]
pos = sapply(vcf.S,function(x) return((x[2])))
pos = pos[-1]
chr = sapply(vcf.S,function(x) return((x[1])))[-1]
snp = paste(chr, pos, sep="_")

s.allele = data.frame(s.allele); names(s.allele)=paste("SNP", snp,sep="")
s.allele$V2 = samps
s.allele = merge(s.allele, pheno, by = "V2")



boxplot(s.allele$V4~s.allele$SNP13.12_1978613)









#Identify major allele in Selected pop ----------------------
s.allele  = sapply(vcf.S,function(x)  substr(unlist(x)[10:39],1 ,3 )) 
s.allele = s.allele[,-1]
alt.freq = apply(s.allele, 2, function(x) (2*length(grep("1/1",x)) + length(grep("0/1",x)))/(2*length(x))) #freq of alt allele
alt.freq[alt.freq>=0.5] = 1
alt.freq[alt.freq<0.5] = 0 #now highlights the S major
#alt.freq = paste(alt.freq, alt.freq, sep = "/")













