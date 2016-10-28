###
# Assess correlations
###


# Uses the regions of High FST to see if SNPs within them are, on average, associated with phenotype	
	
	
	# Takes in SNPs found within "rang" df 	-> "highFSTSNPs.recode.vcf"
	# creates a dataframe of allele frequencies for all SNPs (df)
	# correlates phenotypic value to all SNPs within a region, averages the correlation coefficient (ignores significance_
	# outputs average correlation coefficient for all groups (df.cor)
	# then, takes in random SNPs (not in high FST windows) and does the same analysis, isolating clusters of SNPs and associating them to phenotype, averaging. Does this thousands of times to get an expected sitribution (pnorm.cor)
	
	
# For high SNPs -------------------------
library("FactoMineR")

source("/media/data1/forty3/brock/scripts/VCFFunctions.r")



pheno = read.table(file="/media/data1/forty3/drone/git/data/DroneSamps.txt")
system("vcftools --vcf Drone.Hap.recode.vcf --bed high.bed --keep controlBeesF.txt --recode --out HYGHIGHHAPS --max-alleles 2" ) #change later
i = "HYGHIGHHAPS.recode.vcf" 
vcf = Read.VCF()	



system("vcftools --vcf Drone.Hap.recode.vcf --bed high.bed --remove controlBeesF.txt --recode --out HYGHIGHHAPS --max-alleles 2" ) #change later
i = "HYGHIGHHAPS.recode.vcf" 
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
tt = c()
for(i in 2:(ncol(allele)-2)){
	#afit <- lm(allele$V4~as.factor(unlist(allele[i])))
	ps=c(ps, summary(aov(allele$V4~as.factor(unlist(allele[i]))))[[1]][["Pr(>F)"]][[1]])
	tt=c(tt, TukeyHSD(aov(allele$V4~as.factor(unlist(allele[i]))))[[1]][1])
}




x = data.frame(SNP  = names(allele)[2:(ncol(allele)-2)])
x$p = ps
x$Tuk = tt
x$SNP = gsub("SNP_","", x$SNP)
x$grp = paste("Group",gsub("_.*","",x$SNP),sep="")
x$POS = gsub(".*_","",x$SNP)



#load "regions"

Pi.Fst = c()
chrom = intersect(regions$GRP, x$grp)
for(i in chrom){
	win.temp = regions[regions$GRP==i,]
	deg.temp = x[which(as.character(x$grp)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$GPOS)), ">=") 
	blah1=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$end)), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	Pi.Fst = rbind(temp,Pi.Fst)
	print(i)
}

#Pi.Fst = Pi.Fst[Pi.Fst$Tuk>0,]
test = aggregate(Pi.Fst$p, by = list(Pi.Fst$GRP,Pi.Fst$group2 ), function(x) sd(round(x,2)))
test$meanP = aggregate(Pi.Fst$p, by = list(Pi.Fst$GRP,Pi.Fst$group2 ), mean)$x
test$minP = aggregate(Pi.Fst$p, by = list(Pi.Fst$GRP,Pi.Fst$group2 ), min)$x
test$maxP = aggregate(Pi.Fst$p, by = list(Pi.Fst$GRP,Pi.Fst$group2 ), max)$x
test$N = aggregate(Pi.Fst$p, by = list(Pi.Fst$GRP,Pi.Fst$group2 ), length)$x
test$T = aggregate(Pi.Fst$Tuk, by = list(Pi.Fst$GRP,Pi.Fst$group2 ), mean)$x






# Permutation -------------------------------------------------
	#randomly draw N snps in a group, associate them to hygiene in control population

system("vcftools --vcf Drone.Hap.recode.vcf --exclude-bed high.bed --keep controlBeesF.txt --recode --out HYGHIGHHAPS --max-alleles 2" ) #change later
i = "HYGHIGHHAPS.recode.vcf" 
vcf = Read.VCF()	


system("vcftools --vcf Drone.Hap.recode.vcf --exclude-bed high.bed --remove controlBeesF.txt --recode --out HYGHIGHHAPS --max-alleles 2" ) #change later
i = "HYGHIGHHAPS.recode.vcf" 
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

#ended here, grab random snps in groups and associate.!


ps = c()
for(i in 2:(ncol(allele)-2)){
	#afit <- lm(allele$V4~as.factor(unlist(allele[i])))
	ps=c(ps, summary(aov(allele$V4~as.factor(unlist(allele[i]))))[[1]][["Pr(>F)"]][[1]])
}



#Draw N snps from the Test distribution and get a hist of the "max" and "mean" in that draw --------------------------------

N=10000
perm.df = c()
for(i in unique(test$N)){
	len = c()
	mins = c()
	means = c()
	for(k in 1:N){
		ind = sample(1:length(ps), 1)
		mins = c(mins, min(ps[c(ind:(ind+(i-1)))]))
		means = c(means, mean(ps[c(ind:(ind+(i-1)))]))
		len = rep(i,length(1:N) ) 
	}
	temp = cbind(mins, means, len)
	perm.df=(rbind(perm.df, temp))
	
}


perm.df = data.frame(perm.df)
perm.means = aggregate(perm.df$means, by = list(perm.df$len),mean, na.rm=T)
perm.means$sd = aggregate(perm.df$means, by = list(perm.df$len),sd, na.rm=T)$x
names(perm.means) = c("N","mn", "sd")


#get pvalues -------------
test = merge(test, perm.means, by = "N")
test$P = apply(test, 1, function(x) pnorm(as.numeric(x[5]), as.numeric(x[8]), as.numeric(x[9])))
test = test[test$N>5,]


sigs = rbind(test[test$P>0.95 & test$minP<0.05,],test[test$P<0.05 & test$minP<0.05,] )
	#there is an interesting effect going on where fixed alleles are not necessarily associated positively in control population.
sigs = test[test$P<0.05 & test$minP<0.05,] 






#Plotting Set up--------------------------

#Get genotypes at these SNPs
sig.snps = Pi.Fst[Pi.Fst$grp %in% sigs$Group.1,]
sig.snps = sig.snps[sig.snps$group2 %in% sigs$Group.2,]
sig.snps = sig.snps[sig.snps$p<0.05,]

#              SNP           p      Tuk       grp    POS       GRP group2   GPOS
#18839  6.2_522891 0.006333605 19.10468  Group6.2 522891  Group6.2     19 522891
#18840  6.2_522892 0.006333605 19.10468  Group6.2 522892  Group6.2     19 522891
#18842  6.2_522904 0.006333605 19.10468  Group6.2 522904  Group6.2     19 522891
#18843  6.2_522941 0.006333605 19.10468  Group6.2 522941  Group6.2     19 522891
#18844  6.2_522955 0.006333605 19.10468  Group6.2 522955  Group6.2     19 522891
#23438  7.6_202402 0.035262890 14.75880  Group7.6 202402  Group7.6      1 201914
#23440  7.6_202614 0.035262890 14.75880  Group7.6 202614  Group7.6      1 201914
#23441  7.6_202631 0.035262890 14.75880  Group7.6 202631  Group7.6      1 201914
#23442  7.6_202644 0.035262890 14.75880  Group7.6 202644  Group7.6      1 201914
#9290  14.5_185769 0.047318543 14.09358 Group14.5 185769 Group14.5      4 185733
#9291  14.5_185800 0.047318543 14.09358 Group14.5 185800 Group14.5      4 185733
#9298  14.5_185881 0.047318543 14.09358 Group14.5 185881 Group14.5      4 185733
#9308  14.5_186865 0.047318543 14.09358 Group14.5 186865 Group14.5      4 185733
#9309  14.5_187098 0.047318543 14.09358 Group14.5 187098 Group14.5      4 185733
#9310  14.5_187139 0.047318543 14.09358 Group14.5 187139 Group14.5      4 185733
#9311  14.5_187142 0.047318543 14.09358 Group14.5 187142 Group14.5      4 185733
#9312  14.5_187219 0.047318543 14.09358 Group14.5 187219 Group14.5      4 185733


#Extract Alleles fro each Control Bee ---------------------
	#this is, for now, hard-coded
samps = unlist(vcf[1])[10:20]
pos = sapply(vcf,function(x) return((x[2])))
pos = pos[-1]
chr = sapply(vcf,function(x) return((x[1])))[-1]
snp = paste(chr, pos, sep="_")
allele  = sapply(vcf,function(x)  substr(unlist(x)[10:20],1 ,3 )) 
allele = allele[,-1]
allele = data.frame(allele)
names(allele) = paste("SNP", snp, sep="_")
ref = sapply(vcf,function(x) return((x[4])))
alt = sapply(vcf,function(x) return((x[5])))

allele = allele[,names(allele) %in% c("SNP_6.2_522904", "SNP_7.6_202614", "SNP_14.5_186865")]
allele$V2 = samps
allele = merge(allele, pheno, by = "V2")

ref = ref[snp %in% c("6.2_522904", "7.6_202614", "14.5_186865")]
alt = alt[snp %in% c("6.2_522904", "7.6_202614", "14.5_186865")]
snp = snp[snp %in% c("6.2_522904", "7.6_202614", "14.5_186865")]

allele$V2 = samps
allele = merge(allele, pheno, by = "V2")
num = apply(allele, 2, function(x) length(unique(x)))
allele = allele[,which(num>1)]
#added alleles by hand


































###

#Going to add a plot of correlation via creeper ----------------
 #change later




###











 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 




#####################################################












