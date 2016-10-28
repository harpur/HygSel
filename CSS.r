###
#CSS with GWAs and FST?
###


#original SP set?

load


vcftools --vcf DroneSelection.vcf --plink --out outDSNP

Rscript /media/data1/forty3/drone/git/ScaffMaptoChr.r outDSNP.map

mv outDSNP.map outDSNP_scaff.map
mv NA.chrom.map outDSNP.map

plink --noweb --file outDSNP --recode  --out AllSamplere  

#plink --noweb --file   AllSamplere   --indep-pairwise 50 5 0.4 --maf 0.05
#plink --noweb \
	#--file AllSamplere   \
	#--extract plink.prune.in \
	#--recode \
	#--out gwasINDEP2


plink --noweb --file AllSamplere  --remove /media/data1/forty3/drone/vcf_drone/controlBees.txt --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --assoc --adjust --qq-plot


R
load("HYGRESULTS_FIGS.RDATA")
#
ass = read.table("plink.qassoc.adjusted", header=T, nrow=100000)
plot(-log10(ass$QQ), -log10(ass$UNADJ))
abline(0,1)

ass = read.table("plink.qassoc", header=T)
plot(-log10(ass$P),col=ass$CHR, main = "SNPs_Trim", pch=19)



#Try the join statist?


test = All.Data
ass$SNP = gsub("[:]","_",ass$SNP)
test = merge(test, ass, by = "SNP")
p = test$P
p[is.na(p)] = 1
test$P = p

#composite index:

test$FST
test$P = -log10(test$P)
test$FSTP = -log10(test$FSTP)



p.rank = rank(test$P, na.last ="keep")
f.rank = rank(test$FSTP, na.last ="keep")


p.rank01 = p.rank/ (max(p.rank,    na.rm = T) + 1)
f.rank01 = f.rank/ (max(f.rank,    na.rm = T) + 1)

p.z = qnorm(p.rank01) #I THINK?
f.z = qnorm(f.rank01)


test.mz <- (p.z + f.z)/2
p.test.mz <- pnorm(test.mz, 0, 1/sqrt(2), lower.tail = F) #  zbar ~ N(0, 1/sqrt(3))

css = -log10(p.test.mz))

test$css = css
test = test[order(test$POS),]
test = test[order(c(test$CHR, test$POS)),]



### add results to the input file
p.test.mz=as.data.frame(p.test.mz); test.mz=as.data.frame(test.mz);
results= cbind(test.mz, p.test.mz, -log10(p.test.mz)); names(results)=c('mz','p_mz','css');
ouput.data=cbind(input.data, results)
ouput.data=as.data.frame(ouput.data)
print("MeanZ, P of meanZ, -log10 of P of meanZ has been calculated and added to the output file as: 'mz', 'p_mz', 'css'")
ouput.data




input.data, tests=c('fst','daf','xpehh'), ihs=F
test.stat=input.data
gt=tests; gt1=length(gt)  ## default

for (i in 1:length(gt)){
test.rank[,gt[i]]    <- rank(test.rank[,gt[i]],       na.last = "keep")


}
for (i in 1:length(gt)){
test.rank[,gt[i]]    <- test.rank[,gt[i]]    / (max(test.rank[,gt[i]],    na.rm = T) + 1)
}

test.z = test.stat
for (i in 1:length(gt)){
test.z[,gt[i]]    <- qnorm(test.rank[,gt[i]])
}
test.mz <- apply(test.z[,gt], 1, mean)
p.test.mz <- pnorm(test.mz, 0, 1/sqrt(gt1), lower.tail = F) #  zbar ~ N(0, 1/sqrt(3))
### add results to the input file
p.test.mz=as.data.frame(p.test.mz); test.mz=as.data.frame(test.mz);
results= cbind(test.mz, p.test.mz, -log10(p.test.mz)); names(results)=c('mz','p_mz','css');
ouput.data=cbind(input.data, results)
ouput.data=as.data.frame(ouput.data)
print("MeanZ, P of meanZ, -log10 of P of meanZ has been calculated and added to the output file as: 'mz', 'p_mz', 'css'")
ouput.data














p.rank = rank(test$P, ties.method = "min")
f.rank = rank(test$FSTP, ties.method = "min")


p.rank01 = p.rank/(length(p.rank) +1)
f.rank01 = f.rank/(length(f.rank) +1)



#convert to Z:

p.z = qnorm(p.rank01) #I THINK?
f.z = qnorm(f.rank01)

mn.z = (p.z + f.z)/2

p = 1-pnorm(mn.z, mean=1, sd = 1)


#http://bmcgenet.biomedcentral.com/articles/10.1186/1471-2156-15-34#Fig1

###
#
###












