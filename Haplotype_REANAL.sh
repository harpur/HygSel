###
# Drone Hap assoc
###



#Calculate association between haplotype (in window size) within a population for a phenotype
	#http://www.nature.com/mp/journal/v14/n8/full/mp200943a.html
	#I've stored these in /hapstats

	
	
	
vcftools --vcf candidates.recode.vcf  --keep /media/data1/forty3/drone/vcf_drone/controlBees.txt --plink 
#Convert scaff to chrom --------------------------
Rscript /media/data1/forty3/drone/git/ScaffMaptoChr.r out.map
mv out.map outSCAFF.map
mv NA.chrom.map out.map
plink --noweb --file out --recode  --out out

plink --noweb --file out   --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --hap-window 12 --hap-assoc --out hap
	
cut -d'|' -f5 hap.qassoc.hap > test
paste hap.qassoc.hap test > hap.qassoc.hap.snp
	
	
vcftools --vcf candidates.recode.vcf  --remove /media/data1/forty3/drone/vcf_drone/controlBees.txt --plink --out S
#Convert scaff to chrom --------------------------
Rscript /media/data1/forty3/drone/git/ScaffMaptoChr.r S.map
mv S.map SSCAFF.map
mv S.map S.map
plink --noweb --file out --recode  --out S

plink --noweb --file S   --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --hap-window 12 --hap-freq --out S	
	
	
	
P = 0.05
#i = 5
hap = read.table(file="hap.qassoc.hap.snp", header=F, skip = 1)
hap.df = aggregate(hap$V7, by = list(hap$V1, hap$V9), function(x) length(x[which(as.numeric(x)<P)])/length(x))
hap.df$len = aggregate(hap$V7, by = list(hap$V1, hap$V9), function(x) length(x[!is.na(as.numeric(x))]))$x
hap.df$sig.len = aggregate(hap$V7, by = list(hap$V1, hap$V9), function(x) length(x[which(as.numeric(x)<P)]))$x
hap.df$WIN = gsub("WIN","",hap.df$Group.1)
hap.df = hap.df[order(as.numeric(hap.df$WIN)),]



write.table(hap.df, file=paste("haplenC",Pvalue,".summary",sep=""), append=T, col.names=F, row.names=F, quote=F)
print(i)



haps.C = read.table(file=paste("haplenC",Pvalue,".summary",sep=""), header=F)
haps.C$chr = gsub("[.].*","",haps.C$V2) 
 
 

run.haps = c()
chr = c()
for(i in c(1:chrs)){
	x = runmed(haps.C$V3[haps.C$chr==i],k=win.size)
	ch = rep(i, length(x))
	chr = c(chr, ch)
	run.haps = c(x, run.haps)
	
}
 
haps.C$runmed = run.haps
 