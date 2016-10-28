###
# Drone Hap assoc
###



#Calculate association between haplotype (in window size) within a population for a phenotype
	#http://www.nature.com/mp/journal/v14/n8/full/mp200943a.html
	#I've stored these in /hapstats

for K in  16 2 3 4 5 6 7 8 9 10 11 12 13 14 15 1; \
do
	vcftools --vcf gwas_$K.phased.vcf --plink --out gwas_$K.phased
	plink --noweb --file gwas_$K.phased  --remove controlBees.txt --pheno DronePhenoHB.txt --hap-window 10 --hap-assoc --out hap_$K; done


#run control bees?
for K in  16 2 3 4 5 6 7 8 9 10 11 12 13 14 15 1; \
do
	vcftools --vcf gwas_$K.phased.vcf --plink --out gwas_$K.phased
	plink --noweb --file gwas_$K.phased  --keep controlBees.txt --pheno DronePhenoHB.txt --hap-window 10 --hap-assoc --out hapC_$K; done
	

for K in  16 2 3 4 5 6 7 8 9 10 11 12 13 14 15 1; \
do	
cut -d'|' -f5 hap_$K.qassoc.hap > test
paste hap_$K.qassoc.hap test > hap_$K.qassoc.hap.snp;done

 
 
for K in  16 2 3 4 5 6 7 8 9 10 11 12 13 14 15 1; \
do
	plink --noweb --file gwas_$K.phased  --remove controlBees.txt --pheno DronePhenoHB.txt --assoc --qq-plot --qt-means --out hap_$K; done
	 
  
 
 
#Assemble data ------------------------------------
Rscript <HERE> 16 /media/data1/forty3/drone/PLINK/hapflk-1.3.0/chr/


 
 
 
 
 
 
 
 
 
 