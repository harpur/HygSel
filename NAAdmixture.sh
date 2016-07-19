###
# ADMIXTURE for NA bees in selected areas versus control
###




i = "ALL.recode.vcf" 
C.vcf = Read.VCF() 

system("vcftools --vcf /media/data1/forty3/brock/align/N.raw.vcf --bed VERYhigh.bed --recode --out NA")
i = "NA.recode.vcf" 


1. Create tabix index for merging
bgzip /media/data1/forty3/brock/align/AllBees_SNPs.raw.vcf
tabix -p vcf /media/data1/forty3/brock/align/AllBees_SNPs.raw.vcf.gz
bgzip /media/data1/forty3/brock/align/N.raw.vcf.gz
tabix -p vcf /media/data1/forty3/brock/align/N.raw.vcf.gz




2. Merge VCFs
vcf-merge /media/data1/forty3/brock/align/N.raw.vcf.gz /media/data1/forty3/brock/align/AllBees_SNPs.raw.vcf.gz | bgzip -c > NABEES.vcf.gz



3. admixture with all sites
vcftools --vcf NABEES.vcf  --plink  --out NAall --thin 0.2 --remove /media/data1/forty3/drone/git/data/o.txt 
plink --file NAall --noweb  --thin 0.2 --maf 0.05 --make-bed 

for K in 3 ; \
do /home/brock/admixture/admixture  --cv=20 plink.bed $K -j7 | tee log${K}.out; done
mv plink.3.Q plink.3.Q.NAall



4. admixture with high sites
vcftools --vcf NABEES.vcf  --bed VERYhigh.bed --plink  --out NAhigh --remove /media/data1/forty3/drone/git/data/o.txt --remove /media/data1/forty3/drone/git/data/s.txt
plink --file NAhigh  --noweb   --maf 0.05 --make-bed 

for K in 2 ; \
do /home/brock/admixture/admixture  --cv=20 plink.bed $K -j7 | tee log${K}.out; done
mv plink.2.Q plink.2.Q.NAhigh




5. admixture with high sites
vcftools --vcf NABEES.vcf  --bed VERYhigh.bed --plink  --out NAhigh3 --remove /media/data1/forty3/drone/git/data/o.txt
plink --file NAhigh  --noweb   --maf 0.05 --make-bed 

for K in 3 ; \
do /home/brock/admixture/admixture  --cv=20 plink.bed $K -j7 | tee log${K}.out; done
mv plink.3.Q plink.3.Q.NAhigh

