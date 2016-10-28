pheno = read.table(file="/media/data1/forty3/drone/git/data/DroneSamps.txt")
system("vcftools --vcf Drone.Hap.recode.vcf --bed high.bed --keep controlBeesF.txt --recode --out HYGHIGHHAPS --max-alleles 2" ) #change later
i = "HYGHIGHHAPS.recode.vcf" 
vcf = Read.VCF()	



vcftools --vcf Drone.Hap.recode.vcf --bed high.bed --keep controlBeesF.txt --recode --plink --max-alleles 2

Rscript /media/data1/forty3/drone/git/ScaffMaptoChr.r out.map

plink --noweb --file out --recode  --out AllSamplere  

plink --noweb --file AllSamplere --indep-pairwise 1000 10 0.1 #convert the file

sed -i -e 's/:/\t/g' plink.prune.in

vcftools --vcf Drone.Hap.recode.vcf --positions plink.prune.in --keep controlBeesF.txt --recode --max-alleles 2 --out highPruned
vcftools --vcf Drone.Hap.recode.vcf --positions plink.prune.in --remove controlBeesF.txt --recode --max-alleles 2 --out ShighPruned


#going to try PERMUTATION too -----------------------

vcftools --vcf highPruned.recode.vcf --plink --max-alleles 2 

plink --noweb --file out --pheno DronePhenoHB.txt --assoc --perm --qt-means

	

