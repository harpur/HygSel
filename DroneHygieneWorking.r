#####
# Hygiene SNP Analyses
#####



#Creation of VCF Files can be found in VCFCreation_DroneSelection.txt
	#The final VCF file is /media/data1/forty3/drone/vcf_drone/DroneSelectionFinal.recode.vcf


	
	
#found High FST SNPs in DroneFST.r	
	
	
	
######	
#GWAS WORK (In Brief):
		#I say in Brief, because it is spread amongst DroneEIG.txt and DronePLINK.txt	

# From the final VCF file I will be looking for associations
	#http://scholar.google.ca/scholar?q=artificial+selection+GWAS+Fst&btnG=&hl=en&as_sdt=0%2C5
	#http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0040736#pone-0040736-g008
	#http://reports-archive.adm.cs.cmu.edu/anon/ml2011/CMU-ML-12-104.pdf
	#https://www.cs.cmu.edu/~ggordon/10725-F12/scribes/10725_Lecture2.pdf
	#http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1002316#s4
	#http://www.biomedcentral.com/1471-2164/13/48/

#Vaysse et al. did a permutation to control for pop stratificiation. Going to try this with my high FST SNPs 	
	#http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1002316
	#Doesn't work BAH!!!
	
	
#I want to try EIG with my high FST SNPs to see if it can fix the stratification issue. 
	#I took my high FST SNPs and a random set of SNPs then ran PLINK linear with COVARs from EIG
	#in EIG (from DroneEIG.txt):
		covar=read.table(file="DRONE.pca.evec",skip=1)
		covar=covar[,c(1:3)];covar=cbind(covar[,1], covar)
		write.table(covar, file="/media/data1/forty3/drone/vcf_drone/HBPCA.txt",quote=F,col.names=F,row.names=F)
		
		perl bahexample.perl
		#Saved the two significant PCAs as HBPCA and will run PLINK linear with them as covariates.
		#Supplemental Plot:
		pop=read.table(file="pops.txt")
		ind=read.table(file="Drone.ind",header=F)
		pc=read.table(file="DRONE.pca",header=F,skip=11)
		col=pop$V1
		col=as.factor(pop$V1)
		pdf("Drone_PCA.pdf")
		plot(pc$V1,pc$V2,pch=19, col=col, xlab="PCA1",ylab="PCA2")
		text(pc$V1, pc$V2,ind$V1)
		dev.off()


	
plink --file DroneSelection --mpheno 1 --pheno DronePhenoHB.txt --extract GenicHighestFSTSNPsRAND.snp --noweb --allow-no-sex --linear  --covar HBPCA.txt --adjust 
	R
	source("/media/data1/forty3/brock/scripts/qq.r")
	snps=read.table(file="plink.assoc.linear.adjusted",header=T)
	x11();qq(snps$UNADJ)
plink --file DroneSelection --mpheno 1 --pheno DronePhenoHB.txt --extract GenicHighestFSTSNPsRAND.snp --noweb --allow-no-sex --linear --adjust --out noCoVar 
	UNsnps=read.table(file="noCoVar.assoc.linear.adjusted",header=T)
	x11();qq(UNsnps$UNADJ)
	
	
	#While this works, I get boned on the FDR corrections. 
		#Here, testing is sample size fixes that:
for K in {1..100}
do
	shuf -n 10 GenicHighestFSTSNPs.snp > output$K
	plink --file DroneSelection --mpheno 1 --pheno DronePhenoHB.txt --extract output$K --noweb --allow-no-sex --linear  --covar HBPCA.txt --adjust --out RAND$K
done
	
	
	
plink --file DroneSelection --mpheno 1 --pheno DronePhenoHB.txt --extract GenicHighestFSTSNPsRAND.snp --noweb --allow-no-sex --assoc  --qt-means 








######	
#Runs of Homozygosity in selected lines?
	
	
plink --file DroneSelection --homozyg --keep SelnBees.txt --noweb --homozyg-group --homozyg-snp 50


	#Homozygous segment criteria:
	#  length (kb)       = 1000
	#  # SNPs (N)        = 100
	#  density (kb/SNP)  = 50
	#  largest gap (kb)  = 1000
	#27 of 27 individuals
	#








	



#I am going to find genes near to the SNPs of interest
#This picks up from DroneFST.r

























load(file="/media/data1/forty3/drone/R/SelectionScansPop12.RData")

#a list of genes that are of high FST and associated with hygiene.
candidates = GenicFST[which(GenicFST$SNP %in% test$SNP),]
	#67 unique genes
	#none of them pass Benjamini, but the take-away is that they are all involved in neuronal growth and guidance
	
	



	
candsnp = read.table(file="HLCandidate.assoc.linear.adjusted",header=T)
candsnp = candsnp[which(candsnp$FDR_BY>0.05),]
	
	
candidates = GenicFST[which(GenicFST$SNP %in% test$SNP),]
candidates = candidates[which(candidates$CHROM %in% candsnp$CHR),]	
	
	
	
	
	
	
	
	
	
	
	
	
	