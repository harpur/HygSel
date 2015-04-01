#######
# Drone Associatino Workflow
#######
















#####
#Fst Analyses



#VCFTools
	#Calculate FST  all selected VS control (2 vs 1+3=Selpop.txt)
vcftools --maf 0.05 --vcf DroneSelectionFinal.recode.vcf --weir-fst-pop pop2.txt --weir-fst-pop Selpop.txt --out pop2_vs_sel 




#######
# Output High FST Regions:



R

###
#FST Functions
source("/media/data1/forty3/brock/scripts/GOgetter.r")
source("/media/data1/forty3/brock/scripts/VarFunct.r")
source("/media/data1/forty3/brock/scripts/movingavg.r")

#####

load(file="/media/data1/forty3/drone/FST/SelvsCon/SelectedvsControlFST.RData") 
highFstWindow98=quantile(unlist(Fst$x),0.98)
highFstWindow95=quantile(unlist(Fst$x),0.95)
#highFstWindow90=quantile(unlist(Fst$x),0.90)


HighFstPositions98=as.matrix(apply(Fst, 1, function(x) unlist(as.numeric(unlist(x[3]))[which((as.numeric(unlist(x[2]))>=highFstWindow98))])))
HighFstPositions95=as.matrix(apply(Fst, 1, function(x) unlist(as.numeric(unlist(x[3]))[which((as.numeric(unlist(x[2]))>=highFstWindow95))])))
#HighFstPositions90=as.matrix(apply(Fst, 1, function(x) unlist(as.numeric(unlist(x[3]))[which((as.numeric(unlist(x[2]))>=highFstWindow90))])))


dataset=HighFstPositions98
HighSNPs98=c()
chr = apply(dataset, 1,function(x) length(as.numeric(unlist(x))))
chr = names(chr[chr>0])
for (i in chr){
	test=fst23[fst23$CHROM==i,]
	low=as.numeric(unlist(dataset[i,1]))-1000
	high=as.numeric(unlist(dataset[i,1]))+1000
	blah=outer(as.numeric(unlist(low)), as.numeric(unlist(test$POS)), "<=") 
	blah1=outer(as.numeric(unlist(high)), as.numeric(unlist(test$POS)), ">=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	blah=(test[blah[,2],])
	blah=blah[!duplicated(blah),]
	HighSNPs95=rbind(HighSNPs98,blah)
	print(i)
	#x11();plot(as.numeric(unlist(Fst$Pos[Fst$Group.1==i])),as.numeric(unlist(Fst$x[Fst$Group.1==i])), pch=19)
	#abline(v=blah$PO,col="red")
}



#With these SNPs, create "regions list" for association work:
dataset = HighSNPs98
scaff = gsub(":.*","",dataset$SNP)
Unscaff = unique(scaff) # 19 unique
for (i in Unscaff){
	print(i)
	snps = dataset$SNP[which(scaff==i)]
	num=which(unique(scaff)==i)
	sink(file=paste("Candidates98_", "ALL", ".set",sep=""),append=T)
	cat("\n")
	cat(paste("Set",num,sep=""))
	cat("\n")
	write.table(snps,col.names=F,quote=F,row.names=F)
	cat("END")
	cat("\n")
	sink()
}


#Made 2 files Candiates98.set, and Candiates95.set, formatted as below, for PLINK SET tests

#Set1
#1.32:131761
#1.32:131040
#1.32:131037
#1.32:131008
#...
#END

#Also wrote out candiadates for selective association
write.table(HighSNPs98$SNP,file="CandidateSNPs98.snps",col.names=F,quote=F,row.names=F)


save.image(file="DroneSelection.RData")




#Fuck, have to redo this -> do the same thing I did with RPART trees.
###
#Association tests with a  SET TEST
#http://pngu.mgh.harvard.edu/~purcell/plink/anal.shtml#set

#Control bees only:
plink --file DroneSelection --out ControlHighFST --noweb --make-bed --keep controlBees.txt  --extract CandidateSNPs98.snps
plink --bfile ControlHighFST  --recode --out ControlHighFST --noweb 

plink --bfile ControlHighFST --set-test --set Candidates98_ALL.set --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --linear --mperm 1000 --out Candidates98_ALL01 --set-p 0.01 --set-max 5


#Testing one scaffold only.
#for fil in Candidates98_*;
#do plink --bfile ControlHighFST --set-test --set $fil --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/#DronePhenoHB.txt --noweb --linear --mperm 1000 --out $fil --set-p 0.01 --set-max 5; done


#--assoc does not work with -set-test
#plink --file DroneSelection --set-test --set highFST.set --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --allow-no-sex --assoc --mperm 1000 --out controlset --set-p 0.01 --set-max 10 --qt-means --keep controlBees.txt
	#If I sub in "Candidates.set" it stops working.
	#I can't use candidate and "keep" in the same run.

plink --bfile ControlHighFST --set-test --set Candidates98_ALL.set --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --linear --mperm 10000 --out Candidates98_ALLPHENO011000015 --set-p 0.01 --set-max 15

plink --bfile ControlHighFST --set-test --set Candidates98_ALL.set --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --linear --mperm 10000 --out Candidates98_ALLPHENO01100005 --set-p 0.01 --set-max 5

plink --bfile ControlHighFST --set-test --set Candidates98_ALL.set --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --linear --mperm 10000 --out Candidates98_ALLPHENO01100005 --set-p 1 --set-max 1



#expression candidates?
	#control
for i in {1..9}
do
   plink --bfile ControlHighFST --set-test --set Candidates98_ALL.set --mpheno $i --pheno /media/data1/forty3/drone/vcf_drone/HBexpression.txt --noweb --linear --mperm 1000 --out Candidates98_ALLPHENO$i --set-p 0.01 --set-max 5
   
done

	
	   plink --bfile ControlHighFST  --mpheno 7 --pheno /media/data1/forty3/drone/vcf_drone/HBexpression.txt --noweb --assoc --qt-means 	
	#7.12:417504
	
#Selected bees only: 
plink --file DroneSelection --out SelHighFST --noweb --make-bed --remove controlBees.txt
plink --bfile SelHighFST  --recode --out SelHighFST --noweb 

plink --bfile SelHighFST --set-test --set Candidates98_ALL.set --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --linear --mperm 1000 --out Candidates98_ALLSEL --set-p 0.05 --set-max 5 

	#The files Candidates98_*.assoc.linear.set.mperm contain my associated SNPs
	#means for eahc SNP can also be found in Candidates98_ALL*.qassoc.means


plink --bfile SelHighFST --set-test --set Candidates98_ALL.set --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --linear --mperm 10000 --out SELCandidates98_ALLPHENO011000015 --set-p 0.01 --set-max 15

plink --bfileSelHighFST --set-test --set Candidates98_ALL.set --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --linear --mperm 10000 --out SELCandidates98_ALLPHENO01100005 --set-p 0.01 --set-max 5

plink --bfile SelHighFST --set-test --set Candidates98_ALL.set --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --linear --mperm 10000 --out SELCandidates98_ALLPHENO01100005 --set-p 1 --set-max 1
	
	
	
conset=c()	
for(i in 1:9){
CONset=read.table(paste("Candidates98_ALLPHENO", i,".assoc.linear.set.mperm",sep=""),header=T)
CONset$pheno=rep(i, nrow(CONset))
CONset=CONset[CONset$EMP1<0.01,]
conset=rbind(CONset,conset)	
}


selset=c()	
for(i in 1:9){
SELset=read.table(paste("Candidates98_ALLSELPHENO", i,".assoc.linear.set.mperm",sep=""),header=T)
SELset$pheno=rep(i, nrow(SELset))
SELset=SELset[SELset$EMP1<0.01,]
selset=rbind(SELset,selset)	
}
	
	
	
	
	
	
#Select and Control together
plink --file DroneSelection --set-test --set Candidates98_ALL.set --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --linear --mperm 1000 --out AllBees98Candidates --set-p 0.01 --set-max 5

plink --file DroneSelection --set-test --set Candidates98_ALL.set --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --linear --mperm 1000 --out AllBees98Candidates001 --set-p 0.001 --set-max 5








	
	
	
	
	
	
	
##I'd like to try these independant snps with PLINK assoc
#1) Get out independant SNPs within each segment 
plink --file DroneSelection  --out CandidatesControl --indep-pairwise 50 5 0.5 --noweb --keep controlBees.txt   --extract CandidateSNPs98.snps

plink --file DroneSelection  --out CandidatesSelect --indep-pairwise 50 5 0.5 --noweb --remove controlBees.txt   --extract CandidateSNPs98.snps



plink --file DroneSelection --out Candidates --noweb --make-bed  --extract Candidate.prune.in
plink --bfile Candidates  --recode --out Candidates --noweb 

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
#Plotting these data:	
source(DronePlot1_ChromosomeFST.r)
		#outputs FSTSelectTest.Chromosome*_*_updated.pdf
	
	
	
	
	
	
	
###
#Association tests with RPART
#http://pngu.mgh.harvard.edu/~purcell/plink/anal.shtml#set	
	

#As above, get list of independant SNPs in each chromosome
#Run RPART with all those SNPs
	
#1) Get out independant SNPs within each segment 
plink --file DroneSelection  --out CandidatesControl --indep-pairwise 50 5 0.5 --noweb --keep controlBees.txt   --extract CandidateSNPs98.snps

plink --file DroneSelection  --out CandidatesSelect --indep-pairwise 50 5 0.5 --noweb --remove controlBees.txt   --extract CandidateSNPs98.snps



plink --file DroneSelection --out Candidates --noweb --make-bed  --extract Candidate.prune.in
plink --bfile Candidates  --recode --out Candidates --noweb 



#2) Upload these SNPs and run RPART

R
library(randomForest)
library(foreach)
library(rattle)
library(rpart.plot)
library(RColorBrewer)
library(rpart)
#library(doParallel)
#registerDoParallel(cores=detectCores(all.tests=TRUE))
#library(bigrf)
library(doMC)
registerDoMC(15)
	
	
	
#Load SNPs:
map=read.table(file="Candidates.map")
dt = read.table("Candidates.ped", colClasses=c(rep("character",nrow(map)+6)))
pheno=read.table(file="DronePhenoHB.txt",header=F)
pheno=pheno[,-c(2)]
vec=merge(pheno,dt,by="V1")
pop=read.table(file="DroneSelection201Candidates.fam",header=F)
vec=vec[,-c(3:7)]
names(vec)[2]="pheno"
nms=vec$V1
vec=vec[,-1]
phens=vec$pheno
vec=vec[,-1]


c1=vec[,seq(1,ncol(vec),2)]
c2=vec[,seq(2,ncol(vec),2)]

pfun <- function(x, y) paste(x, y, sep = "")
datNew <- vector("list", ncol(c1)/2)
for (i in 1:ncol(c1)) {
    datNew[[i]] <- pfun(c1[[i]], c2[[i]])
}

vec <- as.data.frame(datNew);names(vec)=seq(1:ncol(vec))
vec$pheno=phens
vec$fam=pop$V3
names(vec)=paste("SNP",names(vec),sep="")





#Random Forest:
set.seed(213)
RFfit <- randomForest(SNPpheno ~ .,  data=vec, importance=TRUE, ntree=200, mtry=100)
 pdf(file="ImportanceSNPTree.pdf")
 varImpPlot(RFfit)
dev.off()
#Most important variables, based on %incMSE
map[c(1637,3622,3476,1623,2744,259,1628,2676,2757,3428),]
#     V1            V2 V3       V4
#1637  7    7.13:20637  0  4455504
#3622 15 15.19:2076143  0  8086642
#3476 14   14.5:207937  0  1737995
#1623  7     7.13:7421  0  4442288
#2744 11 11.18:1311002  0 10215622
#259   3   3.8:1748491  0  5715317
#1628  7    7.13:13474  0  4448341
#2676 11 11.18:1271860  0 10176480
#2757 11 11.18:1318408  0 10223028
#3428 14   14.5:187139  0  1717197
#The importance measures show how much MSE or Impurity increase when that variable is randomly permuted.  If you randomly permute a variable that does not gain you anything in prediction, then predictions won't change much and you will only see small changes in impurity and mse.  On the other hand the important variables will change the predictions by quite a bit if randomly permuted, so you will see bigger changes.  Turn this around and you see big changes indicate important variables.


#Rpart test:
fit <- rpart(SNPpheno~., data=vec, method="anova")
pdf(file="SNPTree.pdf")
fancyRpartPlot(fit)
dev.off()
summary(fit)
#Call:
#rpart(formula = SNPpheno ~ ., data = vec, method = "anova")
#  n= 41
#
#          CP nsplit rel error    xerror      xstd
#1 0.53472624      0 1.0000000 1.0558412 0.2147797
#2 0.11144099      1 0.4652738 0.9497579 0.2363781
#3 0.05690702      2 0.3538328 1.0061672 0.1906868
#4 0.01000000      3 0.2969257 0.9789176 0.1895580
#
#Variable importance
# SNPfam SNP1665  SNP179  SNP686  SNP725  SNP797 SNP2757 SNP1669 SNP2756 SNP2827
#     18      12      12      12      12      12       4       2       2       2
#SNP1772 SNP1623  SNP443 SNP1765   SNP12 SNP1444 SNP1651 SNP1757
#      2       2       2       1       1       1       1       1
#
#Node number 1: 41 observations,    complexity param=0.5347262
#  mean=84.69563, MSE=191.6899
#  left son=2 (12 obs) right son=3 (29 obs)
#  Primary splits:
#      SNPfam  splits as  LR,   improve=0.5347262, (0 missing)
#      SNP2791 splits as  RLL,  improve=0.4891549, (0 missing)
#      SNP463  splits as  RLL,  improve=0.4476991, (0 missing)
#      SNP2342 splits as  LRR,  improve=0.4418003, (0 missing)
#      SNP2741 splits as  LRLL, improve=0.4193188, (0 missing)
#  Surrogate splits:
#      SNP179  splits as  RRL, agree=0.902, adj=0.667, (0 split)
#      SNP686  splits as  LLR, agree=0.902, adj=0.667, (0 split)
#      SNP725  splits as  LR,  agree=0.902, adj=0.667, (0 split)
#      SNP797  splits as  LR,  agree=0.902, adj=0.667, (0 split)
#      SNP1665 splits as  LLR, agree=0.902, adj=0.667, (0 split)
#
#Node number 2: 12 observations
#  mean=68.95675, MSE=120.8659
#
#Node number 3: 29 observations,    complexity param=0.111441
#  mean=91.20827, MSE=76.08027
#  left son=6 (7 obs) right son=7 (22 obs)
#  Primary splits:
#      SNP2757 splits as  RLR, improve=0.3969702, (0 missing)
#      SNP1722 splits as  LLR, improve=0.3636026, (0 missing)
#      SNP365  splits as  RLR, improve=0.3545511, (0 missing)
#      SNP2857 splits as  RLR, improve=0.3496999, (0 missing)
#      SNP198  splits as  RLR, improve=0.3459826, (0 missing)
#  Surrogate splits:
#      SNP1669 splits as  LRR, agree=0.897, adj=0.571, (0 split)
#      SNP2756 splits as  RLR, agree=0.897, adj=0.571, (0 split)
#      SNP2827 splits as  RLR, agree=0.897, adj=0.571, (0 split)
#      SNP443  splits as  RLR, agree=0.862, adj=0.429, (0 split)
#      SNP1623 splits as  RRL, agree=0.862, adj=0.429, (0 split)
#
#Node number 6: 7 observations
#  mean=81.46561, MSE=67.43641
#
#Node number 7: 22 observations,    complexity param=0.05690702
#  mean=94.3082, MSE=39.01939
#  left son=14 (7 obs) right son=15 (15 obs)
#  Primary splits:
#      SNP1772 splits as  RLR, improve=0.5210097, (0 missing)
#      SNP561  splits as  RLR, improve=0.4538726, (0 missing)
#      SNP1621 splits as  LRR, improve=0.4176578, (0 missing)
#      SNP1651 splits as  LLR, improve=0.3918857, (0 missing)
#      SNP2518 splits as  LRL, improve=0.3788216, (0 missing)
#  Surrogate splits:
#      SNP1765 splits as  RLR, agree=0.909, adj=0.714, (0 split)
#      SNP12   splits as  RLR, agree=0.864, adj=0.571, (0 split)
#      SNP1444 splits as  RLL, agree=0.864, adj=0.571, (0 split)
#      SNP1651 splits as  RLR, agree=0.864, adj=0.571, (0 split)
#      SNP1757 splits as  RLR, agree=0.864, adj=0.571, (0 split)
#
#Node number 14: 7 observations
#  mean=87.70796, MSE=40.4887
#
#Node number 15: 15 observations
#  mean=97.38831, MSE=8.517137
#







#Ctree
 library(party)
 vecSEL = vec[vec$SNPfam=="SEL",]
fitSEL <- ctree(SNPpheno~., data=vecSEL)






###
#Genes in associated regions:

#1) Genes in high FST regions
load(file="DroneSelection.RData")

gff=read.table(file="/media/data1/forty3/brock/scripts/NCBIGFF.txt",col.names=c("chrom","GB","st","en","type"))
genes=gff[which(gff$type=="gene"),]
sts=aggregate(genes$st,by=list(genes$chrom),function(x) x)
ends=aggregate(genes$en,by=list(genes$chrom),function(x) x)
scaff.genes=aggregate(as.character(genes$GB),by=list(genes$chrom),function(x) as.character(x))
Find.High.Regions=cbind(gene=scaff.genes,gstart=sts,gend=ends,HighFstPositions98)
Fst.genes=apply(Find.High.Regions, 1, function(x) MultiWayOverlapper(x[7],x[7],x[4],x[6],x[2]) )
Fst.genes=lapply(Fst.genes,unique)
	#contains all genes for every chromosome

	
	
test=data.frame(GB=unlist(Fst.genes))
test=merge(genes,test,by="GB")
	#238 genes in high Fst regions.

	
map=read.table(file="ControlHighFST.map",header=F)
CONset=read.table("Candidates98_ALL.assoc.linear.set.mperm",header=T)
CONset=unlist(strsplit(as.vector((CONset$SNPS)),"[|]"));CONset=CONset[which(CONset!="NA")]
CONpos=map[map$V2 %in% CONset,]
SELset=read.table("Candidates98_ALLSEL.assoc.linear.set.mperm",header=T)
SELset=unlist(strsplit(as.vector((SELset$SNPS)),"[|]"));SELset=SELset[which(SELset!="NA")]
SELpos=map[map$V2 %in% SELset,]
	



















#2) Genes in/near associated SNPs










