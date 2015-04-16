###
# Analyses with Rpart (Teal window).
###









#Along with DronePLINKset.sh, this identifies associated SNPs with hygiene and selection candidates
	#Uses CandidateSNPs.snp 
	#And XXXXXXXXXXX the set of independant SNPs in that list







#Ran this to get out candidate high-FST SNPs in Control population.
#plink --file DroneSelection --noweb --allow-no-sex --recode --out HLCandidateCAT  --extract GenicHighestFSTSNPs.snp --keep controlBees.txt


#Recall, data #plink --file DroneSelection --mpheno 1 --pheno DronePhenoHB.txt --noweb --allow-no-sex --linear --adjust  --out HLCandidate --keep controlBees.txt --extract GenicHighestFSTSNPs.snp
write.table(test[c("CHROM","POS")],file="CandidateSNPsp10.snp",quote=F,row.names=F,col.names=F)
write.table(test$SNP,file="CandidateSNPs.snp",quote=F,row.names=F,col.names=F)





plink --file DroneSelection --extract CandidateSNPs.snp --out DroneSelection201Candidates --noweb --make-bed 
	#--keep controlBees.txt 

plink --bfile DroneSelection201CandidatesP50  --recode --out DroneSelection201Candidates --noweb


map=read.table(file="DroneSelection201Candidates.map")
dt = read.table("DroneSelection201Candidates.ped", colClasses=c(rep("character",nrow(map)+6)))
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

library("rpart")
fit <- rpart(pheno~., data=vec, method="class", control=rpart.control(minsplit=2, minbucket=1, cp=0.01))
plot(fit);text(fit)


printcp(fit)

#Regression tree:
#rpart(formula = pheno ~ ., data = numvec, method = "anova", control = rpart.control(minsplit = 2,
#    minbucket = 1, cp = 0.001))
#
#Variables actually used in tree construction:
#[1] 1   103 107 13  16  25  26  33
#
#Root node error: 1450.4/12 = 120.87
#
#n= 12
#
#          CP nsplit  rel error xerror    xstd
#1  0.7490576      0 1.00000000 1.1350 0.28759
#2  0.1003379      1 0.25094235 1.7973 0.54135
#3  0.0716332      2 0.15060446 1.7281 0.52035
#4  0.0247426      3 0.07897123 1.5655 0.53045
#5  0.0244270      4 0.05422865 1.5614 0.60837
#6  0.0224909      5 0.02980163 1.5614 0.60837
#7  0.0034890      6 0.00731076 1.5387 0.61259
#8  0.0021909      7 0.00382171 1.5621 0.61006
#9  0.0013473      8 0.00163083 1.5621 0.61006
#10 0.0010000      9 0.00028349 1.5716 0.61405
#









#Functional Now:



#numvec=vec[-203]
#id <- c(1:ncol(numvec))
#numvec[,id] <- as.numeric(as.character(unlist(numvec[,id])))
fit <- rpart(pheno~., data=vec, method="anova", control=rpart.control(minsplit=2, minbucket=1, cp=0.001))
x11();plot(fit);text(fit)
pfit<- prune(fit, cp=fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])

plot(pfit)
text(pfit)
summary(pfit)



dfit <- rpart(pheno~., data=numvec, method="anova", control=rpart.control(minsplit=2, minbucket=1, cp=0.001),maxdepth=2)












#FOLLOWING THIS BLOG
	#http://trevorstephens.com/post/72916401642/titanic-getting-started-with-r

#With the 201 candidates:
write.table(test$SNP,file="CandidateSNPs.snp",quote=F,row.names=F,col.names=F)
plink --file DroneSelection --extract CandidateSNPs.snp --out DroneSelection201Candidates --noweb --make-bed 
	#--keep controlBees.txt 
plink --bfile DroneSelection201Candidates  --recode --out DroneSelection201Candidates --noweb


####################################################
# While I've got these files, I want to try some GWAS statsy stuff
#######
#file="FSTDroneSelectionCandidatesALL.map"

plink --file FSTDroneSelectionCandidatesALL --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --remove controlBees.txt --allow-no-sex --assoc --qt-means --adjust 


R
test=read.table("plink.qassoc.adjusted",header=T)
source("/media/data1/forty3/brock/scripts/qq.r")
qq(test$UNADJ)


###############################################




#With just HighSNPs:
write.table(HighSNPs$SNP,file="FSTCandidateSNPs.snp",quote=F,row.names=F,col.names=F)
plink --file DroneSelection --extract FSTCandidateSNPs.snp --out FSTDroneSelectionCandidatesALL --noweb --make-bed #--keep controlBees.txt 
plink --bfile FSTDroneSelectionCandidatesALL  --recode --out FSTDroneSelectionCandidatesALL --noweb





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



#All SNPs
map=read.table(file="FSTDroneSelectionCandidatesALL.map")
dt = read.table("FSTDroneSelectionCandidatesALL.ped", colClasses=c(rep("character",nrow(map)+6)))
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



#Selected SNPs:
map=read.table(file="DroneSelection201Candidates.map")
dt = read.table("DroneSelection201Candidates.ped", colClasses=c(rep("character",nrow(map)+6)))
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

#Rpart test:
fit <- rpart(SNPpheno~., data=vec, method="anova")
x11();fancyRpartPlot(fit)



fit <- rpart(SNPpheno~.*famam, data=vec, method="anova", control=rpart.control(minsplit=2, minbucket=1, cp=0.001))
x11();fancyRpartPlot(fit)
pfit<- prune(fit, cp=fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])
x11();fancyRpartPlot(pfit)











#Random Forest:
names(vec)=paste("SNP",names(vec),sep="")
set.seed(213)




RFfit <- randomForest(SNPpheno ~ .,  data=vec, importance=TRUE, ntree=200, mtry=100)
 x11();varImpPlot(RFfit)

 
 
#In Parallel:
vec.SEL=vec[vec$SNPfam=="SEL",]


names(vec)=paste("SNP",names(vec),sep="")
set.seed(213)

impot=c()
for(i in 1:100){
	vec.test=vec[c(sample(c(1:36510), 500, replace=T), 36511, 36512)]
	rf <- foreach(ntree=rep(300, 12), .combine='combine', .packages='randomForest', .multicombine=TRUE) %dopar% randomForest(pheno ~ .,  data=vec.test, importance=TRUE, ntree=ntree, mtry=100)
impot=rbind(impot,importance(rf))
}


impot=impot[impot[,2]>100,]
#unique(row.names(impot))


pheno = vec.SEL$SNPpheno
vec.test = vec.SEL[ , names(vec.SEL) %in% c(unique(row.names(impot)))]
vec.test$pheno =  pheno
 

fitSEL <- rpart(SNPpheno~., data=vecSEL, method="anova", control=rpart.control(minsplit=2, minbucket=1, cp=0.001))
fitCON <- rpart(SNPpheno~., data=vecCON, method="anova", control=rpart.control(minsplit=2, minbucket=1, cp=0.001))

 

 
 
 
 
 
 
 
 
 
 
 #CTREE:
 library(party)
 
fitSEL <- ctree(SNPpheno~., data=vecSEL)
fitCON <- ctree(SNPpheno~., data=vecCON)
 
 
 
 
 
 
 
 
 
 



