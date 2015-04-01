####################################################
# RPart with RF and Linked SNPs
#######
#file="ControlHighFST.map"
#http://pngu.mgh.harvard.edu/~purcell/plink/anal.shtml#set



#for each High Fst Region, I'm going to 
	#1) output those that are independant
	#2) run RF to see which are important
	#3) Using Important SNPs, RPART


###	
#1) Get out LINKED SNPs within Sets:	

plink --file ControlHighFST  --ld-snp-list highFST.list  --noweb --allow-no-sex --out CONTROL

plink --file SelHighFST  --ld-snp-list highFST.list  --noweb --allow-no-sex --out SELECT



plink --file DroneSelection --indep-pairwise 50 5 0.5 --extract highFST.list --noweb --allow-no-sex  --keep controlBees.txt --out CONTROL
plink --file DroneSelection --indep-pairwise 50 5 0.5 --extract highFST.list --noweb --allow-no-sex  --remove controlBees.txt --out SELECT

plink --file DroneSelection --indep-pairwise 50 5 0.5 --extract highFST.list --noweb --allow-no-sex --out ALL

plink --file DroneSelection --extract ALL.prune.in --out CandidateSET --noweb --make-bed 
plink --bfile CandidateSET  --recode --out  CandidateSET --noweb



#2) Using these candidates (CandidateSET), run RF 

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
map=read.table(file="CandidateSET.map")
dt = read.table("CandidateSET.ped", colClasses=c(rep("character",nrow(map)+6)))
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
#RFfit <- randomForest(SNPpheno ~ .,  data=vec, importance=TRUE, ntree=200, mtry=100)
#x11();varImpPlot(RFfit)

rf <- foreach(ntree=rep(600, 15), .combine='combine', .packages='randomForest', .multicombine=TRUE) %dopar% randomForest(SNPpheno ~ .,  data=vec, importance=TRUE, ntree=ntree, mtry=100)
 x11();varImpPlot(rf)
 
 #I kept all with Inc%MSE >1 
row.names(test[test[,1]>1,])
pheno = vec$SNPpheno
vec.test = vec[ , names(vec) %in% c(row.names(test[test[,1]>1,]))]
vec.test$pheno =  pheno


fit <- rpart(pheno~., data=vec.test, method="anova")
x11();fancyRpartPlot(fit)


#Rpart test:
fit <- rpart(SNPpheno~., data=vec, method="anova")
x11();fancyRpartPlot(fit)



fit <- rpart(SNPpheno~.*famam, data=vec, method="anova", control=rpart.control(minsplit=2, minbucket=1, cp=0.001))
x11();fancyRpartPlot(fit)
pfit<- prune(fit, cp=fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])
x11();fancyRpartPlot(pfit)





#Test with ONLY control bees:
vecCON = vec[vec$SNPfam=="CON", ]; vecCON$SNPfam=NULL
fitCON <- rpart(SNPpheno~., data=vecCON, method="anova")
x11();fancyRpartPlot(fitCON)
 rf <- foreach(ntree=rep(600, 15), .combine='combine', .packages='randomForest', .multicombine=TRUE) %dopar% randomForest(SNPpheno ~ .,  data=vecCON, importance=TRUE, ntree=ntree, mtry=100)
 x11();varImpPlot(rf)



#Test with ONLY SEL bees:
vecSEL = vec[vec$SNPfam=="SEL", ]; vecSEL$SNPfam=NULL
fitSEL <- rpart(SNPpheno~., data=vecSEL, method="anova")
x11();fancyRpartPlot(fitSEL)




#Random Forest:
set.seed(213)


#RFfit <- randomForest(SNPpheno ~ .,  data=vec, importance=TRUE, ntree=200, mtry=100)
#x11();varImpPlot(RFfit)


 rf <- foreach(ntree=rep(600, 15), .combine='combine', .packages='randomForest', .multicombine=TRUE) %dopar% randomForest(SNPpheno ~ .,  data=vec, importance=TRUE, ntree=ntree, mtry=100)
 x11();varImpPlot(rf)
 
 
 
 
#In Parallel:
vec.SEL=vec[vec$SNPfam=="SEL",]


names(vec)=paste("SNP",names(vec),sep="")
set.seed(213)

impot=c()
for(i in 1:100){
	vec.test=vec[c(sample(c(1:2784), 200, replace=T), 2785, 2786)]
	rf <- foreach(ntree=rep(300, 12), .combine='combine', .packages='randomForest', .multicombine=TRUE) %dopar% randomForest(SNPpheno ~ .,  data=vec.test, importance=TRUE, ntree=ntree, mtry=100)
impot=rbind(impot,importance(rf))
}


impot=impot[impot[,2]>300,]
x=row.names(impot);x=gsub(".1","",x)
unique(x)


pheno = vec$SNPpheno
vec.test = vec[ , names(vec) %in% c(unique(x))]
vec.test$pheno =  pheno
 

 
 #fit <- randomForest(pheno ~ .,  data=vec.test, importance=TRUE, ntree=200, mtry=100)
 #varImpPlot(fit)
 
fit <- rpart(pheno~., data=vec.test, method="anova")
#pdf("fitout.pdf")
fancyRpartPlot(fit)
#dev.off()
x11();fancyRpartPlot(fit)
 
 
 

#
x11();plot(fit);text(fit)
pfit<- prune(fit, cp=fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])

plot(pfit)
text(pfit)
summary(pfit)

























##OK

R
load(file="/media/data1/forty3/drone/FST/SelvsCon/SelectedvsControlFST.RData")


map=read.table(file="ControlHighFST.map",header=F)
mns=read.table("CONTROLSET.qassoc.means",header=T)


CONset=read.table("CONTROLSET.qassoc.set.mperm",header=T)
CONset=unlist(strsplit(as.vector((CONset$SNPS)),"[|]"));CONset=CONset[which(CONset!="NA")]
CONpos=map[map$V2 %in% CONset,]
chroms=unique(CONpos$V1)
# 3  5  6  7  9 11 12 14 15 16

SELset=read.table("SelSET.qassoc.set.mperm",header=T)
SELset=unlist(strsplit(as.vector((SELset$SNPS)),"[|]"));SELset=SELset[which(SELset!="NA")]
SELpos=map[map$V2 %in% SELset,]


