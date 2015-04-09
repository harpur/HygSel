###
# Analyses with Rpart 
###

#Along with DronePLINKset.sh, this identifies associated SNPs with hygiene and selection candidates
	#Uses CandidateSNPs.snp 
	#And CandidateSET.xxx the set of independant SNPs in that list (--indep-pairwise)


#for each High Fst Region, I'm going to 
		#1) output those that are independant - DronePLINKset.sh
		#2) run RF to see which are important - here
		#3) Using Important SNPs, RPART - here
		
	

	
	
	
	
	
	
###	
#Libraries and Functions
library(randomForest)
library(foreach)
library(rattle)
library(rpart.plot)
library(RColorBrewer)
library(rpart)
library(doMC)
registerDoMC(15)
###






#Load in SNP set 
	#Selected SNPs:DroneSelection201Candidates.XXX"
	#Candidates (independant SNPs): CandidateSET 
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
 
 
 



