###
# Analyses with Rpart 
###

#Along with DronePLINKset.sh, this identifies associated SNPs with hygiene and selection candidates
	#Uses CandidateSNPs.snp 
	#And CandidateSET.xxx the set of independant SNPs in that list (--indep-pairwise)


#for each High Fst Region, I'm going to 
		#1) output those that are independant - DronePLINKset.sh
		#2) run RF to see which are important - here

		
#RF notes:
#http://www.statistik.uni-dortmund.de/useR-2008/slides/Strobl+Zeileis.pdf
#http://www.biomedcentral.com/1471-2105/10/78
#http://bib.oxfordjournals.org/content/early/2012/07/10/bib.bbs034.full.pdf+html
#http://www.bios.unc.edu/~dzeng/BIOS740/randomforest.pdf (Liaw and Wiener, 2002)
		#Ensemble learnign method - generates many classifiers and aggregates their results
		#each tree is independant 
		#Node is split based on best of a subset of data
		#At each iteration of the bootstrap (i.e. after tree building), predict the remaining (OOB data)
		
		
	

	
	
	
	
#RF analyses	
	
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




###
#For Hygiene only

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

#remove missing data
test=apply(vec,2,function(x) length(x[x=="00"]))
	#46 SNPs with missing data
test=which(test>0)
vec=vec[-(test)]




#1) optimize MTRY	
set.seed(1234)
bestmtry = tuneRF(vec[-2739], vec$SNPpheno,  ntreeTry=1000, stepFactor=0.5,improve=0.01, trace=TRUE, plot=TRUE, dobest=FALSE)
	#best 456

####
#2) RF
set.seed(22453)
RFfit <- randomForest(SNPpheno ~ .,  data=vec, importance=TRUE, ntree=1000, mtry=456, proximity=T,localImp=T)

##
#3a) accuracy of predictions
Predicted=predict(RFfit, vec[-2739])
cor.test(vec$SNPpheno, Predicted)
	#data:  vec$SNPpheno and Predicted
	#t = 32.4195, df = 39, p-value < 2.2e-16
	#alternative hypothesis: true correlation is not equal to 0
	#95 percent confidence interval:
	# 0.9661754 0.9904012
	#sample estimates:
	#      cor
	#0.9819475
	#
	
#x11();varImpPlot(RFfit,scale=FALSE)
pdf(file="SNPRFImportance.pdf")
plot(vec.CON$SNPpheno, Predicted.CON,
	pch = 19,
	col="black",
	xlab = "Observed Hygienic %",
	ylab = "Predicted Hygienic %"
	)
dev.off()





##
#3b) accuracy of predictions on SELECTED only
vec.CON=vec[vec$SNPfam=="CON",]
vec.SEL=vec[vec$SNPfam=="SEL",]
RFfit.CON <- randomForest(SNPpheno ~ .,  data=vec.CON, importance=TRUE, ntree=1000, mtry=456, proximity=T,localImp=T)
Predicted.CON=predict(RFfit.CON, vec.CON[-2739])
cor.test(vec.CON$SNPpheno, Predicted.CON)

#        Pearson's product-moment correlation
#
#data:  vec.SEL$SNPpheno and Predicted.SEL
#t = 35.7862, df = 27, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9777479 0.9951756
#sample estimates:
#      cor
#0.9896223
#


#Con predicting SEL
RFfit.CON <- randomForest(SNPpheno ~ .,  data=vec.CON, importance=TRUE, ntree=1000, mtry=456, proximity=T,localImp=T)
Predicted.CON=predict(RFfit.CON, vec.SEL[-2739])
cor.test(vec.SEL$SNPpheno, Predicted.CON)


#Sel predicting CON

RFfit.SEL <- randomForest(SNPpheno ~ .,  data=vec.CON, importance=TRUE, ntree=1000, mtry=456, proximity=T,localImp=T)
Predicted.SEL=predict(RFfit.SEL, vec.CON[-2739])
cor.test(vec.CON$SNPpheno, Predicted.SEL)
x11();varImpPlot(RFfit.SEL,scale=FALSE)


#plotting:
pdf(file="SelpredictingCon.pdf")
plot(vec.CON$SNPpheno, Predicted.SEL,pch=19, 
	ylab="Observed Control Population Hygienic Behaviour (%)",
	xlab="Predicted Control Population Hygienic Behaviour (%)"
)
cor.test(vec.CON$SNPpheno, Predicted.SEL)
dev.off()


#importance plot
imps=(RFfit$localImportance)
imps=imps[,c(1:29)]
test=apply(imps,1,mean)
test=test[-2739]

imps = data.frame(test)
names(imps)="IncMSE"
imps$SNP=row.names(imps)
imps=imps[order(-imps$IncMSE),]

#plotting
pdf(file="SNPRFImportance.pdf")
ord = rev(order(imps$IncMSE, decreasing = TRUE)[1:30])
main = row.names(imps)[ord]	
main[main=="SNPfam"]="Control vs Selected"
dotchart(imps$IncMSE[ord],
	labels=main,
	xlab = "Random Forest Importance",
	pch=19
	)
abline(v=1,lty=2)
dev.off()

	

##
#4) Importance plot

imps = data.frame(importance(RFfit.SEL, type=1, scale=FALSE))
names(imps)="IncMSE"
imps$SNP=row.names(imps)
imps=imps[order(-imps$IncMSE),]

#plotting
pdf(file="SNPRFImportance.pdf")
ord = rev(order(imps$IncMSE, decreasing = TRUE)[1:30])
main = row.names(imps)[ord]	
main[main=="SNPfam"]="Control vs Selected"
dotchart(imps$IncMSE[ord],
	labels=main,
	xlab = "Random Forest Importance",
	pch=19
	)
abline(v=0.5,lty=2)
dev.off()
#I also plot these on my FST figure...
candidates = row.names(imps)[ord][imps$IncMSE[ord]>1]	
candidates = as.numeric(gsub("SNP","",candidates))
candidates = candidates[!is.na(candidates)]
map[candidates,]
#     V1            V2 V3       V4
#1491  9  9.10:2424861  0  6844829
#1124  7    7.13:37175  0  4472042
#1132  7    7.13:43100  0  4477967
#1972 11 11.18:1366430  0 10271050
#2066 12 12.13:2337671  0  5368151
#232   3   3.8:1814662  0  5781488
#1982 11 11.18:1370550  0 10275170
#1950 11 11.18:1349361  0 10253981
#1485  9  9.10:2419016  0  6838984
#1933 11 11.18:1340362  0 10244982
#80    1    1.32:56464  0 19745368
#1595  9   9.12:738270  0  9987666

#plotted these with DronePlot1_ChromosomeFST.r









##############################
#local importance later...
imps = data.frame(RFfit$localImportance)

names(imps)="IncMSE"
imps$SNP=row.names(imps)
imps=imps[order(-imps$IncMSE),]

#plotting
ord = rev(order(imps$IncMSE, decreasing = TRUE)[1:30])
main = row.names(imps)[ord]	
main[main=="SNPfam"]="Control vs Selected"
dotchart(imps$IncMSE[ord],
	labels=main,
	xlab = "Random Forest Importance",
	pch=19
	)


importance(RFfit, type=1, scale=FALSE)
#importance(RFfit, type=2, scale=FALSE)
x11();varImpPlot(RFfit,scale=FALSE)
x11();boxplot(vec$SNPpheno~vec$SNP1905*vec$SNPfam)








##
#4)Getting significance of imp values
impot=c()
for(i in 1:10000){
	vec.test=vec
	vec.test$SNPpheno=sample(vec$SNPpheno, nrow(vec.test),replace=F) 
	rf <- foreach(ntree=rep(300, 12), .combine='combine', .packages='randomForest', .multicombine=TRUE) %dopar% randomForest(SNPpheno ~ .,  data=vec.test, importance=TRUE, ntree=ntree, mtry=100)
impot=rbind(impot,importance(rf))

}














#Call:
# randomForest(formula = SNPpheno ~ ., data = vec, importance = TRUE,      ntree = 1000, mtry = 456, proximity = T, localImp = T)
#               Type of random forest: regression
#                     Number of trees: 1000
#No. of variables tried at each split: 456
#
#          Mean of squared residuals: 147.8011
#                    % Var explained: 22.9
#


test=as.matrix(RFfit$localImportance)
test=t(test)
test=as.matrix(test)
test=data.frame(test)
test$fam=vec$SNPfam
test$SNPfam=NULL
#made "imps"

plot(imps$con,imps$sel,col="white")
text(imps$con,imps$sel,labels=imps$SNP)
imps[imps$sel>1.5,]
imps[imps$con>1.5,]

x11();boxplot(vec$SNPpheno~vec$SNP1905*vec$SNPfam)


#what to include? http://iospress.metapress.com/content/r36932421476k620/fulltext.pdf


############
RFfit <- randomForest(SNPpheno ~ .,  data=vec, importance=TRUE, ntree=10000, mtry=200, keep.forest=FALSE)
x11();varImpPlot(RFfit)
x11();boxplot(vec$SNPpheno~vec$SNP1905*vec$SNPfam)
	#Call:
	# randomForest(formula = SNPpheno ~ ., data = vec, importance = TRUE,      ntree = 1000, mtry = 232)
	#               Type of random forest: regression
	#                     Number of trees: 1000
	#No. of variables tried at each split: 232
	#
	#          Mean of squared residuals: 144.3785
	#                    % Var explained: 24.68
	#
#OOB error rate, model accuracy


	
	

#In Parallel:
vec.SEL=vec[vec$SNPfam=="SEL",]
set.seed(213)

impot=c()
for(i in 1:100){
	vec.test=vec[c(sample(c(1:2784), 200, replace=T), 2785, 2786)]
	rf <- foreach(ntree=rep(300, 12), .combine='combine', .packages='randomForest', .multicombine=TRUE) %dopar% randomForest(SNPpheno ~ .,  data=vec.test, importance=TRUE, ntree=ntree, mtry=100)
impot=rbind(impot,importance(rf))
}



rf <- foreach(ntree=rep(25000, 6), .combine=combine, .multicombine=TRUE,
              .packages='randomForest') %dopar% {
    randomForest(x, y, ntree=ntree)
}






#For both sets together:

 
 
 
 
 
impot = importance(RFfit)
impot=impot[impot[,2]>50,]
vec.test = vec[ , names(vec) %in% c(unique(row.names(impot)))]
vec.test$SNPpheno =  pheno$V3
#Rpart test with most important SNPs
fit <- rpart(SNPpheno~., data=vec.test, method="anova")
x11();fancyRpartPlot(fit)



pheno = vec.SEL$SNPpheno
vec.test = vec.SEL[ , names(vec.SEL) %in% c(unique(row.names(impot)))]
vec.test$pheno =  pheno






###
#For Protein Markers:
#Load in SNP set 
	#Selected SNPs:DroneSelection201Candidates.XXX"
	#Candidates (independant SNPs): CandidateSET 
map=read.table(file="CandidateSET.map")
dt = read.table("CandidateSET.ped", colClasses=c(rep("character",nrow(map)+6)))
pheno=read.table(file="/media/data1/forty3/drone/vcf_drone/HBexpression.txt",header=F)
pheno=pheno[,c(1,5)] #chaaange
vec=merge(pheno,dt,by="V1")
pop=read.table(file="DroneSelection201Candidates.fam",header=F)
pop=pop[,c(1,3)];names(pop)=c("nms","fam")
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
vec$nms=nms
vec=merge(vec,pop,by="nms")
names(vec)=paste("SNP",names(vec),sep="")

#Rpart test:
fit <- rpart(SNPpheno~., data=vec, method="anova")
x11();fancyRpartPlot(fit)

















































 
#CTREE:
library(party)
fitC<- ctree(SNPpheno~., data=vec)
x11(); plot(fitC) 
 
 
 
 
################################### 
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
 
 
 



