#!/usr/bin/Rscript

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
		
		
	


args = commandArgs(trailingOnly = TRUE)
pedfile = args[1]	
mapfile = args[2]	
chr = 	args[3]	

#pedfile="DroneSamps_16.phasedMAF.ped"
#mapfile="DroneSamps_16.phasedMAF.map"
#
#

	
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
source(file="/media/data1/forty3/brock/scripts/VarFunct.r")

###

print("ready to run!")


###
#For Hygiene only

#Load in SNP set 
	#Selected SNPs:DroneSelection201Candidates.XXX"
	#Candidates (independant SNPs): CandidateSET 

map=read.table(file=mapfile)
dt = read.table(pedfile, colClasses=c(rep("character",nrow(map)+6)))
pheno=read.table(file="DronePhenoHB.txt",header=F)
pheno=pheno[,-c(2)]
vec=merge(pheno,dt,by="V1")
pop=read.table(file="SELCON.txt",header=F);pop=pop[,-c(2)]
vec=vec[,-c(3:7)]
names(vec)[2]="pheno"
nms=vec$V1
vec=vec[,-1]
phens=vec$pheno
vec=vec[,-1]

print("loaded") #I will do this again with haplotypes and NOT genotypes later (i.e. ranomdly select columns)

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



print("SFS")

#Optimize mtry
remo=ncol(vec)-1
bestmtry = tuneRF(vec[-remo], vec$SNPpheno,  ntreeTry=1000, stepFactor=0.5,improve=0.01, trace=TRUE, plot=FALSE, dobest=FALSE)
bt=data.frame(bestmtry)
mtry=as.numeric(bt[,1][which(bt[,2]==min(bt[,2]))])

####
#2) RF #need to optimize....
RFfit <- randomForest(SNPpheno ~ .,  data=vec, importance=TRUE, ntree=1000, mtry=mtry, proximity=T,localImp=T)


save.image(file=paste("testSFS", chr, ".RData",sep=""))



# Load the model and check importance ------------------

fils = list.files(pattern="testSFS*")
fils = fils[-16]
assoc = c()
for(i in fils){
	print(i)
	load(file=i)
	df = data.frame(RFfit$importance)
	snp.name.index = as.numeric(gsub("SNP", "", row.names(df[which(df$X.IncMSE>=5),])))
	snp.name.index = snp.name.index[!is.na(snp.name.index)] 
	out = map[snp.name.index,]
	out$imp = df$X.IncMSE[snp.name.index]
	assoc = rbind(out, assoc)
	}

#x11();varImpPlot(RFfit,scale=FALSE)
write.list(assoc, file="RPARTAssocSNPs")









