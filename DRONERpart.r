###
# Analyses with Rpart (Teal window).
###
























#Ran this to get out candidate high-FST SNPs in Control population.
#plink --file DroneSelection --noweb --allow-no-sex --recode --out HLCandidateCAT  --extract GenicHighestFSTSNPs.snp --keep controlBees.txt


#Recall, data #plink --file DroneSelection --mpheno 1 --pheno DronePhenoHB.txt --noweb --allow-no-sex --linear --adjust  --out HLCandidate --keep controlBees.txt --extract GenicHighestFSTSNPs.snp
write.table(test[c("CHROM","POS")],file="CandidateSNPs.snp",quote=F,row.names=F,col.names=F)

plink --file DroneSelection --extract CandidateSNPs2.snp --out DroneSelection201Candidates --noweb --make-bed --keep controlBees.txt 
	#

plink --bfile DroneSelection201Candidates --recode --out DroneSelection201CandidatesP --noweb





plink --file DroneSelection --mpheno 1 --pheno DronePhenoHB.txt --noweb --allow-no-sex --linear --adjust  --out HLCandidate --keep controlBees.txt --extract CandidateSNPs.snp


test=(HighSNPs[HighSNPs$SNP %in% beta$SNP,])


#It takes too long for all SNPs, going to try set of 99th percentile FST SNPs.
			#dt = read.table("HLCandidateCAT.ped")
		#	x = grep("^3", dt)
		#	x = x[seq(1,80,2)];x=c(x, length(dt)+1) #this value is length of dt+1
		#	vec = matrix(nr = 40,nc =  26141) #nc=number of lines in test+1
		#	for(i in 1:40){
		#		#i =40
		#		test = dt[x[i]:(x[i+1]-1)]
		#		#remove 2:6,the comment BS
		#		test = test[-c(2:6)]
		#		vec[i,] = test
		#		}
		#	vec[,1]

map=read.table(file="HLCandidateCAT.map")
dt = read.table("HLCandidateCAT.ped", colClasses=c(rep("character",432)))
pheno=read.table(file="DronePhenoHB.txt",header=F)
pheno=pheno[,-c(2)]
vec=merge(pheno,dt,by="V1")
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


library("rpart")
fit <- rpart(pheno~., data=vec, method="class", control=rpart.control(minsplit=2, minbucket=1, cp=0.01))
plot(fit);text(fit)





numvec=vec
id <- c(1:ncol(numvec))
numvec[,id] <- as.numeric(as.character(unlist(numvec[,id])))
fit <- rpart(pheno~., data=numvec, method="anova", control=rpart.control(minsplit=2, minbucket=1, cp=0.001))
x11();plot(fit);text(fit)
pfit<- prune(fit, cp=fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])
plot(pfit)
text(pfit)
summary(pfit)






