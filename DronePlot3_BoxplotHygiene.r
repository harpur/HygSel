####
# Plotting Selected HB vs Control
####




#Data
#From DroneRPARTv2.r
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



#Plotting:
#boxplot(SNPpheno~SNPfam, data=vec)
pdf("SelectionResult.pdf")
boxplot(SNPpheno~SNPfam, data=vec,
	#notch = T,
	col = c("gold","lightgreen"),
	ylab="",
	names=c("Control", "Selected"),
	boxwex = 0.75,
	staplewex = 0.5,
	whisklwd = 2,
	medlwd = 2,
	whisklty = "solid",
	staplewd = 2,
	outlwd = 1,
	boxlwd = 2,
	frame=FALSE,
	cex.lab=1.5, cex=1.5, cex.axis=1.5)
mtext("Percent Hygienic",side=2,line=2.4,cex=1.5)

means = aggregate((SNPpheno)~SNPfam,FUN="mean",data=vec)
points(c(1,2),c(means[1,2],means[2,2]),pch=23,cex=2,lwd=2,bg="white")
text(c(1,2),c(means[1,2]-1,means[2,2]-1),
    labels=format(c(means[1,2],means[2,2]),format="f",digits=3),
    pos=1,cex=1.2,col="black")	
dev.off()















