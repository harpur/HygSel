######
# Drone Association Plotting Fst, and Association Set Test 
######











##PLOTTING

##OK
#R
#library
#library(ggplot2)


#Datasets
load(file="/media/data1/forty3/drone/FST/SelvsCon/SelectedvsControlFST.RData")
map=read.table(file="ControlHighFST.map",header=F)
map=read.table(file="ControlHighFST.map",header=F)
CONset=read.table("Candidates98_ALL.assoc.linear.set.mperm",header=T)
CONset=unlist(strsplit(as.vector((CONset$SNPS)),"[|]"));CONset=CONset[which(CONset!="NA")]
CONpos=map[map$V2 %in% CONset,]
chromcon=unique(CONpos$V1)
SELset=read.table("Candidates98_ALLSEL.assoc.linear.set.mperm",header=T)
SELset=unlist(strsplit(as.vector((SELset$SNPS)),"[|]"));SELset=SELset[which(SELset!="NA")]
SELpos=map[map$V2 %in% SELset,]
chromsel=unique(SELpos$V1)

#Covnerted Oxley's hygien QTLs to loaction in AMEL:
#	chr2	14249515	15407932
#2	chr5	9783962	10814860
#3	chr16	2058269
#4	chr16	3127046





setwd("/media/data1/forty3/drone/R")
####################
##Plotting:
pdf(file="FSTSelectTest.Chromosome1_8_updated.pdf")
par(mfrow = c(8, 1),
	cex = 0.4, 
	mar = c(2, 2, 2, 2), 
	oma = c(1, 1, 1, 1)
	)

for(i in 1:8){
	#Main Plot
	plot(as.numeric(unlist(Fst$Pos[Fst$Group.1==i])),as.numeric(unlist(Fst$x[Fst$Group.1==i])), 
	type="l",
	axes=FALSE,
	xlab=NA,
	ylab=NA,
	xlim=c(0,  30000000),
	ylim=c(0,0.35)
	)
	xpos <- seq(0, max(as.numeric(unlist(Fst$Pos[Fst$Group.1==i]))), by=1000000)
	axis(1, at=xpos,labels=xpos/1000)
	axis(2, pos=0, at=c(0,0.17, 0.35), labels=c("0", "0.175", "0.35"))
	mtext(i,side = 3, line = -1, adj = 0.05, cex = 1)
	lines(y=rep(highFstWindow,length(as.numeric(unlist(Fst$Pos[Fst$Group.1==i])))),
		x=as.numeric(unlist(Fst$Pos[Fst$Group.1==i])), 
		col="red",
		lty=2)
	
	#Add Points 
	xpos=CONpos$V4[CONpos$V1==i]
	ypos=rep(0.2,length(xpos));ypos=ypos+(1:length(xpos)/100)
	points(xpos,ypos,cex=2,bg="coral3",col="black",pch=21) #control 
	
	#Add in Selection-associated
	xpos=SELpos$V4[SELpos$V1==i]
	ypos=rep(0.25,length(xpos));ypos=ypos+(1:length(xpos)/100)
	points(xpos,ypos,cex=2,bg="coral4",col="black",pch=24) #selected
	
	#Add in QTLs:
	if(i==5){
		lines(y=c(0.03,0.03),x=c(9783962,10814860),lwd=3)
	}else{
	if(i==2){
		lines(y=c(0.03,0.03),x=c(14249515,15407932),lwd=3)
	}else{
	if(i==16){
		lines(y=c(0.03,0.03),x=c(2058269,3127046),lwd=3)
	}else{		
		
			}
		}
	}
}

dev.off()

