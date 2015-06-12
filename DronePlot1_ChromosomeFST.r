######
# Drone Association Plotting Fst, and Association Set Test 
######











##PLOTTING

##OK
#R
#library
library(ggplot2)


#Datasets
load(file="/media/data1/forty3/drone/FST/SelvsCon/SelectedvsControlFST.RData")
map=read.table(file="ControlHighFST.map",header=F)
map=read.table(file="ControlHighFST.map",header=F)

#Load in significantly associated SNPs
	#CONset=read.table("Candidates98_ALL.assoc.linear.set.mperm",header=T)
	#CONset=unlist(strsplit(as.vector((CONset$SNPS)),"[|]"));CONset=CONset[which(CONset!="NA")]
	#CONpos=map[map$V2 %in% CONset,]
	#chromcon=unique(CONpos$V1)
	#SELset=read.table("Candidates98_ALLSEL.assoc.linear.set.mperm",header=T)
	#SELset=unlist(strsplit(as.vector((SELset$SNPS)),"[|]"));SELset=SELset[which(SELset!="NA")]
	#SELpos=map[map$V2 %in% SELset,]
	#chromsel=unique(SELpos$V1)
	#
RFset=read.table("AssociatedRFSNPs",header=F)



#Covnerted Oxley's hygien QTLs to loaction in AMEL:
#	chr2	14249515	15407932
#2	chr5	9783962	10814860
#3	chr16	2058269
#4	chr16	3127046





setwd("/media/data1/forty3/drone/R")
####################
##Plotting:
pdf(file="FSTSelectTest.Chromosome1_8_RFSNPs.pdf")
par(mfrow = c(8, 1),
	cex = 0.4, 
	mar = c(2, 2, 2, 2), 
	oma = c(1, 1, 1, 1)
	)

for(i in 1:6){
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
	
#	#Add in RFSNPs
#	xpos=RFset$V4[RFset$V1==i]
#	ypos=rep(0.2,length(xpos));ypos=ypos+(1:length(xpos)/100)
#	points(xpos,ypos,cex=2,bg="coral3",col="black",pch=21) #control 
#	
#	##Add Points 
#	#xpos=CONpos$V4[CONpos$V1==i]
#	#ypos=rep(0.2,length(xpos));ypos=ypos+(1:length(xpos)/100)
#	#points(xpos,ypos,cex=2,bg="coral3",col="black",pch=21) #control 
#	#
#	##Add in Selection-associated
#	#xpos=SELpos$V4[SELpos$V1==i]
#	#ypos=rep(0.25,length(xpos));ypos=ypos+(1:length(xpos)/100)
#	#points(xpos,ypos,cex=2,bg="coral4",col="black",pch=24) #selected
#	#
#	
#	
#	#Add in QTLs:
#	if(i==5){
#		lines(y=c(0.03,0.03),x=c(9783962,10814860),lwd=3)
#	}else{
#	if(i==2){
#		lines(y=c(0.03,0.03),x=c(14249515,15407932),lwd=3)
#	}else{
#	if(i==16){
#		lines(y=c(0.03,0.03),x=c(2058269,3127046),lwd=3)
#	}else{		
#		
			}
		}
	}
}

dev.off()










#CHR 5 only:

pdf(file="CHR5.pdf")
g5=genes[genes$chrom=="5",]
for(i in 5){
	#Main Plot
	plot(as.numeric(unlist(Fst$Pos[Fst$Group.1==i])),as.numeric(unlist(Fst$x[Fst$Group.1==i])), 
	type="l",
	xaxt="n",
	xlab=NA,
	ylab=NA,
	xlim=c(9000000,  13000000),
	ylim=c(0,0.35)
	)
	xpos <- seq(0, max(as.numeric(unlist(Fst$Pos[Fst$Group.1==i]))), by=1000000)
	axis(1, at=xpos,labels=xpos/1000)
	axis(2, pos=0, at=c(0,0.17, 0.35), labels=c("0", "0.175", "0.35"))
	mtext("Chromosome 5 (kb)",side = 1, line = 3,  cex = 1)
	mtext("Fixation Index (Fst)",side = 2, line = 3,  cex = 1)
	lines(y=rep(highFstWindow,length(as.numeric(unlist(Fst$Pos[Fst$Group.1==i])))),
		x=as.numeric(unlist(Fst$Pos[Fst$Group.1==i])), 
		col="red",
		lty=2)
	lines(y=c(0.03,0.03),x=c(9783962,10814860),lwd=3)
	points(x=g5$st, y=rep(0.2,length(g5$st)),cex=0.1)
	
	#add in perm SNPs:
	xpos=c(10867357)
	ypos=rep(0.2,length(xpos));ypos=ypos+(1:length(xpos)/100)
	points(xpos,ypos,cex=2,bg="coral3",col="black",pch=21)
	
	#add in perm SNPs (Selected):
	xpos=c(10904460)
	ypos=rep(0.2,length(xpos));ypos=ypos+(1:length(xpos)/100)
	points(xpos,ypos,cex=2,bg="coral4",col="black",pch=24) 
	
	
}
dev.off()

#All
#582524 5.14:1610939 0.1822910     5 10867357
#582580 5.14:1612227 0.3507480     5 10868645
#582597 5.14:1614317 0.0589269     5 10870735
#583827 5.14:1734865 0.5027480     5 10991283
#584099 5.14:1759660 0.5333270     5 11016078



#CON
#5.14:1610939

#SEL
#5.14:1648042



