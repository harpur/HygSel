###
# Plot significant Boxplot of admixture between selected and control bees and high FST SNPs
###

admix = read.table(file="SelvsConADMIXTURE.txt",header=T)
	#significantly shifted to be less M and C so either 1) A is more hygienic or 2) in this mixed background, A at these sites is more hygienic.
		#http://www.tandfonline.com/doi/pdf/10.1080/0005772X.1998.11099408	

		
pdf(file="HighFSTAdmix.pdf")	
#initiate plot, choose colors.		
par(mfrow=c(1,3),
oma = c(5,4,0,0) + 0.1,
mar = c(0,0,0,0) + 0.1)

cols = terrain.colors(2)

#plot C lineage 
boxplot(admix$C~admix$SelCon,
	ylim = c(0,1),
	names=FALSE,
	axes=FALSE,
	ylab="",
	boxwex = 0.75,
	staplewex = 0.5,
	whisklwd = 2,
	medlwd = 2,
	whisklty = "solid",
	staplewd = 2,
	outlwd = 1,
	boxlwd = 2,
	show.names=FALSE,
	frame=FALSE,
	col = c(cols),
	cex.lab=1.5, cex=1.5, cex.axis=1.5)
mtext("Proportion Admixed",side=2,line=2.4,cex=1.5)
mtext("C Lineage",side=1,line=1.5,cex=1.5)
means = aggregate(C~SelCon,FUN="mean",data=admix)
points(c(1,2),c(means[1,2],means[2,2]),pch=23,cex=2,lwd=2,bg="white")
axis(2, at=c(0,0.2,0.4,0.6,0.8,1), labels=c(0,0.2,0.4,0.6,0.8,1),cex.lab=1.5, cex=1.5, cex.axis=1.5)

#plot M lineage 
boxplot(admix$M~admix$SelCon,
	ylim = c(0,1),
	names=FALSE,
	axes=FALSE,
	ylab="",
	boxwex = 0.75,
	staplewex = 0.5,
	whisklwd = 2,
	medlwd = 2,
	whisklty = "solid",
	staplewd = 2,
	outlwd = 1,
	boxlwd = 2,
	show.names=FALSE,
	frame=FALSE,
	col = c(cols),
	cex.lab=1.5, cex=1.5, cex.axis=1.5)
means = aggregate(M~SelCon,FUN="mean",data=admix)
points(c(1,2),c(means[1,2],means[2,2]),pch=23,cex=2,lwd=2,bg="white")
mtext("M Lineage",side=1,line=1.5,cex=1.5)

#plot A lineage 
boxplot(admix$A~admix$SelCon,
	ylim = c(0,1),
	names=FALSE,
	axes=FALSE,
	ylab="",
	boxwex = 0.75,
	staplewex = 0.5,
	whisklwd = 2,
	medlwd = 2,
	whisklty = "solid",
	staplewd = 2,
	outlwd = 1,
	boxlwd = 2,
	show.names=FALSE,
	frame=FALSE,
	col = c(cols),
	cex.lab=1.5, cex=1.5, cex.axis=1.5)
means = aggregate(A~SelCon,FUN="mean",data=admix)
points(c(1,2),c(means[1,2],means[2,2]),pch=23,cex=2,lwd=2,bg="white")	
mtext("A Lineage",side=1,line=1.5,cex=1.5)

#plot legend
legend(0.5,0.5, 
	c("Control", "Selected"),
	fill =  c(cols),
	bty = "n",
	cex=2)
	
dev.off()