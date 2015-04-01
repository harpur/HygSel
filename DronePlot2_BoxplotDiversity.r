############
# Drone Selection BoxPlot
############









#From "Diversity Estimates" in DroneAssocWorking.r



#Data:
head(div.window) 




x=aov(Pi~con+wind+con*wind, data=div.window)
summary(x)
#                Df Sum Sq   Mean Sq F value   Pr(F)
#con              1 0.0001 0.0000973   14.49 0.000141 ***
#wind             1 0.0004 0.0004331   64.50 9.67e-16 ***
#con:wind         1 0.0007 0.0006791  101.16  < 2e-16 ***
#Residuals   385163 2.5860 0.0000067
#---
#Signif. codes:  0 â***â 0.001 â**â 0.01 â*â 0.05 â.â 0.1 â â 1
#2353 observations deleted due to missingness
# TukeyHSD(x)
#  Tukey multiple comparisons of means
#    95% family-wise confidence level
#
#Fit: aov(formula = Pi ~ con + wind + con * wind, data = div.window)
#
#$con
#               diff          lwr          upr     p adj
#SEL-CON 3.17822e-05 1.541624e-05 4.814817e-05 0.0001411
#
#$wind
#                   diff          lwr          upr p adj
#TRUE-FALSE 0.0001473551 0.0001113943 0.0001833158     0
#
#$`con:wind`
#                             diff           lwr           upr     p adj
#SEL:FALSE-CON:FALSE  5.198391e-05  2.991927e-05  7.404855e-05 0.0000000
#CON:TRUE-CON:FALSE   3.321049e-04  2.654058e-04  3.988041e-04 0.0000000
#SEL:TRUE-CON:FALSE   1.502240e-05 -5.159910e-05  8.164391e-05 0.9383297
#CON:TRUE-SEL:FALSE   2.801210e-04  2.134225e-04  3.468195e-04 0.0000000
#SEL:TRUE-SEL:FALSE  -3.696151e-05 -1.035824e-04  2.965938e-05 0.4834280
#SEL:TRUE-CON:TRUE   -3.170825e-04 -4.087355e-04 -2.254295e-04 0.0000000

boxplot(div.window$Pi~div.window$con+div.window$win)
#I wonder if the CONSEL>CON is a result of sequencing oddities?
	#I'm going to restrict analyis only the chromosomes with SNPs with association.

		#CONset
		#seld=paste("Group", CONset,sep="")
		#seld=gsub(":.*","",seld)
		#test=(div.window[div.window$scaff %in% seld,])
		#
		#x11();boxplot(Pi~con+wind,data=test)
		#
		#
		#
		#x=aov(Pi~con+wind+con*wind, data=test)
		#summary(x)
		#

#Summary:
	#Candidate windows had less diversity in SEL than CON (Two-Way Anova F1,1=101.16, p< 2e-16 Tukey HSD P<0.0001)

#Same with K test
	#div.window$Pi=as.factor(div.window$Pi)
	#div.window$con=as.factor(div.window$con)
	#kruskal.test(Pi[wind=="TRUE"]~con[wind=="TRUE"], data=div.window)
	#	#TRUE v TRUE is different.
	#kruskal.test(Pi[wind=="TRUE"]~con[wind=="TRUE"], data=div.window)
	#kruskal.test(Pi[wind=="FALSE"]~con[wind=="FALSE"], data=div.window)


	
	
	
##Tajima's D
> x=aov(TD~con+wind+con*wind, data=div.window)
> summary(x)
                Df Sum Sq Mean Sq F value Pr(>F)
#con              1   1811  1810.7  1349.9 <2e-16 ***
#wind             1    232   232.3   173.2 <2e-16 ***
#con:wind         1    574   574.2   428.1 <2e-16 ***
#Residuals   376661 505243     1.3
#---
#Signif. codes:  0 â***â 0.001 â**â 0.01 â*â 0.05 â.â 0.1 â â 1
#10855 observations deleted due to missingness

TukeyHSD(x)
#  Tukey multiple comparisons of means
#    95% family-wise confidence level
#
#Fit: aov(formula = TD ~ con + wind + con * wind, data = div.window)
#
#$con
#             diff       lwr       upr p adj
#SEL-CON 0.1386695 0.1312722 0.1460668     0
#
#$wind
#               diff        lwr       upr p adj
#TRUE-FALSE 0.108707 0.09251551 0.1248985     0
#
#$`con:wind`
#                           diff         lwr         upr p adj
#SEL:FALSE-CON:FALSE  0.15753527  0.14755977  0.16751076 0e+00
#CON:TRUE-CON:FALSE   0.28009519  0.25004042  0.31014996 0e+00
#SEL:TRUE-CON:FALSE   0.09578193  0.06580622  0.12575765 0e+00
#CON:TRUE-SEL:FALSE   0.12255992  0.09250756  0.15261228 0e+00
#SEL:TRUE-SEL:FALSE  -0.06175333 -0.09172664 -0.03178003 7e-07
#SEL:TRUE-CON:TRUE   -0.18431325 -0.22557071 -0.14305580 0e+00
#
#	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
#################
###PLOTTING:

#I'd like a boxplot of CON vs SEL beside SEL.TRUE and CON.TRUE
#Pi
pdf("DiversityBoxplot.pdf")
par(mfrow = c(1, 2),
	mar = c(0, 0, 1, 1) + 0.1, 
	oma = c(5,4, 0, 0) +0.1
	)
#Left Panel, Sel vs Con:	
boxplot((div.window$Pi*10^3)~div.window$con,
	notch = T,
	col = c("gold","lightgreen"),
	ylab="",
	ylim=c(-0.1, 25),
	#names=c("Control", "Selected"),
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
	cex.lab=1.5, cex=1.5, cex.axis=1.5)
mtext(expression(paste("Nucleotide Diversity (Pi x", "10"^"3",")",sep=" ")),side=2,line=2.2,cex=1.5)
mtext("Genome Wide",side=1,line=1,cex=1.5)		
#legend("topright",inset=.0,title=NULL,c("Queen associated","Worker associated"),
#    fill=c("gold","lightgreen"),horiz=FALSE,cex=1.2)
means = aggregate((Pi*10^3)~con,FUN="mean",data=div.window)
points(c(1,2),c(means[1,2],means[2,2]),pch=23,cex=2,lwd=2,bg="white")
text(c(1,2),c(means[1,2]-0.8,means[2,2]-0.8),
    labels=format(c(means[1,2],means[2,2]),format="f",digits=3),
    pos=1,cex=1.2,col="black")	


boxplot((div.window$Pi[div.window$wind=="TRUE"]*10^3)~div.window$con[div.window$wind=="TRUE"],
	notch = T,
	col = (c("gold","lightgreen")),
	yaxt = "n",
	#names=c("Control", "Selected"),
	ylim=c(-0.1, 25),
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
	cex.lab=1.5, cex=1.5, cex.axis=1.5)
mtext("Within Candidate",side=1,line=1,cex=1.5)
mtext("Windows",side=1,line=2,cex=1.5)			
legend("topright",inset=.0,title=NULL,c("Control","Selected"),
    fill=c("gold","lightgreen"),horiz=FALSE,cex=1.2)
means = aggregate((Pi[wind=="TRUE"]*10^3)~con[wind=="TRUE"],FUN="mean",data=div.window)
points(c(1,2),c(means[1,2],means[2,2]),pch=22,cex=2,lwd=2,bg="white")
text(c(1,2),c(means[1,2]-0.8,means[2,2]-0.8),
    labels=format(c(means[1,2],means[2,2]),format="f",digits=3),
    pos=1,cex=1.2,col="black")	

dev.off()	
	
		

		
		
		
#Theta:
pdf("DiversityBoxplot.pdf")
par(mfrow = c(1, 2),
	mar = c(0, 0, 1, 1) + 0.1, 
	oma = c(5,4, 0, 0) +0.1
	)
#Left Panel, Sel vs Con:	
boxplot((div.window$Theta*10^3)~div.window$con,
	notch = T,
	col = c("gold","lightgreen"),
	ylab="",
	#ylim=c(-0.1, 25),
	#names=c("Control", "Selected"),
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
	cex.lab=1.5, cex=1.5, cex.axis=1.5)
mtext(expression(paste("Theta x", "10"^"3","",sep=" ")),side=2,line=2.2,cex=1.5)
mtext("Genome Wide",side=1,line=1,cex=1.5)		
#legend("topright",inset=.0,title=NULL,c("Queen associated","Worker associated"),
#    fill=c("gold","lightgreen"),horiz=FALSE,cex=1.2)
means = aggregate((Theta*10^3)~con,FUN="mean",data=div.window)
points(c(1,2),c(means[1,2],means[2,2]),pch=23,cex=2,lwd=2,bg="white")
text(c(1,2),c(means[1,2]-0.8,means[2,2]-0.8),
    labels=format(c(means[1,2],means[2,2]),format="f",digits=3),
    pos=1,cex=1.2,col="black")	


boxplot((div.window$Theta[div.window$wind=="TRUE"]*10^3)~div.window$con[div.window$wind=="TRUE"],
	notch = T,
	col = (c("gold","lightgreen")),
	yaxt = "n",
	#names=c("Control", "Selected"),
	ylim=c(-0.1, 25),
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
	cex.lab=1.5, cex=1.5, cex.axis=1.5)
mtext("Within Candidate",side=1,line=1,cex=1.5)
mtext("Windows",side=1,line=2,cex=1.5)			
legend("topright",inset=.0,title=NULL,c("Control","Selected"),
    fill=c("gold","lightgreen"),horiz=FALSE,cex=1.2)
means = aggregate((Theta[wind=="TRUE"]*10^3)~con[wind=="TRUE"],FUN="mean",data=div.window)
points(c(1,2),c(means[1,2],means[2,2]),pch=22,cex=2,lwd=2,bg="white")
text(c(1,2),c(means[1,2]-0.8,means[2,2]-0.8),
    labels=format(c(means[1,2],means[2,2]),format="f",digits=3),
    pos=1,cex=1.2,col="black")	

dev.off()	




#TD:
pdf("TajimaDBoxplot.pdf")
par(mfrow = c(1, 2),
	mar = c(0, 0, 1, 1) + 0.1, 
	oma = c(5,4, 0, 0) +0.1
	)
#Left Panel, Sel vs Con:	
boxplot((div.window$TD)~div.window$con,
	notch = T,
	col = c("gold","lightgreen"),
	ylab="",
	#ylim=c(-0.1, 25),
	#names=c("Control", "Selected"),
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
	cex.lab=1.5, cex=1.5, cex.axis=1.5)
mtext("Tajima's D",side=2,line=2.2,cex=1.5)
mtext("Genome Wide",side=1,line=1,cex=1.5)		
#legend("topright",inset=.0,title=NULL,c("Queen associated","Worker associated"),
#    fill=c("gold","lightgreen"),horiz=FALSE,cex=1.2)
means = aggregate((TD)~con,FUN="mean",data=div.window)
points(c(1,2),c(means[1,2],means[2,2]),pch=23,cex=2,lwd=2,bg="white")
text(c(1,2),c(means[1,2]-0.1,means[2,2]-0.1),
    labels=format(c(means[1,2],means[2,2]),format="f",digits=3),
    pos=1,cex=1,col="black")	


boxplot((div.window$TD[div.window$wind=="TRUE"])~div.window$con[div.window$wind=="TRUE"],
	notch = T,
	col = (c("gold","lightgreen")),
	#yaxt = "n",
	#names=c("Control", "Selected"),
	#ylim=c(-0.1, 25),
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
	cex.lab=1.5, cex=1.5, cex.axis=1.5)
mtext("Within Candidate",side=1,line=1,cex=1.5)
mtext("Windows",side=1,line=2,cex=1.5)			
legend("topright",inset=.0,title=NULL,c("Control","Selected"),
    fill=c("gold","lightgreen"),horiz=FALSE,cex=1.2)
means = aggregate((TD[wind=="TRUE"])~con[wind=="TRUE"],FUN="mean",data=div.window)
points(c(1,2),c(means[1,2],means[2,2]),pch=22,cex=2,lwd=2,bg="white")
text(c(1,2),c(means[1,2]-0.1,means[2,2]-0.1),
    labels=format(c(means[1,2],means[2,2]),format="f",digits=3),
    pos=1,cex=1,col="black")	

dev.off()	






