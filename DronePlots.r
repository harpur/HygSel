###
# Figures for Hygiene Data
###












#Load Functions and Packages -------------------------------
	#load(file="Hygiene.RData")
	#load(file="HYGRESULTS_FIGS.RDATA") #this is FSTP>3
library(ggplot2)
library(ggthemes) #https://github.com/jrnold/ggthemes
library(plyr)
library(reshape2)
library(RColorBrewer)
library(wesanderson)
library(grid)
library(gridExtra)
#http://www.ucl.ac.uk/~zctpep9/Archived%20webpages/Cookbook%20for%20R%20%C2%BB%20Colors%20(ggplot2).htm
#http://www.stat.ubc.ca/~jenny/STAT545A/block14_colors.html


#Plot Hygiene between selected and control populations ---------------
#build df
samps = read.table(file="/media/data1/forty3/drone/git/data/DroneSamps.txt")
samps = samps[which(samps$V3!="U"),]

ids = as.character(samps$V3)
ids[ids=="C"] = "Control"
ids[ids=="S"] = "Selected"
samps$V3 = ids


#plot 
samp.box<-ggplot(samps, aes(factor(V3), V4)) +  
	geom_violin(aes(fill=V3),colour = "white") +
	geom_boxplot(width=.2, lwd=1.1) +
	theme_classic() + 
	#scale_fill_brewer("Chevrolet") +
	scale_fill_manual(values=wes_palette(n=3, name="Chevalier")) +
	theme(
	axis.line.x = element_line(colour = "black", size=1.1),
    axis.line.y = element_line(colour = "black", size=1.1),
	strip.text.x = element_text(size = 14.5),
	axis.title.y = element_text(size = 14.5),
	axis.text.x = element_text(size=14),
	axis.text = element_text(size=14),
	legend.position="none"
	) +
	labs(
	x= "",
	y = "Hygienic Performance(% Cells Removed)") 

pdf(file="SelvsControboxplot.pdf")	
samp.box	
dev.off()
	
wilcox.test(samps$V4~samps$V3)
# 
#        Wilcoxon rank sum test with continuity correction
#
#data:  samps$V4 by samps$V3
#W = 19, p-value = 0.00001444
#alternative hypothesis: true location shift is not equal to 0
#
#Warning message:
#In wilcox.test.default(x = c(77.6, 58.59901089, 75.92039801, 59.38177445,  :
#  cannot compute exact p-value with ties
#

#Plot FST Results -------------------------------------

#Histogram
pdf(file="FSTHistogram.pdf")
options(scipen=5)
hist(as.numeric(All.Data$FST[as.numeric(All.Data$FST)>0]),
	breaks=60,
	main ="",
	col = terrain.colors(50),
	xlab = "Pairwise Fst",
	cex.lab = 1.2,
	axes = FALSE
	)
axis(1,at = seq(0,1,0.2),labels = TRUE,pos = 0, lwd=2, cex.axis=1.2)
axis(2,pos = 0, lwd=2, cex.axis=1.2)	
segments(x0 = 0.27, y0 = 0, x1 = 0.27, y1 = 100000,lty=2, lwd =2)
text(0.28,112500, "P<0.0001" )
dev.off()
		 

		 
		 
		 
		 
#Plot Creeper Results -------------------------------------

#create plot df
df1 = creeper[c(3,5,7)] #I'm going to plot -log10(pFST) from the dataframe PFST and colour significant Q values as red, or something.
df1 = df1[with(df1, order(start)),]
df1$chrom = factor(df1$chrom,  levels = unique(df1$chrom[order(as.numeric(df1$chrom))]))
		
#Create a function to generate a continuous color palette
rbPal = colorRampPalette(c("blue", "red"))

#This adds a column of color values
# based on the y values
df1$Col <- rbPal(30)[as.numeric(cut(df1$FSTP,breaks = 30))]

tiff(file="pFSTPlotHygiene.tiff",res=300,width=2000, height=3000)
#Main Plot
par(mfrow = c(8, 2),
	cex = 0.4, 
	mar = c(1, 2, 2, 2), 
	oma = c(1, 1, 1, 1)
	)
for(i in 1:16){

	plot(as.numeric(unlist(df1$start[df1$chrom==i])),as.numeric(unlist(df1$FSTP[df1$chrom==i])), 
		pch=19,
		cex=0.2,
		col = c(df1$Col[df1$chrom==i] ),
		axes=FALSE,
		xlab=NA,
		ylab=NA,
		xlim=c(0,  30000000),
		ylim=c(0,9)
	)
	xpos <- seq(0, max(as.numeric(unlist(df1$start[df1$chrom==i]))), by=1000000)
	axis(1, at=xpos,labels=xpos/1000)
	axis(2, pos=0, at=c(0,3,6,9), labels=c("","3","6","9"))
	lines(x = as.numeric(unlist(df1$start[df1$chrom==i])), 
		 y = rep(3, length(unlist(df1$start[df1$chrom==i]))) ,
		 col="red", 
		 lty=2)
	
	#Add QTL regions
	#Add in Oxley and Spotter  QTLs:
	if(i==1){
		rect(xleft = 3039231, ybottom=-0.05, xright=8453574, ytop=12, 
		border = NA,
		col=adjustcolor( "grey", alpha.f = 0.2) #Spotter
		)
		rect(xleft = 9418717, ybottom=-0.05, xright=16819942, ytop=12, 
		border = NA,
		col=adjustcolor("grey", alpha.f = 0.2) #Spotter
		)
	}else{
	if(i==2){
		rect(xleft = 14249515, ybottom=-0.05, xright=15407932, ytop=12, 
		border = NA,
		col=adjustcolor( "black", alpha.f = 0.2) #oxley
		)
		rect(xleft = 1, ybottom=-0.05, xright=12503099, ytop=12, 
		border = NA,
		col=adjustcolor("grey", alpha.f = 0.2) #Spotter
		)
	
	}else{
	if(i==5){
		rect(xleft = 9783962, ybottom=-0.05, xright=10814860, ytop=12, 
		border = NA,
		col=adjustcolor( "black", alpha.f = 0.2) #oxley
		)
	}else{
	if(i==6){
		rect(xleft = 11206828, ybottom=-0.05, xright=17739083, ytop=12, 
		border = NA,
		col=adjustcolor( "grey", alpha.f = 0.2) #Spotter
		)
	}else{
	if(i==7){
	rect(xleft = 9515998, ybottom=-0.05, xright=12848973, ytop=12, 
	border = NA,
	col=adjustcolor( "grey", alpha.f = 0.2) #Spotter
	)
	}else{
	if(i==12){
	rect(xleft = 1, ybottom=-0.05, xright=4003353, ytop=12, 
	border = NA,
	col=adjustcolor( "grey", alpha.f = 0.2) #Spotter
	)
	}else{
	if(i==13){
	rect(xleft =5247545, ybottom=-0.05, xright=10266737, ytop=12, 
	border = NA,
	col=adjustcolor( "grey", alpha.f = 0.2) #Spotter
	)
	}else{
	if(i==15){
	rect(xleft =1, ybottom=-0.05, xright=6643609, ytop=12, 
	border = NA,
	col=adjustcolor( "grey", alpha.f = 0.2) #Spotter
	)
	}else{
	
	if(i==16){
		rect(xleft = 2058269, ybottom=-0.05, xright=3127046, ytop=12, 
		border = NA,
		col=adjustcolor( "black", alpha.f = 0.2) #oxley
		)
		rect(xleft = 3196393, ybottom=-0.05, xright=6242592, ytop=12, 
		border = NA,
		col=adjustcolor( "grey", alpha.f = 0.2) #spotter
		)
	}else{		
	
}
}
}
}
}
}
}
}
}
mtext(i,side = 3, line = -1, adj = 0.05, cex = 1)
}


dev.off()

# Plot Tajima's D  -------------------


#build DF
Pi.Fst$grp = rep("High", nrow(Pi.Fst)); Pi.Fst = Pi.Fst[c(4, 10)]
Pi.genic$grp = rep("All", nrow(Pi.genic)); Pi.genic = Pi.genic[c(4, 14)]

whole = rbind(Pi.genic, Pi.Fst)
means.sem = ddply(whole , c("grp"), summarise,
	mean=mean(TajimaD,na.rm=T), 
	sem=sd(TajimaD,na.rm=T)/sqrt(length(TajimaD)))
means.sem = transform(means.sem, lower=mean-sem, upper=mean+sem)
means.sem$grp = c("All Windows", "High Fst Windows")


p<-ggplot(means.sem, aes(x = factor(grp), y = mean,fill = factor(grp))) +  
	geom_bar(stat="identity", width=0.5) + 
	#coord_cartesian(ylim = c(0, .7)) + 
	geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), position="dodge", width=.1, size = 1.1) 	+
	scale_fill_manual(values=wes_palette(n=2, name="Chevalier")) +
	theme_classic() + 
	coord_cartesian(ylim=c(0,0.35)) +
	theme(
		axis.line.x = element_line(size=1.1),
		axis.line.y = element_line(size=1.1),
		strip.text.x = element_text(size = 14.5),
		axis.title.y = element_text(size = 14.5),
		axis.text.x=element_text(size=14),
		axis.text=element_text(size=14),
		legend.position="none"	
		) +
	labs(
	x= "",
	y = "Tajima's D") +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_discrete(expand = c(0.1,0))	+
	geom_segment(aes(x = 1, y = .32, xend = 2, yend = .32), linetype=3) +
	geom_segment(aes(x = 1, y = .32, xend = 1, yend = .31), linetype=3) +
	geom_segment(aes(x = 2, y = .32, xend = 2, yend = .31), linetype=3) +
	annotate("text", x = 1.5,  y = 0.34,label = "***" )
	
	
pdf(file="TajimaDWindows.pdf", width = 5, height =10)	
p
dev.off()



# plot gamma -------------------------------


means.sem = ddply(gamma.res , c("qwd"), summarise,
	mean=mean(gamma,na.rm=T), 
	sem=sd(gamma,na.rm=T)/sqrt(length(gamma)))
means.sem = transform(means.sem, lower=mean-sem, upper=mean+sem)
means.sem$grp = c("Hygienic Candidates", "All Genes")


p<-ggplot(means.sem, aes(x = factor(grp), y = mean,fill = factor(grp))) +  
	geom_bar(stat="identity", width=0.5) + 
	#coord_cartesian(ylim = c(0, .7)) + 
	geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), position="dodge", width=.1, size = 1.1) 	+
	scale_fill_manual(values=wes_palette(n=2, name="Chevalier")) +
	theme_classic() + 
	coord_cartesian(ylim=c(0,0.8)) +
	theme(
		axis.line.x = element_line(size=1.1),
		axis.line.y = element_line(size=1.1),
		strip.text.x = element_text(size = 14.5),
		axis.title.y = element_text(size = 14.5),
		axis.text.x=element_text(size=14),
		axis.text=element_text(size=14),
		legend.position="none"	
		) +
	labs(
	x= "",
	y = expression(paste("Selection Coefficient", "(", gamma, ")", sep=""))) +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_discrete(expand = c(0.1,0))	+
	geom_segment(aes(x = 1, y = .65, xend = 2, yend = .65), linetype=3) +
	geom_segment(aes(x = 1, y = .65, xend = 1, yend = .64), linetype=3) +
	geom_segment(aes(x = 2, y = .65, xend = 2, yend = .64), linetype=3) +
	annotate("text", x = 1.5,  y = 0.66,label = "*****" ) +
	annotate("text", x = 1, y = 0.08, label = "9.2%", size = 5) +
	annotate("text", x = 2, y = 0.3, label = "18.9%", size = 5) 

	
pdf(file="Gamma.pdf", width = 5, height =10)	
p
dev.off()


#Plot Tajima's D in NA populations

#td.NA.re, the highest FSTP sites (>4)
means.sem = ddply(td.NA.res , c("grp"), summarise,
	mean=mean(TDC,na.rm=T), 
	sem=sd(TDC,na.rm=T)/sqrt(length(TDC)))
means.sem = transform(means.sem, lower=mean-sem, upper=mean+sem)
means.sem$grp = c( "All Genes","Hygienic Candidates")


p<-ggplot(means.sem, aes(x = factor(grp), y = mean,fill = factor(grp))) +  
	geom_bar(stat="identity", width=0.5) + 
	#coord_cartesian(ylim = c(0, .7)) + 
	geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), position="dodge", width=.1, size = 1.1) 	+
	geom_hline(yintercept=0,size=1.1 ) +
	scale_fill_manual(values=wes_palette(n=2, name="Chevalier")) +
	theme_classic() + 
	coord_cartesian(ylim=c(-1.20,0.5)) +
	theme(
		axis.line.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.line.y = element_line(size=1.1),
		strip.text.x = element_blank(),
		axis.title.y = element_text(size = 14.5),
		axis.text.x = element_blank(),
		axis.text=element_text(size=14),
		legend.position="none"	
		) +
	labs(
	x= "",
	y = "Tajima's D" ) +
	scale_y_continuous(expand = c(0,0)) +
	#scale_x_discrete(expand = c(0,0))	+
	geom_segment(aes(x = 1, y = -1.1, xend = 2, yend = -1.1), linetype=3) +
	geom_segment(aes(x = 1, y = -1.1, xend = 1, yend = -1.07), linetype=3) +
	geom_segment(aes(x = 2, y = -1.1, xend = 2, yend = -1.07), linetype=3) +
	annotate("text", x = 1.5,  y = -1.2,label = "***" ) +
	annotate("text", x = 1,  y = - 0.05,label = "All Genes", size =5 ) +
	annotate("text", x = 2,  y = 0.05,label = "Hygienic Candidates", size =5  ) 
	
	
	
pdf(file="NATajaimsD.pdf", width = 5, height =10)	
p
dev.off()







#Admixture plots --------------------------
#Load and prepare dataframes 
samps = read.table(file="/media/data1/forty3/drone/git/data/DroneSamps.txt"); names(samps)[2] = "ID"
nams = read.table(file="names",header=F); nams = nams[-(grep("^o",nams$V1)),]
whole = read.table(file="plink.3.Q.wholegenome",header=F)
reg = read.table(file="plink.3.Q.selectedsites",header=F)


whole$ID = reg$ID = nams
names(whole)[c(1:3)] = c("A","M","C") #manual, for now.
names(reg)[c(1:3)] = c("a.A","a.M","a.C") #manual, for now.
whole = merge(reg, whole, by = "ID")
whole = merge(whole, samps, by = "ID")


melted = melt(whole, id.vars=c("V3"),measure.vars=c("a.A", "a.M", "a.C", "A", "M", "C"))
BA = as.character(melted$variable)
BA[grep("[.]",BA)] = "Hygienic Loci"
BA[which(BA!="Hygienic Loci")] = "Genome Average"
melted$BA = BA
BA = as.character(melted$variable)
BA = gsub("a[.]","",BA)
melted$variable = BA

means.sem = ddply(melted, c("V3", "variable", "BA"), summarise,
	mean=mean(value,na.rm=T), 
	sem=sd(value,na.rm=T)/sqrt(length(value)))
means.sem = transform(means.sem, lower=mean-sem, upper=mean+sem)
means.sem = means.sem[means.sem$variable=="C",]
CS = as.character(means.sem$V3)
CS[CS=="C"] = "Control"
CS[CS=="S"] = "Selected"
means.sem$V3 = CS

		
	
dodge <- position_dodge(width=0.9)
#plot 
p<-ggplot(means.sem, aes(x = factor(V3), y = mean,fill = factor(BA))) +  
	geom_bar(position = dodge,stat="identity", ) + 
	#coord_cartesian(ylim = c(0, .7)) + 
	geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), position=dodge, width=.1, size = 1.1) 	+
	scale_fill_manual(values=wes_palette(n=2, name="Chevalier")) +
	theme_classic() + 
	coord_cartesian(ylim=c(0,1.1)) +
	theme(
		axis.line.x = element_line(size=1.1),
		axis.ticks.x = element_blank(),
		axis.line.y = element_line(size=1.1),
		strip.text.x = element_blank(),
		axis.title.y = element_text(size = 14.5),
		axis.text.x = element_text(size = 14.5),
		axis.text=element_text(size=14),
		legend.position=c(0.2,0.95),
		legend.title=element_blank()
		) +
	labs(
	x= "",
	y = "C-lineage Ancestry" ) +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_discrete(expand = c(.1,0))	+


	geom_segment(aes(x = 1.75, y = 0.94, xend = 2.25, yend = 0.94), linetype=3) +
	geom_segment(aes(x = 1.75, y = 0.94, xend = 1.75, yend = 0.92), linetype=3) +
	geom_segment(aes(x = 2.25, y = 0.94, xend = 2.25, yend = 0.92), linetype=3) +
	annotate("text", x = 2.0,  y = 0.95,label = "***" ) +
	
	geom_segment(aes(x = 0.75, y = 0.90, xend = 1.25, yend = 0.90), linetype=3) +
	geom_segment(aes(x = 0.75, y = 0.90, xend = 0.75, yend = 0.88), linetype=3) +
	geom_segment(aes(x = 1.25, y = 0.90, xend = 1.25, yend = 0.88), linetype=3) +
	annotate("text", x = 1,  y = 0.91,label = "*****" ) +
			
	geom_segment(aes(x = .75, y = 0.98, xend = 1.75, yend = 0.98), linetype=3) +
	geom_segment(aes(x = .75, y = 0.98, xend = .75, yend = 0.96), linetype=3) +
	geom_segment(aes(x = 1.75, y = 0.98, xend = 1.75, yend = 0.96), linetype=3) +
	annotate("text", x = 1.25,  y = 0.99,label = "***" ) +

	geom_segment(aes(x = 1.25, y = 1.04, xend = 2.25, yend = 1.04), linetype=3) +
	geom_segment(aes(x = 1.25, y = 1.04, xend = 1.25, yend = 1.02), linetype=3) +
	geom_segment(aes(x = 2.25, y = 1.04, xend = 2.25, yend = 1.02), linetype=3) +
	annotate("text", x = 1.75,  y = 1.05,label = "*****" )
	

pdf(file="ClineageSelected.pdf", width = 5, height =10)		
p
dev.off()





#Supplemental Admixture plots --------------------------
#Load and prepare dataframes 
	#ADD FST!
samps = read.table(file="/media/data1/forty3/drone/git/data/DroneSamps.txt"); names(samps)[2] = "ID"
nams = read.table(file="names",header=F); nams = nams[-(grep("^o",nams$V1)),]
whole = read.table(file="plink.3.Q.wholegenome",header=F)
reg = read.table(file="plink.3.Q.selectedsites",header=F)


whole$ID = reg$ID = nams
names(whole)[c(1:3)] = c("A","M","C") #manual, for now.
names(reg)[c(1:3)] = c("a.A","a.M","a.C") #manual, for now.
whole = merge(reg, whole, by = "ID")
whole = merge(whole, samps, by = "ID")


melted = melt(whole, id.vars=c("ID","V3"),measure.vars=c("a.A", "a.M", "a.C", "A", "M", "C"))
BA = as.character(melted$variable)
BA[grep("[.]",BA)] = "Selected"
BA[which(BA!="Selected")] = "Genome"
melted$BA = BA
BA = as.character(melted$variable)
BA = gsub("a[.]","",BA)
melted$variable = BA


melted.Before = melted[melted$BA!="Selected",]
names(melted.Before)[3] = "Lineage"

p<-ggplot(melted.Before, aes(x = ID, y = value, fill = Lineage)) +  
	geom_bar(stat = "identity") +
	#coord_cartesian(ylim = c(0, .7)) + 
		scale_fill_manual(values=wes_palette(n=3, name="FantasticFox")) +
	theme_classic() + 
	theme(
		axis.line.x = element_line(size=1.1),
		axis.ticks.x = element_blank(),
		axis.line.y = element_line(size=1.1),
		strip.text.x = element_blank(),
		axis.title.y = element_text(size = 14.5),
		axis.text.x = element_text(size = 14.5, angle=90),
		axis.text=element_text(size=14),
		legend.title=element_text(size=14)
		) +
	labs(
	x= "",
	y = "Proportion Ancestry" ) +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_discrete(expand = c(0,0)) 

pdf(file="AdmixturePlot.pdf", width = 10, height =5)		
p
dev.off()	
	
	
	
	
	
	
###B/A Comparison -------------------	

names(melted)[3] = "Lineage"	

	p<-ggplot(melted, aes(x = ID, y = value, fill = Lineage)) +  
	geom_bar(stat = "identity") +
	#facet_wrap(~.BA) +
	#coord_cartesian(ylim = c(0, .7)) + 
		scale_fill_manual(values=wes_palette(n=3, name="FantasticFox")) +
	theme_classic() + 
	theme(
		axis.line.x = element_line(size=1.1),
		axis.ticks.x = element_blank(),
		axis.line.y = element_line(size=1.1),
		strip.text.x = element_blank(),
		axis.title.y = element_text(size = 14.5),
		axis.text.x = element_text(size = 14.5, angle=90),
		axis.text=element_text(size=14),
		legend.title=element_text(size=14)
		) +
	labs(
	x= "",
	y = "Proportion Ancestry" ) +
	scale_y_continuous(expand = c(0,0)) 	
	
##############################	
	
		
		

		

#Plot M-allele purge (via Haploid data) ----- 
	#uses M.def, from DroneHapPlot.r


means.sem = M.def
means.sem$x = as.numeric(means.sem$x )

means.sem = ddply(means.sem , c("chr"), summarise,
	mean=mean(x,na.rm=T), 
	sem=sd(x,na.rm=T)/sqrt(length(x)))
means.sem = transform(means.sem, lower=mean-sem, upper=mean+sem)
means.sem$grp = c("Genome Average", "Hygienic Loci")


p<-ggplot(means.sem, aes(x = factor((grp)), y = mean,)) +  
	geom_bar(position = "dodge",stat="identity",color = "black") + 
	#coord_cartesian(ylim = c(0, .7)) + 
	geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), position="dodge", width=.1) 	+
	theme_few() + 
	theme(
		strip.text.x = element_text(size = 12.5),
		axis.title.y = element_text(size = 12.5),
		axis.text.x=element_text(size=12),
		axis.text=element_text(size=11),
		axis.ticks = element_blank(),
		legend.position="none"	
		) +
	labs(
	x= "",
	y = "M allele deficit (Selected versus Control)") +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_discrete(expand = c(0.1,0))	+
	geom_segment(aes(x = 1, y = .6, xend = 2, yend = .6), linetype=3) +
	geom_segment(aes(x = 1, y = .6, xend = 1, yend = .59), linetype=3) +
	geom_segment(aes(x = 2, y = .6, xend = 2, yend = .59), linetype=3) +
	annotate("text", x = 1.5,  y = 0.62,label = "***" )
	
p



















#Plot Lineage Purge Histogram ------------------------------
	#this is a plot of C-lineage ancestry in candidate regions of hygienic-selected vs unselected populations 
	#x11();hist(AMC.high.fst$C.SEL);x11();hist(AMC.high.fst$C.CON)
 
#Build plot data set:
sel = rep("Selected", length(AMC.high.fst$C.SEL))
cont = rep("Control", length(AMC.high.fst$C.CON))
df1 = data.frame(FST = c(AMC.high.fst$C.SEL,AMC.high.fst$C.CON))
df1$pop = as.factor(c(sel, cont))
df1 = df1[which(complete.cases(df1)),]
df1 = df1[which(!is.infinite(df1$FST)),]



#create barplot df
melted = melt(df1, id.vars=c("pop"),measure.vars=c("FST"))
means.sem = ddply(melted, c("pop"), summarise,
	mean=mean(value,na.rm=T), 
	sem=sd(value,na.rm=T)/sqrt(length(value)))
means.sem = transform(means.sem, lower=mean-sem, upper=mean+sem)

#plot 
p<-ggplot(means.sem, aes(x = factor(pop), y = mean,  fill = factor(pop))) +  
	geom_bar(position = "dodge",stat="identity",color = "black") + 
	coord_cartesian(ylim = c(0, .7)) + 
	geom_errorbar(aes(ymin=lower, ymax=upper), position="dodge", width=.1) 	+
	theme_few() + 
	scale_fill_brewer(palette = "Pastel1") +
	theme(
		strip.text.x = element_text(size = 12.5),
		axis.title.y = element_text(size = 12.5),
		axis.text.x=element_text(size=12),
		axis.text=element_text(size=11),
		axis.ticks = element_blank(),
		legend.position="none"	
		) +
	labs(
	x= "",
	y = "C vs M Lineage Introgression") +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_discrete(expand = c(0.1,0))	+
	geom_segment(aes(x = 1, y = .6, xend = 2, yend = .6), linetype=3) +
	geom_segment(aes(x = 1, y = .6, xend = 1, yend = .59), linetype=3) +
	geom_segment(aes(x = 2, y = .6, xend = 2, yend = .59), linetype=3) +
	annotate("text", x = 1.5,  y = 0.62,label = "***" )
	
p




				#p<-ggplot(df1, aes(x = pop, y = FST))+
				#	geom_violin(aes(fill = pop), adjust = 1.5, width = 1.5) +
				#	coord_flip() +
				#	theme_few()  + 
				#	theme(
				#		strip.text.y = element_text(size = 12),
				#		axis.title.x = element_text(size = 12),
				#		axis.text=element_text(size=10),
				#		legend.position="none"		
				#		) +
				#	labs(
				#		x= "",
				#		y = "C-lineage Introgression")
				#	
				#p




# Supplemental Figure of Ontario Hygiene ----------------------------
 hyg.ontario = read.table(file="clipboard",header=T) #it's the ontario data set
 
#Histogram
pdf(file="OntarioHygienicExpression.pdf")
options(scipen=5)
hist(as.numeric( hyg.ontario$V1),
	breaks=30,
	main ="",
	col = terrain.colors(50),
	xlab = "Hygienic Performance(% Cells Removed)",
	cex.lab = 1.2,
	axes = FALSE

	)
axis(1,at = seq(0,100,20),labels = TRUE,pos = 0, lwd=2, cex.axis=1.2)
axis(2,pos = 10, lwd=2, cex.axis=1.2)
dev.off()
		 




# GWAS PLOT ---------------------------------
	#this uses Association.RData and see CorrelationTest.r for data setup
	

SNPs = names(allele)[c(2:4)]
names(allele)[c(2:4)] = c("SNP1", "SNP2", "SNP3")
	
means.sem1 = ddply(allele , c("SNP1"), summarise,
mean=mean(V4,na.rm=T), 
sem=sd(V4,na.rm=T)/sqrt(length(V4)))
means.sem1 = transform(means.sem, lower=mean-sem, upper=mean+sem)

means.sem2 = ddply(allele , c("SNP2"), summarise,
mean=mean(V4,na.rm=T), 
sem=sd(V4,na.rm=T)/sqrt(length(V4)))
means.sem2 = transform(means.sem2, lower=mean-sem, upper=mean+sem)

means.sem3 = ddply(allele , c("SNP3"), summarise,
mean=mean(V4,na.rm=T), 
sem=sd(V4,na.rm=T)/sqrt(length(V4)))
means.sem3 = transform(means.sem3, lower=mean-sem, upper=mean+sem)





p1<-ggplot(means.sem1, aes(x = factor(SNP1), y = mean,fill = factor(SNP1))) +  
	geom_bar(stat="identity", width=0.5) + 
	#coord_cartesian(ylim = c(0, .7)) + 
	geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), position="dodge", width=.1, size = 1.1) 	+
	scale_fill_manual(values=wes_palette(n=2, name="Chevalier")) +
	theme_classic() + 
	coord_cartesian(ylim=c(0,100)) +
	theme(
		axis.line.x = element_line(size=1.1),
		axis.line.y = element_line(size=1.1),
		strip.text.x = element_text(size = 14.5),
		axis.title.y = element_text(size = 14.5),
		axis.text.x=element_text(size=14),
		axis.text=element_text(size=14),
		legend.position="none"	
		) +
	labs(
	x= "",
	y = "Hygienic Performance(% Cells Removed)") +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_discrete(expand = c(0.1,0))	+
	geom_segment(aes(x = 1, y = 85, xend = 2, yend = 85), linetype=3) +
	geom_segment(aes(x = 1, y = 85, xend = 1, yend = 84), linetype=3) +
	geom_segment(aes(x = 2, y = 85, xend = 2, yend = 84), linetype=3) +
	annotate("text", x = 1.5,  y = 86,label = "*" )

p2<-ggplot(means.sem2, aes(x = factor(SNP2), y = mean,fill = factor(SNP2))) +  
	geom_bar(stat="identity", width=0.5) + 
	#coord_cartesian(ylim = c(0, .7)) + 
	geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), position="dodge", width=.1, size = 1.1) 	+
	scale_fill_manual(values=wes_palette(n=2, name="Chevalier")) +
	theme_classic() + 
	coord_cartesian(ylim=c(0,100)) +
	theme(
		axis.line.x = element_line(size=1.1),
		axis.line.y = element_line(size=1.1),
		strip.text.x = element_text(size = 14.5),
		axis.title.y = element_text(size = 14.5),
		axis.text.x=element_text(size=14),
		axis.text=element_text(size=14),
		legend.position="none"	
		) +
	labs(
	x= "",
	y = "Hygienic Performance(% Cells Removed)") +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_discrete(expand = c(0.1,0))	+
	geom_segment(aes(x = 1, y = 85, xend = 2, yend = 85), linetype=3) +
	geom_segment(aes(x = 1, y = 85, xend = 1, yend = 84), linetype=3) +
	geom_segment(aes(x = 2, y = 85, xend = 2, yend = 84), linetype=3) +
	annotate("text", x = 1.5,  y = 86,label = "*" )

p3<-ggplot(means.sem3, aes(x = factor(SNP3), y = mean,fill = factor(SNP3))) +  
	geom_bar(stat="identity", width=0.5) + 
	#coord_cartesian(ylim = c(0, .7)) + 
	geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), position="dodge", width=.1, size = 1.1) 	+
	scale_fill_manual(values=wes_palette(n=2, name="Chevalier")) +
	theme_classic() + 
	coord_cartesian(ylim=c(0,100)) +
	theme(
		axis.line.x = element_line(size=1.1),
		axis.line.y = element_line(size=1.1),
		strip.text.x = element_text(size = 14.5),
		axis.title.y = element_text(size = 14.5),
		axis.text.x=element_text(size=14),
		axis.text=element_text(size=14),
		legend.position="none"	
		) +
	labs(
	x= "",
	y = "Hygienic Performance(% Cells Removed)") +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_discrete(expand = c(0.1,0))	+
	geom_segment(aes(x = 1, y = 85, xend = 2, yend = 85), linetype=3) +
	geom_segment(aes(x = 1, y = 85, xend = 1, yend = 84), linetype=3) +
	geom_segment(aes(x = 2, y = 85, xend = 2, yend = 84), linetype=3) +
	annotate("text", x = 1.5,  y = 86,label = "*" )

	

#grid.arrange(p1, p2, p3 ncol = 3)
	

pdf(file="snp1.pdf", width = 5, height =10)	
p1
dev.off()

	
pdf(file="snp2.pdf", width = 5, height =10)	
p2
dev.off()
	
pdf(file="snp3.pdf", width = 5, height =10)	
p3
dev.off()	
	
	
	
	
	
	
	
	
	
	
	
	
	



