###
# Figures for Hygiene Data
###












#Load Functions and Packages -------------------------------
load(file="Hygiene.RData")
library(ggplot2)
library(ggthemes) #https://github.com/jrnold/ggthemes
library(plyr)
library(reshape2)
library(RColorBrewer)
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
	geom_violin(aes(fill=V3)) +
	geom_boxplot(width=.3) +
	theme_classic() + 
	scale_fill_brewer(palette = "Pastel1") +
	theme(
	axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
	strip.text.x = element_text(size = 12.5),
	axis.title.y = element_text(size = 12.5),
	axis.text.x = element_text(size=12),
	axis.text = element_text(size=11),
	legend.position="none"
	) +
	labs(
	x= "",
	y = "Hygienic (% Cells Removed)") 
samp.box	
	


#Plot FST Results -------------------------------------

#Histogram
pdf(file="FSTHistogram.pdf")
options(scipen=5)
hist(as.numeric(All.Data$FST[as.numeric(All.Data$FST)>0]),
	breaks=80,
	main ="",
	col = terrain.colors(70),
	xlab = "Pairwise Fst",
	lwd = 1.5
	)
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

tiff(file="pFSTPlotHygiene.tiff",res=300,width=1000, height=1500)
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



