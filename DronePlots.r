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



