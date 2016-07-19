###
# Plot admixture results for North American Bees
###












#Load Functions and Packages -------------------------------
#load(file="Hygiene.RData")
library(ggplot2)
library(ggthemes) 
library(plyr)
library(reshape2)
library(RColorBrewer)
library(wesanderson)


#datasets -------------------------------
na.high = read.table(file="plink.3.Q.NAhigh",header=T); names(na.high) = paste("a",names(na.high),sep="." )
na.all = read.table(file="plink.3.Q.NAall",header=T)
na.high = na.high[c(1:6),]
na.all = na.all[c(1:6),]
na.all$ID = na.high$ID = seq(1,nrow(na.all), 1)


#stats-----------------------------------
wilcox.test(na.all$C, na.high$C)

#n307
#n308
#n315
#n318
#n320
#n322

        #Wilcoxon signed rank test

#data:  na.all$C
#V = 21, p-value = 0.03125
#alternative hypothesis: true location is not equal to 0




#Admixture plots --------------------------
#Load and prepare dataframes 
whole = merge(na.all, na.high, by = "ID")


melted = melt(whole, id.vars=c("ID"),measure.vars=c("a.M", "a.C", "A", "M", "C"))

means.sem = ddply(melted, c("variable"), summarise,
	mean=mean(value,na.rm=T), 
	sem=sd(value,na.rm=T)/sqrt(length(value)))
means.sem = transform(means.sem, lower=mean-sem, upper=mean+sem)
means.sem = means.sem[means.sem$variable %in% c("a.C", "C"),]
means.sem$variable = c("Selected", "Genome")	

	
dodge <- position_dodge(width=0.9)
#plot 
p<-ggplot(means.sem, aes(x = factor(variable), y = mean,fill = factor(variable))) +  
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
		legend.position = 'none',
		legend.title=element_blank()
		) +
	labs(
	x= "",
	y = "C-lineage Ancestry" ) +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_discrete(expand = c(.1,0))	+

	geom_segment(aes(x = 1, y = 0.96, xend = 2, yend = 0.96), linetype=3) +
	geom_segment(aes(x = 1, y = 0.96, xend = 1, yend = 0.94), linetype=3) +
	geom_segment(aes(x = 2, y = 0.96, xend = 2, yend = 0.94), linetype=3) +
	annotate("text", x = 1.5,  y = 0.97,label = "**" ) 

	

pdf(file="ClineageSelectedNA.pdf", width = 5, height =10)		
p
dev.off()