###
# FST plot by Chromosomes
##








#I'm going to plot -log10(pFST) from the dataframe PFST and colour significant Q values as red, or something.
df.fst = PFST[with(PFST, order(POS)),]
df.fst$logP = -log10(df.fst$P) #ranges 0-12


#Create a function to generate a continuous color palette
rbPal = colorRampPalette(c("blue", "red"))

#This adds a column of color values
# based on the y values
df.fst$Col <- rbPal(30)[as.numeric(cut(df.fst$logP,breaks = 30))]

pdf(file="pFSTPlotHygiene.pdf")
#Main Plot
par(mfrow = c(8, 2),
	cex = 0.4, 
	mar = c(1, 2, 2, 2), 
	oma = c(1, 1, 1, 1)
	)
for(i in 1:16){

	plot(as.numeric(unlist(df.fst$POS[df.fst$CHROM==i])),as.numeric(unlist(df.fst$logP[df.fst$CHROM==i])), 
	pch=19,
	cex=0.2,
	col = c(df.fst$Col[df.fst$CHROM==i] ),
	axes=FALSE,
	xlab=NA,
	ylab=NA,
	xlim=c(0,  30000000),
	ylim=c(0,12)
	)
	xpos <- seq(0, max(as.numeric(unlist(df.fst$POS[df.fst$CHROM==i]))), by=1000000)
	axis(1, at=xpos,labels=xpos/1000)
	axis(2, pos=0, at=c(0,4,8,12), labels=c("","4","8","12"))
	
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

#Legend - Annoying, have to print it separately fro now.
pdf(file="legend.pdf")
legend_image <- as.raster(matrix(rbPal(30), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5)*13)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()




#Histogram
pdf(file="FSTHistogram.pdf")
options(scipen=5)
hist(df.fst$FST,
	breaks=80,
	main ="",
	col = terrain.colors(70),
	xlab = "Pairwise Fst",
	lwd = 1.5
	)
dev.off()




