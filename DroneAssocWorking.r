###
# Drone Selection Analyses, Working Files
###



#contains rough notes for the work performed on Drone analyses. Specifics can be found in similar files, mentioend throughout.
#Jumps between BASH, R, and other programs.



#Creation of VCF Files can be found in VCFCreation_DroneSelection.txt
	#The final VCF file is /media/data1/forty3/drone/vcf_drone/DroneSelectionFinal.recode.vcf
	







########################
#Fst Analyses
	#DroneFST.r

#VCFTools
	#Calculte FST  all selected VS control (2 vs 1+3=Selpop.txt)
vcftools --maf 0.05 --vcf DroneSelectionFinal.recode.vcf --weir-fst-pop pop2.txt --weir-fst-pop Selpop.txt --out pop2_vs_sel 

#found in cd /media/data1/forty3/drone/FST/SelvsCon




 
 
#LOAD into R 
  
 R
 
#FST Functions
source("/media/data1/forty3/brock/scripts/GOgetter.r")
MultiWayOverlapper = function(win.start,win.end,gene.start,gene.end,gene.list) {
  #this is a monster, but basically, looks within each row for genes  overlapping with whatever you want
  blah=outer(as.numeric(unlist(win.start)), as.numeric(unlist(gene.end)), "<=") 
  blah1=outer(as.numeric(unlist(win.end)), as.numeric(unlist(gene.start)), ">=") 
  blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
  if(!is.null(nrow(blah))){return(as.character(unlist(gene.list)[blah[,2]]))}	
  }

lenunique<-function(x){length(unique(x))}

write.list<-function(x,file){write.table(x,file,col.names=F,row.names=F,quote=F)}


#FST Munging:
#Fst in Pop2 vs 3 (control vs FAS) and then all selected vs control
fst23=read.table(file="/media/data1/forty3/drone/vcf_drone/pop2_vs_sel.weir.fst",header=T, colClasses=c("character","numeric","numeric"))
names(fst23)[3]="FST"
fst23=fst23[which(fst23$FST>=0 & !is.na(fst23$FST)),] 
pdf(file="/media/data1/forty3/drone/FST/HistogramofFASvsControlfst23.pdf")
hist(fst23$FST,col="lightgreen",main=NULL,xlab="fst23")
dev.off()
fst23$SNP=paste(fst23$CHROM,fst23$POS, sep=":")
fst23=fst23[,-c(1,2)]


map=read.table(file="/media/data1/forty3/drone/vcf_drone/DroneSelection.map");names(map)=c("CHROM", "SNP","NAL","POS");map$NAL=NULL
fst23=merge(fst23,map,by="SNP")
fst23=fst23[order(fst23$CHROM, fst23$POS), ]


source("/media/data1/forty3/brock/scripts/movingavg.r")
Fst=with(fst23, aggregate(FST,by=list(CHROM),function(x) movingAverage(x,n=1000,centered=F)))
FstPositions=with(fst23, aggregate(POS,by=list(CHROM),function(x) movingAverage(x,n=1000,centered=F)))
Fst$x=apply(Fst, 1, function(x) as.numeric(unlist(x[2]))[-c(1:100)])
FstPositions$x=apply(FstPositions, 1, function(x) as.numeric(unlist(x[2]))[-c(1:100)])
Fst$Pos=FstPositions$x;rm(FstPositions)
Fst=Fst[with(Fst, order(as.numeric(Group.1))), ]
highFstWindow=quantile(unlist(Fst$x),0.98)









#Creating FST Figures
#Figure of FST distribution
pdf(file = "/media/data1/forty3/drone/R/SelvsConFstPlot.pdf")
hist(fst23$FST,col = "lightgreen",xlab = "Fst",main = NULL, freq = F)
#hst=hist(fst23$FST,plot = F)
#axis(1, pos = hst$mids)#x axis 
#axis(4, pos = seq(1,413822,41382))#y axis 
dev.off()




#Figures and Tables of FST across all chromosomes
		#FST FIGURE:
		
pdf(file="FstChromosome1516.pdf")
par(mfrow=c(2, 1),mar = rep(4, 4))
for(i in c(1,5)){
	plot(as.numeric(unlist(Fst$Pos[Fst$Group.1==i])),as.numeric(unlist(Fst$x[Fst$Group.1==i])), 
	type="l",
	xlab="Genomic Position",
	ylab="Fst",
	ylim=c(0,0.4)
	)
	abline(h=c(median(unlist(Fst[1,2]))-1.5*IQR(unlist(Fst[1,2])),
		median(unlist(Fst[1,2]))+1.5*IQR(unlist(Fst[1,2]))),
		col="darkgrey",lty=2)
	abline(h=highFstWindow,col="red",lty=2)
	#abline(v=c(5484901,	5491534),lty=3,col="blue")
}
dev.off()


#Scribble plot:		
par(mfrow=c(8, 2),mar = rep(2, 4))
for(i in 1:16){
	plot(as.numeric(unlist(Fst$Pos[Fst$Group.1==i])),as.numeric(unlist(Fst$x[Fst$Group.1==i])), 
	type="l",
	axes=FALSE,
	xlab=NA,
	ylab=NA,
	ylim=c(0,0.4)
	)
	abline(h=highFstWindow95,col="red",lty=2)
}

		
#Table 1 Genes List		
	#Fst.genes=lapply(Fst.genes,unique)	
	tab = data.frame("GB" = unlist(Fst.genes))
	tab = GO.getter(tab)
	tab = merge(genes, tab,by="GB"); tab = tab[-5]
	write.list(tab, file="HighFSTGenes")	
		
			
		
		



HighFstPositions=as.matrix(apply(Fst, 1, function(x) unlist(as.numeric(unlist(x[3]))[which((as.numeric(unlist(x[2]))>=highFstWindow))])))
		#Fst$Group.1 for order of Chromosomes.

#Get out all SNPs within high FST Windows:

HighFstPositions=as.matrix(apply(Fst, 1, function(x) unlist(as.numeric(unlist(x[3]))[which((as.numeric(unlist(x[2]))>=highFstWindow))])))


#Get out all SNPs within high outlier windows:

HighSNPs=c()
for (i in c(1,3,5,6,7,9,11:16)){
	test=fst23[fst23$CHROM==i,]
	low=as.numeric(unlist(HighFstPositions[i,1]))-1000
	high=as.numeric(unlist(HighFstPositions[i,1]))+1000
	blah=outer(as.numeric(unlist(low)), as.numeric(unlist(test$POS)), "<=") 
	blah1=outer(as.numeric(unlist(high)), as.numeric(unlist(test$POS)), ">=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	blah=(test[blah[,2],])
	blah=blah[!duplicated(blah),]
	HighSNPs=rbind(HighSNPs,blah)
	print(i)
	#x11();plot(as.numeric(unlist(Fst$Pos[Fst$Group.1==i])),as.numeric(unlist(Fst$x[Fst$Group.1==i])), pch=19)
	#abline(v=blah$PO,col="red")
}


#Genes?
	#First, "Gene", not "CDS"

#Genes in each region of pi.high and TD.high will overlap with values in the "starts" and "ends"
gff=read.table(file="/media/data1/forty3/brock/scripts/NCBIGFF.txt",col.names=c("chrom","GB","st","en","type"))
genes=gff[which(gff$type=="gene"),]
sts=aggregate(genes$st,by=list(genes$chrom),function(x) x)
ends=aggregate(genes$en,by=list(genes$chrom),function(x) x)
scaff.genes=aggregate(as.character(genes$GB),by=list(genes$chrom),function(x) as.character(x))
#Merge everything together, and then find overlap of metrics (TD and PI) and genes.
#output,for each scaffold, the indices to get out genes.





Find.High.Regions=cbind(gene=scaff.genes,gstart=sts,gend=ends,HighFstPositions)
	#function(win.start,win.end,gene.start,gene.end,gene.list)
Fst.genes=apply(Find.High.Regions, 1, function(x) MultiWayOverlapper(x[7],x[7],x[4],x[6],x[2]) )
Fst.genes=lapply(Fst.genes,unique)


GenicSNPs=genes[genes$GB %in% unlist(Fst.genes),]
#write.table(GenicSNPs[,c(1,3,4,2)],file="GenicSNPsHighFst",col.names=F,row.names=F,quote=F)



		
		
		
		
		

SNPsinGenes=c()
for(i in 1:nrow(GenicSNPs)){
	test=cbind(rep(GenicSNPs[i,1],length(c(as.numeric(GenicSNPs[i,3]):as.numeric(GenicSNPs[i,4]))) ),c(as.numeric(GenicSNPs[i,3]):as.numeric(GenicSNPs[i,4])))
	print(i)
	SNPsinGenes=rbind(test,SNPsinGenes)
} 
 
#write.table(snps, file="GenicFstP3.snp", col.names=F,row.names=F,quote=F)







#I'd like to get out the highest FST SNPs from these genes.
chroms=intersect(fst23$CHROM, GenicSNPs$chrom)
GenicFST=c()
for(i in chroms){
	#i=3
	fst=fst23[fst23$CHROM==i,];genic=GenicSNPs[GenicSNPs$chrom==i,]
	blah1=outer(unlist(fst$POS),unlist(genic$st), ">=") 
	blah=outer(unlist(fst$POS), unlist(genic$en), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #row is df1, col is df2
	overlap=fst[(blah[,1]),];overlap$GB=genic[blah[,2],2]
	print(i)
	GenicFST=rbind(GenicFST,overlap)	
	} 
rm(fst,genic,blah,blah1,i)
	#That's all the SNPs in high FST genes.


	
HighestGenicFST=with(GenicFST, aggregate(FST, by=list(GB), max))
names(HighestGenicFST)[1]="GB"
test=merge(HighestGenicFST, GenicFST, by="GB")
test=test[test$x==test$FST,]
test$pos2=gsub(".*_","", test$SNP)
test$pos1=gsub("_.*","", test$SNP)
HighestGenicFST=test
	#205 Genes
pdf(file="/media/data1/forty3/drone/R/SelvsConFstPlot.pdf"); hist(fst23$FST,col="lightgreen",xlab="Fst",main=NULL);dev.off()

	
#took highest FST SNPs (>0.99) in that region and looked for associations
test=test[which(test$FST>quantile(fst23$FST,0.98)),]
table(test$CHROM)
write.table(test[,c(7)], file="GenicHighestFSTSNPs.snp", col.names=F,row.names=F, quote=F)















#Randomly chose one of the highest SNPs (if there are many)
x=test[unlist(tapply(1:nrow(test),test$GB,function(x) {if (length(x)>1) {sample(x,min(length(x),1))}else {x}})),]
x=x[complete.cases(x),]
write.table(test[,c(7,8)], file="GenicHighestFSTSNPs.snp", col.names=F,row.names=F, quote=F)
write.table(x[,c(7,8)], file="GenicHighestFSTSNPs.snp", col.names=F,row.names=F, quote=F)





##########
# Diversity Estimates
	#Using MY TD and Pi Scripts:
	#Run it within each population
	
vcftools --vcf DroneSelectionFinal.recode.vcf --keep pop2.txt --recode --out CON
vcftools --vcf DroneSelectionFinal.recode.vcf --keep pop3.txt --recode --out FAS 
	
	#FAS.recode.vcf and CON.recode.vcf can be scanned for Pi and TD
	#in /media/data1/forty3/drone/vcf_drone/FAS/split
R
load(file="/media/data1/forty3/drone/FST/SelvsCon/SelectedvsControlFST.RData") 
source("/media/data1/forty3/brock/scripts/VCFDiploidScans.r")
fils=list.files(pattern="*.vcf")
div=list.files(pattern="*.divers");div=gsub(".divers",".vcf",div)
fils=fils[!(fils %in% div)]

for(i in fils){
	VCF.TDThetaPiWindows(i, w=1000)
}
	
	


#convert to chromosome
divs=list.files(pattern="*.divers")

for(i in divs){

	div=read.table(file=i,header=T)
	#LocusID Scaffold Pos1 Post2*optional
	scaff=gsub(".divers","",paste("Group",i,sep=""))
	div$Scaff=rep(scaff,nrow(div))	
	div$locus=paste(div$Scaff, div$start,sep=":")
	write.table(div[c(8,7,1,2)],file="/media/data1/forty3/drone/ScaffConvert/convSEL",col.names=F,row.names=F,append=T,quote=F)
	write.table(div,file="convDiversSEL",col.names=F,row.names=F,append=T,quote=F)
	print(i)
	}
	
#/media/data1/forty3/drone/ScaffConvert/scaffold_to_chr.pl

select.div=read.table(file="/media/data1/forty3/drone/vcf_drone/FAS/split/convSEL_onChr.txt",header=F);names(select.div)=c("SNP","chr","st","en")
div=read.table(file="/media/data1/forty3/drone/vcf_drone/FAS/split/convDiversSEL");names(div)=c("start","end","pSNPs","Theta","Pi","TD","scaff","SNP")
select.div=merge(div,select.div,by="SNP")	
select.div=select.div[order(select.div$chr, select.div$start), ]

	
con.div=read.table(file="/media/data1/forty3/drone/vcf_drone/CONTROL/split/convCON_onChr.txt",header=F);names(con.div)=c("SNP","chr","st","en")
div=read.table(file="/media/data1/forty3/drone/vcf_drone/CONTROL/split/convDivers");names(div)=c("start","end","pSNPs","Theta","Pi","TD","scaff","SNP")
con.div=merge(div,con.div,by="SNP")	
con.div=con.div[order(con.div$chr, con.div$start), ]


	
#If selected, FST high windows should have low diversity in SEL and high diversity in CON

#e.g. chr 1: 19669445 19820509 and chr11: 10161028 10308643

#To get out the clumps of data:

matr=c()
for(i in 1:16){

	test=as.vector(unlist(HighFstPositions[i,1])) #Two regions in CHR3
	if(length(test)==0){print(i)}else{
		test2=c(test,test[length(test)]);test2=test2[-1]
		x=test2-test
	if(length(which(x>5000))<=1){
		matr=data.frame(rbind(matr,(cbind(i,min(test), max(test)))))
	}else{
		st=which(x>5000)[1]
		en=which(x>5000)[2]
		#First range (1:868)
		matr=data.frame(rbind(matr, cbind(i, min(test[1:st]), max(test[1:st]))))
		#write.table(matr,file="HighRangesbyCHR", col.names=F,row.names=F,append=T)
		#Second Range (869:1710)
		matr=data.frame(rbind(matr, cbind(i,min(test[st+2:en-1]),max(test[st+2:en-1]) )))
		#Third Range (1711:end)
		matr=data.frame(rbind(matr, cbind(i,min(test[en:length(test)]),max(test[en:length(test)]) )))
		
}
}
}

#i=14 did it.
names(matr)=c("chr","st","en");matr$chr=paste("chr", matr$chr,sep="");matr=matr[complete.cases(matr),]

con.div$con=rep("CON",nrow(con.div)); select.div$con=rep("SEL",nrow(select.div))
div=rbind(con.div,select.div)

#Now, I want windo by non-window ANOVA.

#I'd like to get out the highest FST SNPs from these genes.
div.window=c()
for(i in 1:nrow(matr)){
	#i=3
	chrom=matr[i,1]
	fst=matr[i,];genic=div[div$chr==chrom,]
	blah1=outer(unlist(fst$st),unlist(genic$en), "<=") 
	blah=outer(unlist(fst$en), unlist(genic$st), ">=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #row is df1, col is df2
	overlap=genic[blah[,2],]
	print(i)
	div.window=rbind(div.window,overlap)	
	} 
rm(fst,genic,blah,blah1,i)
div.window$wind=rep("TRUE",nrow(div.window))

div$wind=rep("FALSE",nrow(div))
test=div[!(div$SNP %in% div.window$SNP),]
div.window=rbind(test,div.window)
boxplot(div.window$Pi~div.window$con* div.window$win)

aggregate(div.window$Pi,by=list(div.window$con, div.window$win),mean,na.rm=T)
 # Group.1 Group.2           x
#1     CON   FALSE 0.003135762
#2     SEL   FALSE 0.003187745
#3     CON    TRUE 0.003467866
#4     SEL    TRUE 0.003150784

 # Group.1 Group.2            x
#1     CON   FALSE 6.040007e-06
#2     SEL   FALSE 6.067260e-06
#3     CON    TRUE 2.611912e-05
#4     SEL    TRUE 2.429507e-05



x=aov(Pi~con+wind+con*wind, data=div.window)
summary(x)
TukeyHSD(x)
boxplot(div.window$Pi~div.window$con+div.window$win)

test=div.window[div.window$wind=="FALSE",]
x=aov(test$Pi~test$con)
summary(x)

#               Df  Sum Sq   Mean Sq F value Pr(>F)
#test$con        1 0.00053 0.0005304   78.61 <2e-16 ***
#Residuals   21098 0.14234 0.0000067
#---
#Signif. codes:  0 â***â 0.001 â**â 0.01 â*â 0.05 â.â 0.1 â â 1
#110 observations deleted due to missingness
#> TukeyHSD(x)
#  Tukey multiple comparisons of means
#    95% family-wise confidence level
#
#Fit: aov(formula = test$Pi ~ test$con)
#
#$`test$con`
#                 diff           lwr           upr p adj
#SEL-CON -0.0003170825 -0.0003871804 -0.0002469846     0
#
#Control windwos have higher diversity than the same selected windows.


#Non-selected windows are higher diversity in SEL population. 
#Minutely different
	#Group.1           x
	#1     CON 0.003135762
	#2     SEL 0.003187745




test$wind=(test$st.x <= test$en.y & test$en.x >= test$st.y )
boxplot(test$Pi~test$con* test$win)
aov(test$Pi~test$con* test$win)




ch1=con.div[which(con.div$chr=="chr1"),]
x11();boxplot(ch1$Pi, ch1$Pi[ch1$st>19669445 & ch1$en<19820509])





sel1=select.div[which(select.div$chr=="chr1"),]
boxplot(ch1$Pi, ch1$Pi[ch1$st>19669445 & ch1$en<19820509],sel1$Pi, sel1$Pi[sel1$st>19669445 & sel1$en<19820509])



ch1=con.div[which(con.div$chr=="chr11"),]
sel1=select.div[which(select.div$chr=="chr11"),]
x11();boxplot(ch1$Pi, ch1$Pi[ch1$st>10161028 & ch1$en<10308643],sel1$Pi, sel1$Pi[sel1$st>10161028 & sel1$en<10308643])








	


source("/media/data1/forty3/brock/scripts/movingavg.r")	
x11();plot(movingAverage(con.div$start[con.div$chr=="chr7"],n=1000),movingAverage(con.div$Pi[con.div$chr=="chr7"],n=1000))
x11();plot(movingAverage(select.div$start[select.div$chr=="chr7"],n=100),movingAverage(select.div$Pi[select.div$chr=="chr7"],n=100))
	
	
	
	


#I'm curious if these regions are higehr diversity in C b/c of the high phenotypic diversity?
#nope.



	
	
	
	
	
	




##########
# LD Decay Around candidate SNPs
	#most work in DroneHAP.r
#I want to find haploblocks in each population to identify which regions where selected. 
	#Selected versus unselected haploblocks
#/media/data1/forty3/drone/HaploView
#PLINK:
plink --file DroneSelection --ld-window-kb 15 --noweb --allow-no-sex --keep controlBees.txt --out Control
plink --file DroneSelection --ld-window-kb 15 --noweb --allow-no-sex --remove controlBees.txt --out Selected
plink --file DroneSelection --ld-window-kb 15 --noweb --allow-no-sex --keep FASBees.txt --out FASSelected
	

plink --file DroneSelection --ld-window-kb 15 --noweb --allow-no-sex  --out ALL



	
R
sel = read.table(file = "FASSelected.ld",header = T) #15kb windows
con = read.table(file = "Control.ld",header = T)


#for CHR6, 1,1700000
con = con[con$CHR_A =="6",]
#con = con[con$BP_A <= 1700000,]
sel = sel[sel$CHR_A =="6",]
#sel = sel[sel$BP_A <= 1700000,]

dis = abs(sel$BP_A - sel$BP_B)#/1000
brks = hist(dis, breaks=10, plot=F)$breaks
grp = cut(dis, breaks=brks)
r2meansSEL = tapply(sel$R2, grp, FUN = mean)
x11();plot(r2meansSEL,
	ylab = "Average R2",
	xlab = "Binned Distance",
	pch = 19,
	ylim=c(0.1,1),
	xlim=c(0,10))
dis = abs(con$BP_A - con$BP_B)#/1000
brks = hist(dis, breaks=10, plot=F)$breaks
grp = cut(dis, breaks=brks)
r2meansCON = tapply(con$R2, grp, FUN = mean)
points(r2meansCON, col="blue", pch=19)

x11();boxplot(r2meansCON, r2meansSEL)


c(1000000, 1800000)

boxplot(con$R2[con$BP_A<=1000000 & con$BP_A<=1800000], sel$R2[sel$BP_A<=1000000 & sel$BP_A<=1800000])

boxplot(sel$R2, sel$R2[sel$BP_A<=1000000 & sel$BP_A<=1800000])

brks=hist(sel$BP_A,breaks=70,plot=F)$breaks
grp = cut(sel$BP_A, breaks=brks)

source("/media/data1/forty3/brock/scripts/movingavg.r")
selBP = aggregate(sel$R2, by=list(sel$BP_A), FUN = mean)
plot(movingAverage(sel$BP_A,n=1000),movingAverage(sel$R2,n=1000), pch=19)

points(movingAverage(con$BP_A,n=1000),movingAverage(con$R2,n=1000), pch=19)



#Did high R2 chromosomes have higher Fst?
	r2=aggregate(con$R2, by=list(con$CHR_A),mean)
	chrfst=aggregate(fst23$FST, by=list(fst23$CHROM),mean)
	 fst.r2.chrom=merge(r2,chrfst, by="Group.1");names(fst.r2.chrom)=c("chr","r2","fst")
#No relation between FST and R2 in control chromosomes.(t=-0.0143; P=0.988)

con$scaff=gsub("[:].*","",con$SNP_A)
fst23$scaff=gsub("[:].*","",fst23$SNP)
r2=aggregate(con$R2, by=list(con$scaff),mean)
chrfst=aggregate(fst23$FST, by=list(fst23$scaff),mean)
 fst.r2.scaff=merge(r2,chrfst, by="Group.1");names(fst.r2.scaff)=c("chr","r2","fst")

plot(fst.r2.scaff$r2, fst.r2.scaff$fst)
cor.test(fst.r2.scaff$r2[fst.r2.scaff$r2>0.6], fst.r2.scaff$fst[fst.r2.scaff$r2>0.6])
#No relation between FST and R2 in control scaffolds.(t=0.0356; P=0.9716)










####################################################
# I want to try a SET TEST
#######
#file="ControlHighFST.map"
#http://pngu.mgh.harvard.edu/~purcell/plink/anal.shtml#set

plink --file DroneSelection --extract CandidateSNPs.snp --out ControlHighFST --noweb --make-bed --keep controlBees.txt 
plink --bfile ControlHighFST  --recode --out ControlHighFST --noweb

plink --file ControlHighFST --set-test --set highFST.set --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --allow-no-sex --assoc --mperm 1000 --out CONTROLSET --set-p 0.01 --set-max 10 --qt-means


plink --file DroneSelection --extract CandidateSNPs.snp --out SelHighFST --noweb --make-bed --remove controlBees.txt 
plink --bfile SelHighFST  --recode --out SelHighFST --noweb

plink --file SelHighFST --set-test --set highFST.set --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --allow-no-sex --assoc --mperm 1000 --out SelSET --set-p 0.01 --set-max 10 --qt-means




##OK

R
load(file="/media/data1/forty3/drone/FST/SelvsCon/SelectedvsControlFST.RData")


map=read.table(file="ControlHighFST.map",header=F)
mns=read.table("CONTROLSET.qassoc.means",header=T)


CONset=read.table("CONTROLSET.qassoc.set.mperm",header=T)
CONset=unlist(strsplit(as.vector((CONset$SNPS)),"[|]"));CONset=CONset[which(CONset!="NA")]
CONpos=map[map$V2 %in% CONset,]
chroms=unique(CONpos$V1)
# 3  5  6  7  9 11 12 14 15 16

SELset=read.table("SelSET.qassoc.set.mperm",header=T)
SELset=unlist(strsplit(as.vector((SELset$SNPS)),"[|]"));SELset=SELset[which(SELset!="NA")]
SELpos=map[map$V2 %in% SELset,]




###
# Plotting these data:






par(mfrow=c(8, 2),mar = rep(2, 4))
for(i in 1:16){
	plot(as.numeric(unlist(Fst$Pos[Fst$Group.1==i])),as.numeric(unlist(Fst$x[Fst$Group.1==i])), 
	type="l",
	axes=FALSE,
	xlab=NA,
	ylab=NA,
	ylim=c(0,0.4)
	)
	abline(h=highFstWindow,col="red",lty=2)
}









#pdf(file="FstChromosome1516.pdf")
for(i in chroms){
	#pdf(file=paste("SignPointPlotFST", i, ".pdf",sep=""))
	plot(as.numeric(unlist(Fst$Pos[Fst$Group.1==i])),as.numeric(unlist(Fst$x[Fst$Group.1==i])), 
		type="l",
		bty="n",
		lwd=1.32,
		xaxt="n",
		xlab="Genomic Position (Kb)",
		ylab="Fixation Index (Fst)",
		ylim=c(0,0.35)
	)
	xpos <- seq(0, max(as.numeric(unlist(Fst$Pos[Fst$Group.1==i]))), by=1000000)
	axis(1, at=xpos,labels=xpos/1000)
	
	abline(h=c(median(unlist(Fst[1,2]))-1.5*IQR(unlist(Fst[1,2])),
		median(unlist(Fst[1,2]))+1.5*IQR(unlist(Fst[1,2]))),
		col="darkgrey",lty=2)
	abline(h=highFstWindow,col="tomato",lty=2)
	#abline(v=c(5484901,	5491534),lty=3,col="blue")
	
	#Add in Control-associated
	xpos=CONpos$V4[CONpos$V1==i]
	ypos=rep(0.3,length(xpos));ypos=ypos+(1:length(xpos)/100)
	points(xpos,ypos,cex=0.8,col="honeydew3",pch=19) #control 
	
	#Add in Selection-associated
	xpos=SELpos$V4[SELpos$V1==i]
	ypos=rep(0.25,length(xpos));ypos=ypos+(1:length(xpos)/100)
	points(xpos,ypos,cex=0.8,col="honeydew4",pch=17) #selected
	
	#dev.off()
}






##############
#Coverage in selected versus unselected sites:




































###############################################









#############
#From R, Use HighSNPs for association (18K)


#I've tried so many different iterations of this in order to find significance AND control for population stratification. 
#It's been a slog. To that end, I'm trying a new iteration I'm going to get positive associations in control bees 

plink --file /media/data1/forty3/drone/vcf_drone/DroneSelection --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --keep controlBees.txt --allow-no-sex --assoc --qt-means --adjust 


plink --file /media/data1/forty3/drone/vcf_drone/DroneSelection --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --remove controlBees.txt --allow-no-sex --assoc --qt-means --adjust --out SEL

cd /media/data1/forty3/drone/vcf_drone/
R
beta = read.table(file="/media/data1/forty3/drone/vcf_drone/plink.qassoc",header=T)
mns = read.table(file="/media/data1/forty3/drone/vcf_drone/plink.qassoc.means",header=T)
assoc = read.table(file="/media/data1/forty3/drone/vcf_drone/plink.qassoc.adjusted",header=T)
selmns = read.table(file="/media/data1/forty3/drone/vcf_drone/SEL.qassoc.means",header=T)
selassoc = read.table(file="/media/data1/forty3/drone/vcf_drone/SEL.qassoc.adjusted",header=T)

selassoc = selassoc[selassoc$UNADJ<0.05,]


beta=beta[which(beta$P<0.05),]
mns=mns[mns$SNP %in% beta$SNP,]
nrow(beta)


#58018 Most false positives.... (which is most likely, false P or false neg?)


#So, which are high FST?

load(file="/media/data1/forty3/drone/FST/SelvsCon/SelectedvsControlFST.RData")
#recall HighSNPs contains all snps in high FST windows.

nrow(HighSNPs[HighSNPs$SNP %in% beta$SNP,])/nrow(HighSNPs)
nrow(fst23[fst23$SNP %in% beta$SNP,])/nrow(fst23)


test=(HighSNPs[HighSNPs$SNP %in% beta$SNP,])
tight=(test[test$SNP %in% selassoc$SNP,])



test=(beta[beta$SNP %in% selassoc$SNP,])

nrow(HighSNPs[HighSNPs$SNP %in% test$SNP,])/nrow(HighSNPs)

#14 snps with high fst and associated in both pops.

          SNP         FST CHROM      POS
484901   3.16:831420 5.65076e-02     3 11702609
582523  5.14:1610874 5.77411e-02     5 10867292
582525  5.14:1610950 1.04598e-01     5 10867368
582526  5.14:1610955 6.11928e-02     5 10867373
582527  5.14:1610956 6.11928e-02     5 10867374
582582  5.14:1613370 1.86831e-02     5 10869788
582597  5.14:1614317 5.89269e-02     5 10870735
582600  5.14:1614671 6.19718e-03     5 10871089
680306   6.28:277943 1.67152e-01     6 13138943
714210   7.13:107873 1.89668e-02     7  4542740
714212   7.13:107930 1.89668e-02     7  4542797
717800  7.17:1046886 4.85678e-02     7  6415548
849797   9.12:774053 2.24441e-05     9 10023449
72254  11.18:1358314 2.74698e-02    11 10262934






#How many SNPs associated in control AND high FST across both?
nrow(HighSNPs[HighSNPs$SNP %in% beta$SNP,])
#201
test=HighSNPs[HighSNPs$SNP %in% beta$SNP,]

selmns[selmns$SNP %in% test$SNP,]




1.18:1377556
11.18:1377211







mns[mns$SNP=="1.32:17659",]
selmns[selmns$SNP=="1.32:17659",]










#I've tried so many different iterations of this in order to find significance AND control for population stratification. 
#It's been a slog. To that end, I'm trying a new iteration I'm going to get positive associations in control bees 

plink --file /media/data1/forty3/drone/vcf_drone/DroneSelection   --all-pheno --pheno /media/data1/forty3/drone/vcf_drone/HBexpression.txt --noweb --keep controlBees.txt --allow-no-sex --assoc --qt-means --adjust --out CONTROLALLPHENO
#Order: CoAT	CBP	P450	HYP	OBP18	OBP16	SPARC	TaTBP	VAMP


R
beta = read.table(file="CONTROLALLPHENO.P1.qassoc",header=T)
mns = read.table(file="CONTROLALLPHENO.P1.qassoc.means",header=T)
assoc = read.table(file="CONTROLALLPHENO.P1.qassoc.adjusted",header=T)


beta=beta[which(beta$P<0.05),]
assoc=assoc[assoc$FDR_BH<0.05,]

mns=mns[mns$SNP %in% beta$SNP,]

nrow(beta)




load(file="/media/data1/forty3/drone/FST/SelvsCon/SelectedvsControlFST.RData")
#recall HighSNPs contains all snps in high FST windows.

nrow(HighSNPs[HighSNPs$SNP %in% assoc$SNP,])/nrow(HighSNPs)
nrow(fst23[fst23$SNP %in% beta$SNP,])/nrow(fst23)









































































##############################################################


plink --file /media/data1/forty3/drone/vcf_drone/DroneSelection --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt  --noweb --allow-no-sex --assoc --adjust --thin 0.1 --maf 0.1 --out THIN01 --remove controlBees.txt



plink --file /media/data1/forty3/drone/vcf_drone/DroneSelection --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt  --noweb --allow-no-sex --linear --covar FSTPCA7.txt --adjust --out lin7ALL





#/media/data1/forty3/drone/FST/SelvsCon		
#This is my list of SNPs in high FST regions, Candidates for GWAS. 
write.list(HighSNPs$SNP, file="HighFSTRegions.snp")
	#I'm going to try a few things. 
		#First, just look for association with PLINK Assoc.
plink --file /media/data1/forty3/drone/vcf_drone/DroneSelection --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --extract HighFSTRegions.snp --noweb --allow-no-sex --assoc --adjust 
	source("/media/data1/forty3/brock/scripts/qq.r")
		x11();qq(snps$UNADJ)
			#WAY off expected.
	
	
	
	#Following EIG:
	snps=read.table(file="plink.qassoc.adjusted",header=T)
	snps=read.table(file="lin7ALL.assoc.linear.adjusted",header=T)
	
	#after a qassoc
	snps=read.table(file="lin7.qassoc",header=T)
	par(mfrow=c(2, 1),mar = rep(4, 4))
	for(i in c(6,5)){
	plot(as.numeric(unlist(snps$BP[snps$CHR==i])),as.numeric(-log10(snps$P[snps$CHR==i])),
	pch=19,
	ylim=c(0,8)
	)
}

	
	
	
	
	
	
	
	
	
	
	#Second, just look for association with PLINK Perm
	
	
plink --file /media/data1/forty3/drone/vcf_drone/DroneSelection --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --extract HighFSTRegions.snp --noweb --allow-no-sex --assoc --mperm 5000 --out PERM
	snps=read.table(file="PERM.qassoc.mperm",header=T)
	x11();qq(snps$EMP1)
	snps[snps$EMP2<0.05,]
	#      CHR           SNP  EMP1     EMP2
	#1921    5  5.14:1610939 2e-04 0.029390
	#11843   9   9.12:738270 2e-04 0.031990
	#12801   9   9.12:884562 2e-04 0.043790
	#12802   9   9.12:884614 2e-04 0.009998
	#13446  11 11.18:1340140 2e-04 0.024600
	#13447  11 11.18:1340142 2e-04 0.024600
	#16913  15 15.19:2092841 2e-04 0.016000
	#16914  15 15.19:2092921 2e-04 0.016000
	#16915  15 15.19:2093022 2e-04 0.016000
	#16916  15 15.19:2093043 2e-04 0.016000
	#16917  15 15.19:2093057 2e-04 0.032790
	#16918  15 15.19:2093102 2e-04 0.032790
	#16919  15 15.19:2093104 2e-04 0.032790
	#16920  15 15.19:2093147 2e-04 0.032790

plink --file /media/data1/forty3/drone/vcf_drone/DroneSelection --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt  --noweb --allow-no-sex --make-perm-pheno 1000 --out PHENPERM 
		
plink --file /media/data1/forty3/drone/vcf_drone/DroneSelection --pheno PHENPERM.pphe --extract HighFSTRegions.snp --noweb --allow-no-sex --assoc --out PHENPERM --all-pheno

cat PHENPERM.P* >> allPHEN.txt	
	

	
fils=list.files(pattern="PHENPERM.P[0-9]*")
norms=c()
for(i in fils){
	x = read.table(file=i,header=T)
	x = x[-c(4:8)]
	norms = rbind(x, norms)
	print(i)
}


#The significance will then be  the significance of getting the p value given distribution of p values
agg.norm=aggregate(norms$P, by=list(norms$SNP),min); names(agg.norm)[1]="SNP"
test=apply(snps[3],1,function(x) pnorm(x, mean(agg.norm$x), sd(agg.norm$x)))
snps=read.table(file="PHENPERM.qassoc",header=T)
	#I get no significant SNPs with this method (test>>0.05). 



	
	
	#Third, EIG time:
#From DroneEIG.txt:
#Jan2015 edit: 
#I also will try this with the set of High FST SNPs (300 high FST and 300 random selected)
vcftools --vcf /media/data1/forty3/drone/vcf_drone/DroneSelectionFinal.recode.vcf  --positions  HighFSTRegions.vcfsnp --plink


cp out.map /media/data1/forty3/drone/EIG5.0.2/bin
cp out.ped /media/data1/forty3/drone/EIG5.0.2/bin


#Moved that file to Eig/bin
cd /media/data1/forty3/drone/EIG5.0.2/bin

	
R
map=read.table(file="out.map",header=F)
map$V1=gsub("[.].*","",map$V2)
write.table(map,file="out.map",col.names=F,row.names=F,quote=F)
q()

./convertf -p parbahFST
#outputs drone.snp, drone.ind, and drone.eigenstratgeno
	#Update the .ind file to have the appropriate phenotypes (change ??? to the phenotype of interest)
	
#Made bahexampleFST.perl, input files 
	#
	#For full FST dataset, no 17s and all individuals, I found below:
	 #N    eigenvalue  difference    twstat      p-value effect. n
  # 1      8.559713          NA     9.166  7.48565e-10    22.589
  # 2      2.728468   -5.831245     1.977     0.010809    65.946
  # 3      2.461373   -0.267095     2.347    0.0057201    72.521
  # 4      2.146637   -0.314736     1.957    0.0111794    80.909
  # 5      2.004183   -0.142454     3.049    0.0015508    89.472
  # 6      1.773608   -0.230575     3.157   0.00125738   102.804
  # 7      1.495642   -0.277966     1.770    0.0152108   118.897
  # 8      1.274201   -0.221441    -0.066     0.180453   132.014
  # 9      1.202507   -0.071694     0.225     0.130763   140.055
  #10      1.100066   -0.102441    -0.271       0.2224   150.086



   
perl bahexample.perl
#Saved the  significant PCAs as FSTPCA7 and FSTPCA2 and will run PLINK linear with them as covariates.


#Now, linear time:

plink --file /media/data1/forty3/drone/vcf_drone/DroneSelection --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt  --noweb --allow-no-sex --linear --covar FSTPCA7.txt --adjust --out lin7ALL
	snps=read.table(file="lin7.assoc.linear.adjusted",header=T)
	x11();qq(snps$UNADJ)
	
	
	
	
	
	
	--extract HighFSTRegions.snp
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	plink --file DroneSelection --mpheno 1 --pheno DronePhenoHB.txt --extract GenicHighestFSTSNPsRAND.snp --noweb --allow-no-sex --linear --adjust --out noCoVar 
	
	UNsnps=read.table(file="noCoVar.assoc.linear.adjusted",header=T)
	x11();qq(UNsnps$UNADJ)
	
for K in {1..1000}
do
	shuf -n 10 GenicHighestFSTSNPsRAND.snp > output$K
	plink --file DroneSelection --mpheno 1 --pheno DronePhenoHB.txt --extract output$K --noweb --allow-no-sex --linear  --covar HBPCA.txt --adjust --out RAND$K
done


		
		
		
		





		
		
		
		
##
#To get gene soverlapping these regions, I'll load the GFF file.
#Genes in each region of pi.high and TD.high will overlap with values in the "starts" and "ends"
gff=read.table(file="/media/data1/forty3/brock/scripts/NCBIGFF.txt",col.names=c("chrom","GB","st","en","type"))




genes=gff[which(gff$type=="gene"),]
sts=aggregate(genes$st,by=list(genes$chrom),function(x) x)
ends=aggregate(genes$en,by=list(genes$chrom),function(x) x)
scaff.genes=aggregate(as.character(genes$GB),by=list(genes$chrom),function(x) as.character(x))
#Merge everything together, and then find overlap of metrics (TD and PI) and genes.
#output,for each scaffold, the indices to get out genes.





Find.High.Regions=cbind(gene=scaff.genes,gstart=sts,gend=ends,HighFstPositions)
	#function(win.start,win.end,gene.start,gene.end,gene.list)
Fst.genes=apply(Find.High.Regions, 1, function(x) MultiWayOverlapper(x[7],x[7],x[4],x[6],x[2]) )
Fst.genes=lapply(Fst.genes,unique)


GenicSNPs=genes[genes$GB %in% unlist(Fst.genes),]
#write.table(GenicSNPs[,c(1,3,4,2)],file="GenicSNPsHighFst",col.names=F,row.names=F,quote=F)



snps=c()
for(i in 1:nrow(GenicSNPs)){
	test=cbind(rep(GenicSNPs[i,1],length(c(as.numeric(GenicSNPs[i,3]):as.numeric(GenicSNPs[i,4]))) ),c(as.numeric(GenicSNPs[i,3]):as.numeric(GenicSNPs[i,4])))
	print(i)
	snps=rbind(test,snps)
} 
 
#write.table(snps, file="GenicFstP3.snp", col.names=F,row.names=F,quote=F)




#I'd like to get out the highest FST SNPs from these genes.
chroms=intersect(fst23$CHROM, GenicSNPs$chrom)
GenicFST=c()
for(i in chroms){
	#i=3
	fst=fst23[fst23$CHROM==i,];genic=GenicSNPs[GenicSNPs$chrom==i,]
	blah1=outer(unlist(fst$POS),unlist(genic$st), ">=") 
	blah=outer(unlist(fst$POS), unlist(genic$en), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #row is df1, col is df2
	overlap=fst[(blah[,1]),];overlap$GB=genic[blah[,2],2]
	print(i)
	GenicFST=rbind(GenicFST,overlap)	
	} 
rm(fst,genic,blah,blah1,i)



HighestGenicFST=with(GenicFST, aggregate(FST, by=list(GB), max))
names(HighestGenicFST)[1]="GB"
test=merge(HighestGenicFST, GenicFST, by="GB")
test=test[test$x==test$FST,]
test$pos2=gsub(".*_","", test$SNP)
test$pos1=gsub("_.*","", test$SNP)
HighestGenicFST=test
	#205 Genes
pdf(file="/media/data1/forty3/drone/R/SelvsConFstPlot.pdf"); hist(fst23$FST,col="lightgreen",xlab="Fst",main=NULL);dev.off()

	
#took highest FST SNPs (>0.99) in that region and looked for associations
test=test[which(test$FST>quantile(fst23$FST,0.98)),]
table(test$CHROM)
write.table(test[,c(7)], file="GenicHighestFSTSNPs.snp", col.names=F,row.names=F, quote=F)















#Randomly chose one of the highest SNPs (if there are many)
x=test[unlist(tapply(1:nrow(test),test$GB,function(x) {if (length(x)>1) {sample(x,min(length(x),1))}else {x}})),]
x=x[complete.cases(x),]
write.table(test[,c(7,8)], file="GenicHighestFSTSNPs.snp", col.names=F,row.names=F, quote=F)
write.table(x[,c(7,8)], file="GenicHighestFSTSNPs.snp", col.names=F,row.names=F, quote=F)






save.image(file="/media/data1/forty3/drone/R/SelectionScansPop3.RData")








#FST FIGURE:
par(mfrow=c(8, 2),mar = rep(2, 4))
for(i in 1:16){
plot(as.numeric(unlist(Fst$Pos[Fst$Group.1==i])),as.numeric(unlist(Fst$x[Fst$Group.1==i])), 
	type="l",
	axes=FALSE,
	xlab=NA,
	ylab=NA,
	ylim=c(0,0.4)
	)
	abline(h=quantile(unlist(Fst$x),0.98),col="red",lty=2)
}


#FST FIGURE:
pdf(file="FstOutlierExampleChromosome6.pdf")
par(mfrow=c(2, 1),mar = rep(4, 4))
for(i in c(1,6)){
plot(as.numeric(unlist(Fst$Pos[Fst$Group.1==i])),as.numeric(unlist(Fst$x[Fst$Group.1==i])), 
	type="l",
	xlab="Genomic Position",
	ylab="Fst",
	ylim=c(0,0.4)
	)
	abline(h=c(0.06,0.13),col="darkgrey",lty=2)
	abline(h=0.149,col="red",lty=2)
	#abline(v=c(1331658, 1337245),lty=3,col="blue")
}

#quantile(unlist(Fst[1,2]))

#> median(unlist(Fst[1,2]))+1.5*IQR(unlist(Fst[1,2]))
#[1] 0.12626
#> median(unlist(Fst[1,2]))-1.5*IQR(unlist(Fst[1,2]))
#[1] 0.05564069
dev.off()