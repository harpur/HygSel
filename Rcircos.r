###
# R Circos test
###


# Information from the RCircos README -------------
#http://cran.r-project.org/web/packages/RCircos/vignettes/Using_RCircos.pdf
#data(UCSa.HG19.Human.CytoBandIdeogram);
#head(UCSa.HG19.Human.CytoBandIdeogram);

#Chromosome ChromStart ChromEnd Band Stain
#1 chr1 0 2300000 p36.33 gneg
#2 chr1 2300000 5400000 p36.32 gpos25
#3 chr1 5400000 7200000 p36.31 gneg
#4 chr1 7200000 9200000 p36.23 gpos25
#5 chr1 9200000 12700000 p36.22 gneg
#6 chr1 12700000 16200000 p36.21 gpos50
#


# Testing RCircos Commands on AMEL -------------

# Plot paramaters -------------
require("RCircos")
Amel.CytoBandIdeogram=read.table(file="/media/data1/forty3/brock/balancedSNPs/data/AMELRCircos.txt",header=T)
genes = read.table(file="/media/data1/forty3/brock/balancedSNPs/data/NCBIGenes.txt",header=F);names(genes)[2]="Gene"


# Load FST data for all populations -------------

load(file="pop2_vs_sel.RDATA")

for(i in 1:nrow(Fst)){
	pos = unlist(Fst[i,3])
	fst = unlist(Fst[i,2])
	chr = unlist(Fst[i,1])
	outs = cbind(rep(chr, length(fst)), pos, pos+10, fst)
	write.list(outs,file="SelVvsConFst.df",append=T)
	print(i)
}

ypi = read.table(file="SelVvsConFst.df")
names(ypi) = c("chromosome", "start", "stop", "seg.mean")
ypi$seg.mean = (((ypi$seg.mean- min(ypi$seg.mean)) * (9 - 1)) / (max(ypi$seg.mean) - min(ypi$seg.mean))) + 1





# Map Heat Map of SNPs associated -------------
	#top 90%
	#perm.assoc
	#hap.assoc - not included here
	#cmh.assoc
	rp.assoc = read.table(file="RPARTAssocSNPs");names(rp.assoc)[2]="SNP"
	
perm.assoc = perm.assoc[c(2,1,6)];names(perm.assoc)=c("chr","snp","bp")
cmh.assoc = cmh.assoc[c(1,2,3)]	;names(cmh.assoc)=c("chr","snp","bp")
rp.assoc = rp.assoc[c(1,2,4)];	names(rp.assoc)=c("chr","snp","bp")
	
test = rbind(perm.assoc, rp.assoc, cmh.assoc)
test2=aggregate(test$snp, by=list(test$snp), length)	
	
names(test2) = c("snp", "num")	
test = merge(test, test2, by="snp")
test$stop = as.numeric(test$bp) + 100	
test = test[c(2,3,5,1,4)]	
names(test)=c("Chromosome", "chromStart", "chromEnd", "GeneName", "count")
hd = test




 
 
# Plot -------------


#Essential Commands and set yp
chr.exclude <- NULL;
cyto.info <- Amel.CytoBandIdeogram;
tracks.inside <- 10;
tracks.outside <- 0;
RCircos.Set.Core.Components(cyto.info, chr.exclude,
+ tracks.inside, tracks.outside);

params <- RCircos.Get.Plot.Parameters();
params$heatmap.width = 400
params$base.per.unit = 300
params$point.size = 0.5
RCircos.Reset.Plot.Parameters(params)

#Output information


pdf(file="out.pdf", height=8, width=8, compress=TRUE);
RCircos.Set.Plot.Area()
par(mai=c(0.25, 0.25, 0.25, 0.25));
plot.new();
plot.window(c(-2.5,2.5), c(-2.5, 2.5))
RCircos.Chromosome.Ideogram.Plot();



#Plotting of Tajima's D tracks
#Y TD
data.col <- 4;
track.num <- 5;
side <- "in";
by.fold <- quantile(ypi$seg.mean,0.95);
RCircos.Scatter.Plot(ypi, data.col,
+ track.num, side, by.fold);


#Plotting of High diversity genes
data.col <- 5;
track.num <- 7;
side <- "in";
RCircos.Heatmap.Plot(hd, data.col, track.num, side);

dev.off()


























#################################################

##plot SFS genes within each population
genes = read.table(file="NCBIGenes.txt",header=F);names(genes)[2]="Gene"
c.sfs =  read.table(file="AExonic.recode.vcf.SIGgeneSFS")
names(c.sfs)[1]="Gene"
c.sfs = merge(genes, c.sfs, by="Gene")
c.sfs$value = rep(3, nrow(c.sfs))
c.sfs = c.sfs[c(2,3,4,8,1)]; names (c.sfs)=c("chromosome", "start", "stop", "seg.mean")
c.sfs = c.sfs[with(c.sfs, order(start)),]
#c.sfs = (c.sfs[c.sfs$chromosome==2,])


a.sfs =  read.table(file="AExonic.recode.vcf.SIGgeneSFS")
names(a.sfs)[1]="Gene"
a.sfs = merge(genes, a.sfs, by="Gene")
a.sfs$value = rep(3, nrow(a.sfs))
a.sfs = a.sfs[c(2,3,4,8,1)]; names (a.sfs)=c("chromosome", "start", "stop", "seg.mean")
a.sfs = a.sfs[with(a.sfs, order(start)),]
#a.sfs = (a.sfs[a.sfs$chromosome==2,])

m.sfs =  read.table(file="MExonic.recode.vcf.SIGgeneSFS")
names(m.sfs)[1]="Gene"
m.sfs = merge(genes, m.sfs, by="Gene")
m.sfs$value = rep(3, nrow(m.sfs))
m.sfs = m.sfs[c(2,3,4,8,1)]; names (m.sfs)=c("chromosome", "start", "stop", "seg.mean")
m.sfs = m.sfs[with(m.sfs, order(start)),]
#m.sfs = (m.sfs[m.sfs$chromosome==2,])



y.sfs =  read.table(file="AExonic.recode.vcf.SIGgeneSFS")
names(y.sfs)[1]="Gene"
y.sfs = merge(genes, y.sfs, by="Gene")
y.sfs$value = rep(3, nrow(y.sfs))
y.sfs = y.sfs[c(2,3,4,8,1)]; names (y.sfs)=c("chromosome", "start", "stop", "seg.mean")
y.sfs = y.sfs[with(y.sfs, order(start)),]
#y.sfs = (y.sfs[y.sfs$chromosome==2,])

test=rbind(c.sfs,a.sfs,m.sfs,y.sfs)
test=aggregate(test[,5], by=list(test[,5]),length)




#plotting
chr.exclude <- NULL;
cyto.info <- Amel.CytoBandIdeogram;
tracks.inside <- 10;
tracks.outside <- 0;
RCircos.Set.Core.Components(cyto.info, chr.exclude,
+ tracks.inside, tracks.outside);


pdf(file="out.pdf", height=8, width=8, compress=TRUE);
RCircos.Set.Plot.Area()
par(mai=c(0.25, 0.25, 0.25, 0.25));
plot.new();
plot.window(c(-2.5,2.5), c(-2.5, 2.5))
RCircos.Chromosome.Ideogram.Plot();



data.col <- 4;
track.num <- 2;
side <- "in";
by.fold <- 0;
RCircos.Scatter.Plot(a.sfs, data.col,
+ track.num, side, by.fold);


data.col <- 4;
track.num <- 4;
side <- "in";
by.fold <- 0;
RCircos.Scatter.Plot(m.sfs, data.col,
+ track.num, side, by.fold);


data.col <- 4;
track.num <- 6;
side <- "in";
by.fold <- 0;
RCircos.Scatter.Plot(c.sfs, data.col,
+ track.num, side, by.fold);


data.col <- 4;
track.num <- 8;
side <- "in";
by.fold <- 0;
RCircos.Scatter.Plot(y.sfs, data.col,
+ track.num, side, by.fold);


dev.off()






#I want to plot diversity of all 4 populations
	#I'll start with showing bars where high diversity genes are 
a.div = read.table(file="a.TD.high",header=F)
a.div$V3 = a.div$V2 + 999
a.div$V6 = NULL ; a.div$V5 = NULL
names(a.div) = c("chromosome","start", "stop", "seg.mean")
a.div=a.div[a.div$chromosome=="11",]
a.div$seg.mean=a.div$seg.mean*250

M.div = read.table(file="M.TD.high",header=F)
M.div$V3 = M.div$V2 + 999
M.div$V6 = NULL ; M.div$V5 = NULL
names(M.div) = c("chromosome","start", "stop", "seg.mean")
M.div=M.div[M.div$chromosome=="11",]
M.div$seg.mean=M.div$seg.mean*250

Y.div = read.table(file="Y.TD.high",header=F)
Y.div$V3 = Y.div$V2 + 999
Y.div$V6 = NULL ; Y.div$V5 = NULL
names(Y.div) = c("chromosome","start", "stop", "seg.mean")
Y.div=Y.div[Y.div$chromosome=="11",]
Y.div$seg.mean=Y.div$seg.mean*250

A.div = read.table(file="A.TD.high",header=F)
A.div$V3 = A.div$V2 + 999
A.div$V6 = NULL ; A.div$V5 = NULL
names(A.div) = c("chromosome","start", "stop", "seg.mean")
A.div=A.div[A.div$chromosome=="11",]
A.div$seg.mean=A.div$seg.mean*250

load(file="HighDiv999.RData")
output =  output[output$V1=="11",]
output$V3=output$V2 + 999
output=output[c(1,2,3,4)]
names(output) = c("chromosome","start", "stop", "seg.mean")
output$seg.mean=output$seg.mean*250

#plotting
pdf(file="out.pdf", height=8, width=8, compress=TRUE);
RCircos.Set.Plot.Area()
par(mai=c(0.25, 0.25, 0.25, 0.25));
plot.new();
plot.window(c(-2.5,2.5), c(-2.5, 2.5))
RCircos.Chromosome.Ideogram.Plot();

#A
data.col <- 4;
track.num <- 2;
side <- "in";
by.fold <- 0;
RCircos.Scatter.Plot(A.div, data.col,
+ track.num, side, by.fold);


#M
data.col <- 4;
track.num <- 4;
side <- "in";
by.fold <- 0;
RCircos.Scatter.Plot(M.div, data.col,
+ track.num, side, by.fold);


#C
data.col <- 4;
track.num <- 6;
side <- "in";
by.fold <- 0;
RCircos.Scatter.Plot(a.div, data.col,
+ track.num, side, by.fold);


#Y
data.col <- 4;
track.num <- 8;
side <- "in";
by.fold <- 0;
RCircos.Scatter.Plot(Y.div, data.col,
+ track.num, side, by.fold);

#overlap
data.col <- 4;
track.num <- 10;
side <- "in";
by.fold <- 0;
RCircos.Scatter.Plot(output, data.col,
+ track.num, side, by.fold);



dev.off()
















