###
#
###


#Analysis of Fst following pFST (pFST.sh) and VCFFst (DroneFst.r)




# Load datasets and functions ---------------------------
source("/media/data1/forty3/brock/scripts/VarFunct.r")
source("/media/data1/forty3/brock/scripts/GOgetter.r")
load(file="/media/data1/forty3/drone/FST/SelvsCon/SelectedvsControlFSTPvals.RData")
gff = read.table(file = "/media/data1/forty3/brock/scripts/GFF/NCBIGFF.txt")
gff=gff[gff$V5=="gene",]
degs = read.table(file="/media/data1/forty3/drone/git/boutinDEGs.txt",header=T)

# Analysis ---------------------------
sigFST = PFST[PFST$q < 0.01,]
#sigFST1 = PFST[PFST$q < 0.001,]
#48933 SNPs significant at q<0.01
#1010 at Q<0.001 (sigFST1)

# Genes with high FST SNPs -------------------	
Genic.Fst1 = c()
chrom = intersect(gff$V1,sigFST$CHROM)
for(i in chrom){
	win.temp = gff[gff$V1==i,]
	deg.temp = sigFST[which(as.character(sigFST$CHROM)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$V3)), ">=") 
	blah1=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$V4)), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	Genic.Fst1 = rbind(temp,Genic.Fst1)
	print(i)
}
names(Genic.Fst1)[8] ="GB"
high.genes = as.character(unique(Genic.Fst1$GB))

# Write list for GO analysis -------------------	
FSTFBGN = GO.getter(Genic.Fst1)
FSTFBGN = FSTFBGN[!duplicated(FSTFBGN$FGN),]
FSTFBGN = FSTFBGN[!is.na(FSTFBGN$FGN),]
write.list(FSTFBGN$FGN, file="GOFst1")


#GOstats
#  GOBPID       Pvalue OddsRatio   ExpCount Count Size
#1  GO:0048468 0.0006252067  3.952355  4.4411348    12   31
#2  GO:0009653 0.0009395183  2.528406 10.6014184    21   74
#3  GO:0007610 0.0023563396  2.300438 11.3177305    21   79
#4  GO:0007417 0.0026979840  2.827385  6.4468085    14   45
#5  GO:0030707 0.0029029841       Inf  0.4297872     3    3
#6  GO:0032501 0.0039034210  1.530518 70.6283688    88  493
#7  GO:0016319 0.0039632986  3.708808  3.4382979     9   24
#8  GO:0003006 0.0044435018  7.639594  1.2893617     5    9
#9  GO:0048610 0.0044435018  7.639594  1.2893617     5    9
#10 GO:0007420 0.0054370643  3.474093  3.5815603     9   25
#11 GO:0007626 0.0078625235  6.106599  1.4326241     5   10
#12 GO:0035107 0.0081441255  3.516937  3.1517730     8   22
#13 GO:0048736 0.0081441255  3.516937  3.1517730     8   22
#14 GO:0048737 0.0081441255  3.516937  3.1517730     8   22
#15 GO:0007472 0.0081441255  3.516937  3.1517730     8   22
#16 GO:0007476 0.0081441255  3.516937  3.1517730     8   22
#17 GO:0035120 0.0081441255  3.516937  3.1517730     8   22
#18 GO:0035114 0.0081441255  3.516937  3.1517730     8   22
#19 GO:0030537 0.0100584246  8.114478  1.0028369     4    7
#20 GO:0008345 0.0100584246  8.114478  1.0028369     4    7
#21 GO:0048513 0.0109626026  1.787143 17.6212766    27  123
#22 GO:0032504 0.0114035867  2.480000  6.0170213    12   42
#23 GO:0048609 0.0114035867  2.480000  6.0170213    12   42
#24 GO:0022414 0.0114035867  2.480000  6.0170213    12   42
#25 GO:0000003 0.0114035867  2.480000  6.0170213    12   42
#26 GO:0007154 0.0140194764  1.804173 15.4723404    24  108
#27 GO:0023052 0.0140194764  1.804173 15.4723404    24  108
#28 GO:0048563 0.0144923668  3.072165  3.4382979     8   24
#29 GO:0048569 0.0144923668  3.072165  3.4382979     8   24
#30 GO:0007560 0.0144923668  3.072165  3.4382979     8   24
#31 GO:0050896 0.0158213290  1.593777 26.6468085    37  186
#32 GO:0031175 0.0178835953  6.080808  1.1460993     4    8
#33 GO:0048812 0.0178835953  6.080808  1.1460993     4    8
#34 GO:0007276 0.0203563078  2.341397  5.7304965    11   40
#35 GO:0019953 0.0203563078  2.341397  5.7304965    11   40
#36 GO:0008045 0.0204370083       Inf  0.2865248     2    2
#37 GO:0030703 0.0204370083       Inf  0.2865248     2    2
#38 GO:0007304 0.0204370083       Inf  0.2865248     2    2
#39 GO:0009887 0.0247950692  2.367788  5.1574468    10   36
#40 GO:0051705 0.0283487928  2.855043  3.1517730     7   22
#41 GO:0048707 0.0298314751  2.580575  3.8680851     8   27
#42 GO:0007552 0.0298314751  2.580575  3.8680851     8   27
#43 GO:0009886 0.0298314751  2.580575  3.8680851     8   27
#44 GO:0032502 0.0322091131  1.357814 65.0411348    77  454
#45 GO:0009605 0.0424844345  4.047138  1.4326241     4   10
#                                             Term
#1                                cell development
#2              anatomical structure morphogenesis
#3                                        behavior
#4              central nervous system development
#5               ovarian follicle cell development
#6                multicellular organismal process
#7                       mushroom body development
#8  developmental process involved in reproduction
#9       cellular process involved in reproduction
#10                              brain development
#11                            locomotory behavior
#12                        appendage morphogenesis
#13                          appendage development
#14    imaginal disc-derived appendage development
#15                        wing disc morphogenesis
#16       imaginal disc-derived wing morphogenesis
#17         post-embryonic appendage morphogenesis
#18  imaginal disc-derived appendage morphogenesis
#19                                larval behavior
#20                     larval locomotory behavior
#21                              organ development
#22            multicellular organism reproduction
#23  multicellular organismal reproductive process
#24                           reproductive process
#25                                   reproduction
#26                             cell communication
#27                                      signaling
#28             post-embryonic organ morphogenesis
#29               post-embryonic organ development
#30                    imaginal disc morphogenesis
#31                           response to stimulus
#32                  neuron projection development
#33                neuron projection morphogenesis
#34                              gamete generation
#35                            sexual reproduction
#36                            motor axon guidance
#37                             eggshell formation
#38          chorion-containing eggshell formation
#39                            organ morphogenesis
#40       behavioral interaction between organisms
#41           instar larval or pupal morphogenesis
#42                                  metamorphosis
#43                   post-embryonic morphogenesis
#44                          developmental process
#45                  response to external stimulus
#







# Do I have DEGs with high FST? -------------------	
#1) Boutin's DEGs (they used 11 169 total genes)
head(degs[degs$GB %in% high.genes,])
	#6/81 DEGS at 0.01 and lower
	#4/81 DEGS at 0.01 and lower with permutation (and all show DWN regulation)


#DEGS not enriched for significant SNPs. 
	#Expression at a time versus development....
 
#2) Foster's DEPs
head(deps[deps$GB %in% high.genes,])	
	# none. at 0.01 or lower


#3) No overlap between DEGs and DEPs.	
	
	
#Yes, DEGs have High FST.
	

	
#Permutation....AGAIN....	
	
# Analysis ---------------------------
sigFST = PFST[PFST$q < 0.01,]

#48933 SNPs significant at q<0.01
#1010 at Q<0.001 (sigFST1)



# Genes with high FST SNPs, Permutation tests -------------------	
	#Generate expected distribution to see if number of high FST snps per gene is expected by chance alone.
for(x in 1:100){
Genic.Fstq = c()
chrom = intersect(gff$V1,PFST$CHROM)
PFSTp = transform(PFST, q = sample(q))
for(i in chrom){
	win.temp = gff[gff$V1==i,]
	deg.temp = PFSTp[which(as.character(PFSTp$CHROM)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$V3)), ">=") 
	blah1=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$V4)), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	Genic.Fstq = rbind(temp,Genic.Fstq)
	print(c(i,x))
}
write.list(Genic.Fstq, file="Genic.Fstq",append=T)
}

Genic.Fstq =  fread("Genic.Fstq",header=FALSE,colClasses="character")
Genic.Fstq = data.frame(Genic.Fstq)
#made grp group ings grp[grp>1 & grp<265579]=1
Genic.Fstq$grp = grp
#grp = c(1:nrow(Genic.Fstq))
#grp[grp>=1 & grp<=265579]=1
#grp[grp>=265580 & grp<=531158]=2
Genic.Fstq$grp = grp
Genic.Fstq$V6 = as.numeric(Genic.Fstq$V6)
expect =  aggregate(Genic.Fstq$V6, by=list(Genic.Fstq$V8, Genic.Fstq$grp), function(x) length(x[x<0.01]))
#number of high FST SNPs within eacvh gene for each permutation
mn =  aggregate(expect$x, by = list(expect$Group.1), mean)
sd =  aggregate(expect$x, by = list(expect$Group.1), sd)
expec = merge(mn, sd, by="Group.1")





#Output observed.
obs = aggregate(Genic.Fst1$q, by=list(Genic.Fst1$GB), length ) 	
names(obs)[2] = "obs"	
test = merge(obs, expec, by = "Group.1")	

test$p = apply(test, 1, function(x) pnorm(as.numeric(x[2]),mean=as.numeric(x[3]),sd=as.numeric(x[4]),lower.tail=FALSE))	
names(test)[1]="GB"	
sigFSTperm = test #contains significant SNPs within gene regions.
#save.image(file="/media/data1/forty3/drone/FST/SelvsCon/SelectedvsControlFSTPvals.RData")	
	
# Write list for GO analysis -------------------	
FSTFBGN = GO.getter(test[test$p<0.05,])
FSTFBGN = FSTFBGN[!duplicated(FSTFBGN$FGN),]
FSTFBGN = FSTFBGN[!is.na(FSTFBGN$FGN),]
write.list(FSTFBGN$FGN, file="GOFst.perm")
	
	
	
	
#After permutation:
 (summary(outcomeBP))
#       GOBPID      Pvalue OddsRatio  ExpCount Count Size
#1  GO:0009605 0.002325024 10.313725 0.6312057     4   10
#2  GO:0008045 0.003942236       Inf 0.1262411     2    2
#3  GO:0030703 0.003942236       Inf 0.1262411     2    2
#4  GO:0007304 0.003942236       Inf 0.1262411     2    2
#5  GO:0031175 0.010802103  9.181395 0.5049645     3    8
#6  GO:0048812 0.010802103  9.181395 0.5049645     3    8
#7  GO:0048468 0.010903963  3.747470 1.9567376     6   31
#8  GO:0030707 0.011339527 30.344828 0.1893617     2    3
#9  GO:0034470 0.011339527 30.344828 0.1893617     2    3
#10 GO:0034660 0.011339527 30.344828 0.1893617     2    3
#11 GO:0016319 0.014620276  4.078947 1.5148936     5   24
#12 GO:0007420 0.017383564  3.872024 1.5780142     5   25
#13 GO:0030182 0.027713892  5.725291 0.6943262     3   11
#14 GO:0048666 0.027713892  5.725291 0.6943262     3   11
#15 GO:0006935 0.034770292 10.099617 0.3156028     2    5
#16 GO:0042330 0.034770292 10.099617 0.3156028     2    5
#17 GO:0040011 0.034770292 10.099617 0.3156028     2    5
#18 GO:0007411 0.034770292 10.099617 0.3156028     2    5
#19 GO:0007409 0.034770292 10.099617 0.3156028     2    5
#20 GO:0048699 0.035295719  5.085271 0.7574468     3   12
#21 GO:0009653 0.038679583  2.173846 4.6709220     9   74
#22 GO:0009628 0.043834224  4.573256 0.8205674     3   13
#23 GO:0035107 0.045212051  3.406536 1.3886525     4   22
#24 GO:0048736 0.045212051  3.406536 1.3886525     4   22
#25 GO:0048737 0.045212051  3.406536 1.3886525     4   22
#26 GO:0007472 0.045212051  3.406536 1.3886525     4   22
#27 GO:0007476 0.045212051  3.406536 1.3886525     4   22
#28 GO:0035120 0.045212051  3.406536 1.3886525     4   22
#29 GO:0035114 0.045212051  3.406536 1.3886525     4   22
#                                            Term
#1                  response to external stimulus
#2                            motor axon guidance
#3                             eggshell formation
#4          chorion-containing eggshell formation
#5                  neuron projection development
#6                neuron projection morphogenesis
#7                               cell development
#8              ovarian follicle cell development
#9                               ncRNA processing
#10                       ncRNA metabolic process
#11                     mushroom body development
#12                             brain development
#13                        neuron differentiation
#14                            neuron development
#15                                    chemotaxis
#16                                         taxis
#17                                    locomotion
#18                                 axon guidance
#19                                  axonogenesis
#20                         generation of neurons
#21            anatomical structure morphogenesis
#22                  response to abiotic stimulus
#23                       appendage morphogenesis
#24                         appendage development
#25   imaginal disc-derived appendage development
#26                       wing disc morphogenesis
#27      imaginal disc-derived wing morphogenesis
#28        post-embryonic appendage morphogenesis
#29 imaginal disc-derived appendage morphogenesis
	
# NSYN SNPs in my highpFST list? ------------------------
	#From SNPEFF.sh
	
nsyn = read.table(file="/media/data1/forty3/drone/vcf_drone/exons1.eff",header=F,colClasses = rep("character", 17))
nsyn$SNP = paste(nsyn$V1, nsyn$V2, sep=":")
nsyn = nsyn[c(18,3,4,8,11,13,14)]	
test = merge(sigFST, nsyn, by = "SNP")	
test[test$V13=="NON_SYNONYMOUS_CODING",]
	#494 coding SNPs are nonsynonymous
test$V8[test$V13=="NON_SYNONYMOUS_CODING"]
	
	
	
#GO with genes with NSYN outliers:
summary(outcomeBP)
#       GOBPID       Pvalue OddsRatio   ExpCount Count Size
#1  GO:0031175 0.0002604537 37.636364 0.14184397     3    8
#2  GO:0048812 0.0002604537 37.636364 0.14184397     3    8
#3  GO:0008045 0.0003020099       Inf 0.03546099     2    2
#4  GO:0030182 0.0007407623 23.471591 0.19503546     3   11
#5  GO:0048666 0.0007407623 23.471591 0.19503546     3   11
#6  GO:0032989 0.0008850661 11.800866 0.46099291     4   26
#7  GO:0048699 0.0009761224 20.848485 0.21276596     3   12
#8  GO:0048468 0.0017539635  9.580247 0.54964539     4   31
#9  GO:0006935 0.0029225832 40.057971 0.08865248     2    5
#10 GO:0042330 0.0029225832 40.057971 0.08865248     2    5
#11 GO:0040011 0.0029225832 40.057971 0.08865248     2    5
#12 GO:0007411 0.0029225832 40.057971 0.08865248     2    5
#13 GO:0007409 0.0029225832 40.057971 0.08865248     2    5
#14 GO:0030030 0.0046039296 10.973262 0.35460993     3   20
#15 GO:0048858 0.0046039296 10.973262 0.35460993     3   20
#16 GO:0032990 0.0046039296 10.973262 0.35460993     3   20
#17 GO:0048667 0.0060048728 24.000000 0.12411348     2    7
#18 GO:0000904 0.0060048728 24.000000 0.12411348     2    7
#19 GO:0000902 0.0060749734  9.803828 0.39007092     3   22
#20 GO:0071842 0.0101036289  4.493151 1.38297872     5   78
#21 GO:0071841 0.0112297747  4.366667 1.41843972     5   80
#22 GO:0009605 0.0124538818 14.967391 0.17730496     2   10
#23 GO:0007040 0.0177304965       Inf 0.01773050     1    1
#24 GO:0007033 0.0177304965       Inf 0.01773050     1    1
#25 GO:0051168 0.0177304965       Inf 0.01773050     1    1
#26 GO:0008202 0.0177304965       Inf 0.01773050     1    1
#27 GO:0022412 0.0177304965       Inf 0.01773050     1    1
#28 GO:0006694 0.0177304965       Inf 0.01773050     1    1
#29 GO:0007306 0.0177304965       Inf 0.01773050     1    1
#30 GO:0006611 0.0177304965       Inf 0.01773050     1    1
#31 GO:0010628 0.0351589830 57.666667 0.03546099     1    2
#32 GO:0010604 0.0351589830 57.666667 0.03546099     1    2
#33 GO:0045935 0.0351589830 57.666667 0.03546099     1    2
#34 GO:0051173 0.0351589830 57.666667 0.03546099     1    2
#35 GO:0051254 0.0351589830 57.666667 0.03546099     1    2
#36 GO:0030703 0.0351589830 57.666667 0.03546099     1    2
#37 GO:0031328 0.0351589830 57.666667 0.03546099     1    2
#38 GO:0007304 0.0351589830 57.666667 0.03546099     1    2
#39 GO:0010557 0.0351589830 57.666667 0.03546099     1    2
#40 GO:0009891 0.0351589830 57.666667 0.03546099     1    2
#41 GO:0045893 0.0351589830 57.666667 0.03546099     1    2
#42 GO:0009653 0.0382238962  3.578231 1.31205674     4   74
#                                                                      Term
#1                                            neuron projection development
#2                                          neuron projection morphogenesis
#3                                                      motor axon guidance
#4                                                   neuron differentiation
#5                                                       neuron development
#6                                         cellular component morphogenesis
#7                                                    generation of neurons
#8                                                         cell development
#9                                                               chemotaxis
#10                                                                   taxis
#11                                                              locomotion
#12                                                           axon guidance
#13                                                            axonogenesis
#14                                            cell projection organization
#15                                           cell projection morphogenesis
#16                                                 cell part morphogenesis
#17                   cell morphogenesis involved in neuron differentiation
#18                          cell morphogenesis involved in differentiation
#19                                                      cell morphogenesis
#20                       cellular component organization at cellular level
#21         cellular component organization or biogenesis at cellular level
#22                                           response to external stimulus
#23                                                   lysosome organization
#24                                                    vacuole organization
#25                                                          nuclear export
#26                                               steroid metabolic process
#27     cellular process involved in reproduction in multicellular organism
#28                                            steroid biosynthetic process
#29                                               eggshell chorion assembly
#30                                             protein export from nucleus
#31                                  positive regulation of gene expression
#32                  positive regulation of macromolecule metabolic process
#33 positive regulation of nucleobase-containing compound metabolic process
#34              positive regulation of nitrogen compound metabolic process
#35                            positive regulation of RNA metabolic process
#36                                                      eggshell formation
#37                    positive regulation of cellular biosynthetic process
#38                                   chorion-containing eggshell formation
#39               positive regulation of macromolecule biosynthetic process
#40                             positive regulation of biosynthetic process
#41                     positive regulation of transcription, DNA-dependent
#42                                      anatomical structure morphogenesis
#	
	
	
	
	
	
	
	
	
	
	
# Any NSYN SNPs? ------------------------	
nsyn = read.table(file="/media/data1/forty3/drone/vcf_drone/exons1.eff",header=F,colClasses = rep("character", 17))
nsyn$SNP = paste(nsyn$V1, nsyn$V2, sep=":")
nsyn = nsyn[c(18,3,4,8,11,13,14)]	
test = merge(sigFST, nsyn, by = "SNP")	
test[test$V13=="NON_SYNONYMOUS_CODING",]
	#494 coding SNPs are nonsynonymous
test$V8[test$V13=="NON_SYNONYMOUS_CODING"]
	
	
	
# Low- high and FST between AMC ------------------------
	#Are these genes closer related to A, M, or C?
		#I think I need to do this with STRUCTURE.
	#CommonSNPs.sh 
	#Ran ADMIXTURE and selected bees at High FST SNPs were less "M" and "C" than control bees at the same SNPs
	#SelvsConADMIXTURE.txt	
admix = read.table(file="SelvsConADMIXTURE.txt",header=T)
	#significantly shifted to be less M and C so either 1) A is more hygienic or 2) in this mixed background, A at these sites is more hygienic.
		#http://www.tandfonline.com/doi/pdf/10.1080/0005772X.1998.11099408	
	#plotted with HighFSTAdmixPlot.r
		
		
		
		
		
# QTL Regions ------------------------
#From Oxley et al.
qtl.hyg1 = (sigFST[grep("^2.19",sigFST$SNP),])
qtl.hyg1$chrPOS = gsub(".*:","",qtl.hyg1$SNP)
qtl.hyg1$chrPOS = as.numeric(gsub(".*:","",qtl.hyg1$SNP))
qtl.hyg1[which(qtl.hyg1$chrPOS > 86589 & qtl.hyg1$chrPOS < 145006),]
	#4 snps here
		
qtl.hyg2 = (sigFST[grep("^5.14",sigFST$SNP),])
qtl.hyg2$chrPOS = gsub(".*:","",qtl.hyg2$SNP)
qtl.hyg2$chrPOS = as.numeric(gsub(".*:","",qtl.hyg2$SNP))
qtl.hyg2[which(qtl.hyg2$chrPOS > 527544 & qtl.hyg2$chrPOS < 1558442),]
     #238 snps here
		
qtl.hyg3.1 = (sigFST[grep("^16.2",sigFST$SNP),])
qtl.hyg3.1$chrPOS = gsub(".*:","",qtl.hyg3.1$SNP)
qtl.hyg3.1$chrPOS = as.numeric(gsub(".*:","",qtl.hyg3.1$SNP))
qtl.hyg3.1[which(qtl.hyg3.1$chrPOS > 42885),]
	#250 snps here
	
qtl.hyg3.2 = (sigFST[grep("^16.3",sigFST$SNP),])
qtl.hyg3.2$chrPOS = gsub(".*:","",qtl.hyg3.2$SNP)
qtl.hyg3.2$chrPOS = as.numeric(gsub(".*:","",qtl.hyg3.2$SNP))
	#none here

qtl.hyg3.3 = (sigFST[grep("^16.4",sigFST$SNP),])
qtl.hyg3.3$chrPOS = gsub(".*:","",qtl.hyg3.3$SNP)
qtl.hyg3.3$chrPOS = as.numeric(gsub(".*:","",qtl.hyg3.3$SNP))
qtl.hyg3.3[which(qtl.hyg3.3$chrPOS < 920721),]
	#202 here

		
		

#And in genes?
qtlgene.hyg1 = (Genic.Fst1[grep("^2.19",Genic.Fst1$SNP),])
qtlgene.hyg1$chrPOS = gsub(".*:","",qtlgene.hyg1$SNP)
qtlgene.hyg1$chrPOS = as.numeric(gsub(".*:","",qtlgene.hyg1$SNP))
qtlgene.hyg1[which(qtlgene.hyg1$chrPOS > 86589 & qtlgene.hyg1$chrPOS < 145006),]
	#none
		
qtlgene.hyg2 = (Genic.Fst1[grep("^5.14",Genic.Fst1$SNP),])
qtlgene.hyg2$chrPOS = gsub(".*:","",qtlgene.hyg2$SNP)
qtlgene.hyg2$chrPOS = as.numeric(gsub(".*:","",qtlgene.hyg2$SNP))
qtlgene.hyg2[which(qtlgene.hyg2$chrPOS > 527544 & qtlgene.hyg2$chrPOS < 1558442),]
     #47 snps here, 16 gnes
		
qtlgene.hyg3.1 = (Genic.Fst1[grep("^16.2",Genic.Fst1$SNP),])
qtlgene.hyg3.1$chrPOS = gsub(".*:","",qtlgene.hyg3.1$SNP)
qtlgene.hyg3.1$chrPOS = as.numeric(gsub(".*:","",qtlgene.hyg3.1$SNP))
qtlgene.hyg3.1[which(qtlgene.hyg3.1$chrPOS > 42885),]
	#none
	
qtlgene.hyg3.3 = (Genic.Fst1[grep("^16.4",Genic.Fst1$SNP),])
qtlgene.hyg3.3$chrPOS = gsub(".*:","",qtlgene.hyg3.3$SNP)
qtlgene.hyg3.3$chrPOS = as.numeric(gsub(".*:","",qtlgene.hyg3.3$SNP))
qtlgene.hyg3.3[which(qtlgene.hyg3.3$chrPOS < 920721),]
     #81 snps here, 15 gnes





#UBX is in there, hmm 
			#http://dev.biologists.org/content/135/20/3435.full
			#http://www.sdbonline.org/sites/fly/segment/ultrabt5.htm
			# develops neurons, prevents death.
#and Maxillopadia
			#http://www.sdbonline.org/sites/fly/segment/probosc2.htm
			#http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0014536
			#MAxila have Octomamene receptors in the brain (re: SPivak paper)
			#http://www.researchgate.net/profile/Hans_Smid/publication/40113610_Octopamine-like_immunoreactivity_in_the_brain_and_suboesophageal_ganglion_of_two_parasitic_wasps_Cotesia_glomerata_and_Cotesia_rubecula/links/02e7e52fc87ac6bc48000000.pdf#page=61
			#Seratonin suppresses ant feeding, OA? http://www.sciencedirect.com/science/article/pii/S0022191011002605
		
	




	
		
###########################################################################################################
vcftools --vcf HYGmergedAMC.vcf  --weir-fst-pop  /media/data1/afz/git/M.plink  --weir-fst-pop /media/data1/forty3/drone/vcf_drone/pop2.txt --maf 0.05 --out MvLHB
vcftools --vcf HYGmergedAMC.vcf  --weir-fst-pop  /media/data1/afz/git/M.plink  --weir-fst-pop /media/data1/forty3/drone/vcf_drone/Selpop.txt --maf 0.05 --out MvHHB

vcftools --vcf HYGmergedAMC.vcf  --weir-fst-pop  /media/data1/afz/git/C.plink  --weir-fst-pop /media/data1/forty3/drone/vcf_drone/pop2.txt --maf 0.05 --out CvLHB
vcftools --vcf HYGmergedAMC.vcf  --weir-fst-pop  /media/data1/afz/git/C.plink  --weir-fst-pop /media/data1/forty3/drone/vcf_drone/Selpop.txt --maf 0.05 --out CvHHB


#Isolate all common SNPs
fstLM = read.table(file="MvLHB.weir.fst",header=T,colClasses=rep("character",3));names(fstLM)[3]="FSTLM"
fstLM = fstLM[which(as.numeric(fstLM$FST)>0),]
fstLM$SNP = paste(fstLM$CHROM,fstLM$POS, sep=":")

fstHM = read.table(file="MvHHB.weir.fst",header=T,colClasses=rep("character",3));names(fstHM)[3]="FSTHM"
fstHM = fstHM[which(as.numeric(fstHM$FST)>0),]
fstHM$SNP = paste(fstHM$CHROM,fstHM$POS, sep=":")

fstLC = read.table(file="CvLHB.weir.fst",header=T,colClasses=rep("character",3));names(fstLC)[3]="FSTLC"
fstLC = fstLC[which(as.numeric(fstLC$FST)>0),]
fstLC$SNP = paste(fstLC$CHROM,fstLC$POS, sep=":")

fstHC = read.table(file="CvHHB.weir.fst",header=T,colClasses=rep("character",3));names(fstHC)[3]="FSTHC"
fstHC = fstHC[which(as.numeric(fstHC$FST)>0),]
fstHC$SNP = paste(fstHC$CHROM,fstHC$POS, sep=":")


snps=(Reduce(intersect, list(fstLM$SNP ,fstHM$SNP,fstLC$SNP ,fstHC$SNP )))
fstLM=fstLM[fstLM$SNP %in% snps,]
fstHM=fstHM[fstHM$SNP %in% snps,]
fstLC=fstLC[fstLC$SNP %in% snps,]
fstHC=fstHC[fstHC$SNP %in% snps,]

conv=data.frame(cbind(fstLM[4],fstHM[3],fstHC[3],fstLM[3],fstLC[3] ))



#as per Zayed and Whitfield
conv$Hdis=with(conv,(as.numeric(FSTHM)/(as.numeric(FSTHM) + as.numeric(FSTHC)))) #Higher Hdis indicates more M
conv$Ldis=with(conv,(as.numeric(FSTLM)/(as.numeric(FSTLM) + as.numeric(FSTLC))))
conv$difs =  conv$Hdis - conv$Ldis # high difs indicates more M in selected


#difs in sigFST
boxplot(conv$difs, conv$difs[conv$SNP %in% sigFST$SNP],	notch=T	)
	#quick perm test
	#pnorm(mean(conv$difs[conv$SNP %in% sigFST$SNP]), mean=mean(perm.conv), sd = sd(perm.conv))



boxplot(as.numeric(conv$FSTHM[conv$SNP %in% sigFST$SNP]),as.numeric(conv$FSTLM[conv$SNP %in% sigFST$SNP]))		
		
		
		
		

x11();boxplot(as.numeric(fstHC$FSTHC[fstHC$SNP %in% sigFST$SNP]),as.numeric(fstHM$FSTHM[fstHM$SNP %in% sigFST$SNP]))
	# High bee SNPs become more "C" like.
	
	
boxplot(as.numeric(fstHM$FSTHM[fstHM$SNP %in% sigFST$SNP]),
		as.numeric(fstLM$FSTLM[fstLM$SNP %in% sigFST$SNP]), 
		as.numeric(fstHM$FSTHM),
		as.numeric(fstLM$FSTLM),
		names = c("Selected High FST",
			"Control High FST",
			"Selected All",
			"Control All"),
		notch=T	)

	
boxplot(as.numeric(fstHC$FSTHC[fstHC$SNP %in% sigFST$SNP]),
	as.numeric(fstLC$FSTLC[fstLC$SNP %in% sigFST$SNP]), 
	as.numeric(fstHC$FSTHC),
	as.numeric(fstLC$FSTLC),
	names = c("Selected High FST",
		"Control High FST",
		"Selected All",
		"Control All"),
	notch=T	)

#Selected bees less at those SNPs.
	
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
boxplot(as.numeric(fstHM$FSTHM[fstHM$SNP %in% test$SNP]),
		as.numeric(fstLM$FSTLM[fstLM$SNP %in% test$SNP]), 
		as.numeric(fstHM$FSTHM),
		as.numeric(fstLM$FSTLM),
		names = c("Selected High FST",
			"Control High FST",
			"Selected All",
			"Control All"),
		notch=T	)

	
x11();boxplot(as.numeric(fstHC$FSTHC[fstHC$SNP %in% test$SNP]),
	as.numeric(fstLC$FSTLC[fstLC$SNP %in% test$SNP]), 
	as.numeric(fstHC$FSTHC),
	as.numeric(fstLC$FSTLC),
	names = c("Selected High FST",
		"Control High FST",
		"Selected All",
		"Control All"),
	notch=T	)
	
	
	
	
	
	