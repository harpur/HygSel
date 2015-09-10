###
#
###


#Analysis of Fst following pFST (pFST.sh) and VCFFst (DroneFst.r)
	#RE: DRIFT http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3878089/

#save.image(file="/media/data1/forty3/drone/FST/SelvsCon/SelectedvsControlFSTPvals.RData")	
	


# Load datasets and functions ---------------------------
source("/media/data1/forty3/brock/scripts/VarFunct.r")
source("/media/data1/forty3/brock/scripts/GOgetter.r")
load(file="/media/data1/forty3/drone/FST/SelvsCon/SelectedvsControlFSTPvals.RData")
gff = read.table(file = "/media/data1/forty3/brock/scripts/GFF/NCBIGFF.txt")
gff=gff[gff$V5=="gene",]
degs = read.table(file="/media/data1/forty3/drone/git/boutinDEGs.txt",header=T)

# Analysis ---------------------------
test = read.table("pFstDRONE.counts",header=F,colClasses=c("character","character","character"))
test$SNP = paste(test$V1,test$V2,sep=":")
test = merge(fst23,test, by="SNP")
test = test[-c(5,6)]
names(test)[5]="P"
test$q=qvalue0(as.numeric(test$P))$q
PFST = test
sigFST = PFST[PFST$q < 0.01,]
#sigFST1 = PFST[PFST$q < 0.001,]
#53894 SNPs significant at q<0.01
#1101 at Q<0.001 (sigFST1)

# Genes with high FST SNPs -------------------	
Genic.Fst = c()
chrom = intersect(gff$V1,sigFST$CHROM)
for(i in chrom){
	win.temp = gff[gff$V1==i,]
	deg.temp = sigFST[which(as.character(sigFST$CHROM)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$V3)), ">=") 
	blah1=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$V4)), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	Genic.Fst = rbind(temp,Genic.Fst)
	print(i)
}
names(Genic.Fst)[8] ="GB"
high.genes = as.character(unique(Genic.Fst$GB))
	#1973 genes total.
	


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
obs = aggregate(Genic.Fst$q, by=list(Genic.Fst$GB), length ) 	
names(obs)[2] = "obs"	
test = merge(obs, expec, by = "Group.1")	

test$p = apply(test, 1, function(x) pnorm(as.numeric(x[2]),mean=as.numeric(x[3]),sd=as.numeric(x[4]),lower.tail=FALSE))	
names(test)[1]="GB"	
sigFSTperm = test #contains significant SNPs within gene regions.
high.genes.perm = as.character(sigFSTperm$GB[sigFSTperm$p<0.05])
# Write list for GO analysis -------------------	
FSTFBGN = GO.getter(test[test$p<0.05,])
FSTFBGN = FSTFBGN[!duplicated(FSTFBGN$FGN),]
FSTFBGN = FSTFBGN[!is.na(FSTFBGN$FGN),]
write.list(FSTFBGN$FGN, file="GOFst.perm")
	
	
	
#After permutation, repeat GOSTATS 
 
 
 
 
 
	
# NSYN SNPs in my highpFST list? ------------------------
	#From SNPEFF.sh
	
nsyn = read.table(file="/media/data1/forty3/drone/vcf_drone/exons1.eff",header=F,colClasses = rep("character", 17))
nsyn$SNP = paste(nsyn$V1, nsyn$V2, sep=":")
nsyn = nsyn[c(18,3,4,8,11,13,14)]	
test = merge(sigFST, nsyn, by = "SNP")	
test[test$V13=="NON_SYNONYMOUS_CODING",]
	#494 coding SNPs are nonsynonymous
high.genes.perm.nsyn = as.character(unique(test$V8[test$V13=="NON_SYNONYMOUS_CODING"]))
high.genes.perm.nsyn = high.genes.perm.nsyn[high.genes.perm.nsyn   %in% high.genes.perm]
	#148 genes
	
	
	
	
# GO Analyses ------------------------	
#GO with genes with permuted genes:
test = data.frame(high.genes.perm);names(test)="GB"
FSTFBGN = GO.getter(test)
FSTFBGN = FSTFBGN[!duplicated(FSTFBGN$FGN),]
FSTFBGN = FSTFBGN[!is.na(FSTFBGN$FGN),]
write.list(FSTFBGN$FGN, file="GOFst.high.genes.perm")

#       GOBPID      Pvalue OddsRatio  ExpCount Count Size
#1  GO:0009605 0.002325024 10.313725 0.6312057     4   10
#5  GO:0031175 0.010802103  9.181395 0.5049645     3    8
#6  GO:0048812 0.010802103  9.181395 0.5049645     3    8
#7  GO:0048468 0.010903963  3.747470 1.9567376     6   31
#8  GO:0034660 0.011339527 30.344828 0.1893617     2    3
#9  GO:0034470 0.011339527 30.344828 0.1893617     2    3
#10 GO:0030707 0.011339527 30.344828 0.1893617     2    3
#11 GO:0016319 0.014620276  4.078947 1.5148936     5   24
#12 GO:0007420 0.017383564  3.872024 1.5780142     5   25
#13 GO:0030182 0.027713892  5.725291 0.6943262     3   11
#14 GO:0048666 0.027713892  5.725291 0.6943262     3   11
#15 GO:0006935 0.034770292 10.099617 0.3156028     2    5
#16 GO:0042330 0.034770292 10.099617 0.3156028     2    5
#17 GO:0007411 0.034770292 10.099617 0.3156028     2    5
#18 GO:0007409 0.034770292 10.099617 0.3156028     2    5
#19 GO:0040011 0.034770292 10.099617 0.3156028     2    5
#20 GO:0048699 0.035295719  5.085271 0.7574468     3   12
#21 GO:0009653 0.038679583  2.173846 4.6709220     9   74
#22 GO:0009628 0.043834224  4.573256 0.8205674     3   13
#23 GO:0007472 0.045212051  3.406536 1.3886525     4   22
#24 GO:0007476 0.045212051  3.406536 1.3886525     4   22
#25 GO:0048736 0.045212051  3.406536 1.3886525     4   22
#26 GO:0048737 0.045212051  3.406536 1.3886525     4   22
#27 GO:0035120 0.045212051  3.406536 1.3886525     4   22
#28 GO:0035114 0.045212051  3.406536 1.3886525     4   22
#29 GO:0035107 0.045212051  3.406536 1.3886525     4   22
#                                            Term
#1                  response to external stimulus
#5                  neuron projection development
#6                neuron projection morphogenesis
#7                               cell development
#8                        ncRNA metabolic process
#9                               ncRNA processing
#10             ovarian follicle cell development
#11                     mushroom body development
#12                             brain development
#13                        neuron differentiation
#14                            neuron development
#15                                    chemotaxis
#16                                         taxis
#17                                 axon guidance
#18                                  axonogenesis
#19                                    locomotion
#20                         generation of neurons
#21            anatomical structure morphogenesis
#22                  response to abiotic stimulus
#23                       wing disc morphogenesis
#24      imaginal disc-derived wing morphogenesis
#25                         appendage development
#26   imaginal disc-derived appendage development
#27        post-embryonic appendage morphogenesis
#28 imaginal disc-derived appendage morphogenesis
#29                       appendage morphogenesis





test = data.frame(high.genes.perm.nsyn);names(test)="GB"
FSTFBGN = GO.getter(test)
FSTFBGN = FSTFBGN[!duplicated(FSTFBGN$FGN),]
FSTFBGN = FSTFBGN[!is.na(FSTFBGN$FGN),]
write.list(FSTFBGN$FGN, file="GOFst.high.genes.perm.nsyn")

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
#11 GO:0007411 0.0029225832 40.057971 0.08865248     2    5
#12 GO:0007409 0.0029225832 40.057971 0.08865248     2    5
#13 GO:0040011 0.0029225832 40.057971 0.08865248     2    5
#14 GO:0030030 0.0046039296 10.973262 0.35460993     3   20
#15 GO:0032990 0.0046039296 10.973262 0.35460993     3   20
#16 GO:0048858 0.0046039296 10.973262 0.35460993     3   20
#17 GO:0000904 0.0060048728 24.000000 0.12411348     2    7
#18 GO:0048667 0.0060048728 24.000000 0.12411348     2    7
#19 GO:0000902 0.0060749734  9.803828 0.39007092     3   22
#20 GO:0071842 0.0101036289  4.493151 1.38297872     5   78
#21 GO:0071841 0.0112297747  4.366667 1.41843972     5   80
#22 GO:0009605 0.0124538818 14.967391 0.17730496     2   10
#23 GO:0051168 0.0177304965       Inf 0.01773050     1    1
#24 GO:0007306 0.0177304965       Inf 0.01773050     1    1
#25 GO:0008202 0.0177304965       Inf 0.01773050     1    1
#26 GO:0022412 0.0177304965       Inf 0.01773050     1    1
#27 GO:0006694 0.0177304965       Inf 0.01773050     1    1
#28 GO:0006611 0.0177304965       Inf 0.01773050     1    1
#29 GO:0007040 0.0177304965       Inf 0.01773050     1    1
#30 GO:0007033 0.0177304965       Inf 0.01773050     1    1
#31 GO:0051173 0.0351589830 57.666667 0.03546099     1    2
#32 GO:0007304 0.0351589830 57.666667 0.03546099     1    2
#33 GO:0051254 0.0351589830 57.666667 0.03546099     1    2
#34 GO:0030703 0.0351589830 57.666667 0.03546099     1    2
#35 GO:0045935 0.0351589830 57.666667 0.03546099     1    2
#36 GO:0010557 0.0351589830 57.666667 0.03546099     1    2
#37 GO:0045893 0.0351589830 57.666667 0.03546099     1    2
#38 GO:0010628 0.0351589830 57.666667 0.03546099     1    2
#39 GO:0010604 0.0351589830 57.666667 0.03546099     1    2
#40 GO:0009891 0.0351589830 57.666667 0.03546099     1    2
#41 GO:0031328 0.0351589830 57.666667 0.03546099     1    2
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
#11                                                           axon guidance
#12                                                            axonogenesis
#13                                                              locomotion
#14                                            cell projection organization
#15                                                 cell part morphogenesis
#16                                           cell projection morphogenesis
#17                          cell morphogenesis involved in differentiation
#18                   cell morphogenesis involved in neuron differentiation
#19                                                      cell morphogenesis
#20                       cellular component organization at cellular level
#21         cellular component organization or biogenesis at cellular level
#22                                           response to external stimulus
#23                                                          nuclear export
#24                                               eggshell chorion assembly
#25                                               steroid metabolic process
#26     cellular process involved in reproduction in multicellular organism
#27                                            steroid biosynthetic process
#28                                             protein export from nucleus
#29                                                   lysosome organization
#30                                                    vacuole organization
#31              positive regulation of nitrogen compound metabolic process
#32                                   chorion-containing eggshell formation
#33                            positive regulation of RNA metabolic process
#34                                                      eggshell formation
#35 positive regulation of nucleobase-containing compound metabolic process
#36               positive regulation of macromolecule biosynthetic process
#37                     positive regulation of transcription, DNA-dependent
#38                                  positive regulation of gene expression
#39                  positive regulation of macromolecule metabolic process
#40                             positive regulation of biosynthetic process
#41                    positive regulation of cellular biosynthetic process
#42                                      anatomical structure morphogenesis
#

















	
		
		
# QTL Regions ------------------------
	#This is manual, for now.

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
	
	
	
	
#Spotter's QTLs?		

#And in genes?
qtlgene.hyg1.spot = (Genic.Fst1[which(Genic.Fst1$CHROM=="1"),])
qtlgene.hyg1.spot = qtlgene.hyg1.spot[which(qtlgene.hyg1.spot$POS > 3039231 & qtlgene.hyg1.spot$POS < 8453574),]
	#66 SNPs, 7 genes 
qtlgene.hyg1.1.spot = (Genic.Fst1[which(Genic.Fst1$CHROM=="1"),])
qtlgene.hyg1.1.spot = qtlgene.hyg1.1.spot[which(qtlgene.hyg1.1.spot$POS > 9418717 & qtlgene.hyg1.1.spot$POS < 16819942),]
lenunique(qtlgene.hyg1.1.spot$GB)
nrow(qtlgene.hyg1.1.spot)
	#299 SNPs, 37 genes 
qtlgene.hyg2.spot = (Genic.Fst1[which(Genic.Fst1$CHROM=="2"),])
qtlgene.hyg2.spot = qtlgene.hyg2.spot[which(qtlgene.hyg2.spot$POS > 1 & qtlgene.hyg2.spot$POS < 12503099),]
lenunique(qtlgene.hyg2.spot$GB)
nrow(qtlgene.hyg2.spot)
	#852 SNPs, 94 genes 
		
qtlgene.hyg6.spot = (Genic.Fst1[which(Genic.Fst1$CHROM=="2"),])
qtlgene.hyg6.spot = qtlgene.hyg6.spot[which(qtlgene.hyg6.spot$POS > 11206828 & qtlgene.hyg6.spot$POS < 17739083),]
lenunique(qtlgene.hyg6.spot$GB)
nrow(qtlgene.hyg6.spot)
	#94 SNPs, 7 genes 	
		

qtlgene.hyg7.spot = (Genic.Fst1[which(Genic.Fst1$CHROM=="2"),])
qtlgene.hyg7.spot = qtlgene.hyg7.spot[which(qtlgene.hyg7.spot$POS > 9515998 & qtlgene.hyg7.spot$POS < 12848973),]
lenunique(qtlgene.hyg7.spot$GB)
nrow(qtlgene.hyg7.spot)
	#187 SNPs, 21 genes 			
		
qtlgene.hyg12.spot = (Genic.Fst1[which(Genic.Fst1$CHROM=="2"),])
qtlgene.hyg12.spot = qtlgene.hyg12.spot[which(qtlgene.hyg12.spot$POS > 1 & qtlgene.hyg12.spot$POS < 4003353),]
lenunique(qtlgene.hyg12.spot$GB)
nrow(qtlgene.hyg12.spot)
	#72 SNPs, 9 genes 
	
qtlgene.hyg13.spot = (Genic.Fst1[which(Genic.Fst1$CHROM=="2"),])
qtlgene.hyg13.spot = qtlgene.hyg13.spot[which(qtlgene.hyg13.spot$POS > 5247545 & qtlgene.hyg13.spot$POS < 10266737),]
lenunique(qtlgene.hyg13.spot$GB)
nrow(qtlgene.hyg13.spot)		
		#373 SNPs, 50 genes 

qtlgene.hyg15.spot = (Genic.Fst1[which(Genic.Fst1$CHROM=="2"),])
qtlgene.hyg15.spot = qtlgene.hyg15.spot[which(qtlgene.hyg15.spot$POS > 1 & qtlgene.hyg15.spot$POS < 6643609),]
lenunique(qtlgene.hyg15.spot$GB)
nrow(qtlgene.hyg15.spot)	
		#583 SNPs, 53 genes 		

		
qtlgene.hyg16.spot = (Genic.Fst1[which(Genic.Fst1$CHROM=="2"),])
qtlgene.hyg16.spot = qtlgene.hyg16.spot[which(qtlgene.hyg16.spot$POS > 3196393 & qtlgene.hyg16.spot$POS < 6242592),]
lenunique(qtlgene.hyg16.spot$GB)
nrow(qtlgene.hyg16.spot)			
		#366 SNPs, 39 genes 		



		
		
# testing REHH aginast FST results ---------------------------
	#http://www.broadinstitute.org/scientific-community/science/programs/medical-and-population-genetics/cms/overview-0
	#http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0050171
	#http://www.cs.cmu.edu/~sssykim/teaching/f11/slides/Lecture10.pdf
	#http://haplotter.uchicago.edu/instruction.html
	#
	#OK, so there is definitely overlap
	#Run REHH.sh -> IHS.out list of significant SNPs
rehh = read.table(file="/media/data1/forty3/drone/vcf_drone/IHS.out",header=F)
names(rehh)=c("CHROM", "POS","iHS","P","SNP")	

out=c()
for(i in c(1:16)){
	test1 = sigFST[sigFST$CHROM==i,]
	chr1 = rehh[rehh$CHROM==i,]
	blah=outer(chr1$POS, test1$POS-500, ">=") 
	blah1=outer(chr1$POS, test1$POS+500, "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	out = rbind(out, (test1[blah[,2],]))
	print(i)
}

out = out[!duplicated(out),]
	#6106 SNPs


# Genes with high FST SNPs -------------------	
Genic.Fst.REHH = c()
for(i in chrom){
	win.temp = gff[gff$V1==i,]
	deg.temp = out[which(as.character(out$CHROM)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$V3)), ">=") 
	blah1=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$V4)), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	Genic.Fst.REHH = rbind(temp,Genic.Fst.REHH)
	print(i)
}








	
		
		
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
		
		
		
		
# Longer-term selection ------------------------			
	#Gamma in hygiene-associated genes? TD? Pi?

#load datasets (Pi, TD, gamma)
imm.new =  read.table(file = "EvansNewImm.txt", header=T)
gamm = read.table(file = "/media/data1/forty3/brock/immune/Gene_gamma.txt", header=T)
kaks = read.table(file="/media/data1/forty3/brock/immune/KAKS.out")
names(kaks) =  c("GB", "c.ks", "m.ks", "y.ks", "a.ks", "c.ka", "m.ka", "y.ka", "a.ka" )
kaks$kaks = kaks$a.ka/kaks$a.ks
kaks = kaks[kaks$a.ks>0,]
Pi = read.table(file="/media/data1/forty3/brock/immune/A.snpeff1.eff.Pi")
names(Pi) = c("GB", "NS","NN", "NumS", "NumN", "PiS", "PiN")
source("/media/data1/forty3/brock/scripts/VarFunct.r")
TD = read.table(file="/media/data1/forty3/brock/balancedSNPs/data/A.recode.vcf.ETD")
TD = aggregate(TD$V4, by=list(TD$V1),mean); names(TD)=c("GB","TD")



	
	
#Gamma
test = gamm ; qwd = rep("N", nrow(test)); 
qwd[test$GB %in%  high.genes.perm.nsyn] = "HYG" 
test$qwd = qwd

boxplot(test$gamma~test$qwd)
	#HYG is higher.....
perm.test(test$gamma, 360, test$gamma[(test$qwd=="HYG")], 10000)		
		# P (hyg gene Gamma) is 0.001652618, so it's significant. (with hyg genes)
#So these genes are under selection voer long term.	And likely more genes under purifying selection in rest of genome
		
#divergence and Pi and TD

test2 = merge(Pi, kaks, by="GB")
test2 = test2[which(test2$NumS>5),]
test = merge(test, test2, by = "GB")
test = merge(test, TD, by = "GB")
boxplot(test$PiS~test$qwd,notch=T) #NS
boxplot(test$PiN~test$qwd,notch=T) #higher (and likely higher than IMM too)
boxplot(test$kaks~test$qwd,notch=T) #higher (and likely higher than IMM too)
boxplot(test$TD~test$qwd,notch=T) #slightly higher, NS
boxplot(test$NS~test$qwd,notch=T) #NS
	#positive selection. 

	
	
	

############
#VS Immune genes:	
imm.new =  read.table(file = "/media/data1/forty3/brock/immune/EvansNewImm.txt", header=T) 
qwd = test$qwd
qwd[test$GB %in%  imm.new$GB] = "IMM" 	
test$qwd = qwd	


	
#> summary(aov(test$gamma~test$qwd))
#               Df Sum Sq Mean Sq F value   Pr(>F)
#test$qwd        2     11   5.665   10.64 2.41e-05 ***
#Residuals   12300   6546   0.532
#---
#Signif. codes:  0 â***â 0.001 â**â 0.01 â*â 0.05 â.â 0.1 â â 1
#> TukeyHSD(aov(test$gamma~test$qwd))
#  Tukey multiple comparisons of means
#    95% family-wise confidence level
#
#Fit: aov(formula = test$gamma ~ test$qwd)
#
#$`test$qwd`
#               diff         lwr         upr     p adj
#IMM-HYG  0.09289207 -0.06300664  0.24879079 0.3426527
#N-HYG   -0.11008321 -0.20220150 -0.01796492 0.0141234
#N-IMM   -0.20297528 -0.33070817 -0.07524239 0.0005755
#	
#################	
	
	
	
	
	
	
# QWD biased genes? ------------------------			
	#high.genes.perm in QWD List

load(file="/media/data1/forty3/brock/expression_data/DWQ_expression.RData")
#for exclusive genes:
#follow CK's instructions:
Q<-as.character(Grozinger2007QueenGenes1$GB)
W<-as.character(Grozinger2007WorkerGenes1$GB)
g<-as.character(ZayedDroneWorker$GB); D<-g[ZayedDroneWorker$W.D.1=="D" & (ZayedDroneWorker$casteF1_fdr<.05)]
D1<-setdiff(D,W); D1<-setdiff(D1,Q); W1<-setdiff(W,D); Q1<-setdiff(Q,D)	

test = gff ; qwd = rep("N", nrow(test)); 
qwd[test$V2 %in% Q1] = "Q" 
qwd[test$V2 %in% W1] = "W" 
qwd[test$V2 %in% D1] = "D"   	
test$qwd = qwd

qwd = rep("N", nrow(test)); 
qwd[test$V2 %in% high.genes.perm] = "sig" 
test$sig = qwd
aggregate(test$sig, by = list(test$qwd), table)
#  Group.1   x.N x.sig
#1       D   815    46
#2       N 11140   567
#3       Q   366    14
#4       W   317    20


#no tendancy to be caste-bias in expression.

	

	

	
	
# Do they make up any known networks? Are they TFs? Are they central? ------------------------			
	#high.genes.perm vs TRN 
test = trn 
qwd = rep("N", nrow(test)); 
qwd[test$GB %in% high.genes.perm.nsyn] = "sig" 
test$sig = qwd
aggregate(test$sig, by = list(test$TFTARG), table)
#  Group.1  x.N x.sig
#1    TARG 1555    40
#2      TF  182     7


	#NS (Fisher exact test)
	#no difference in connectedness|TF or targ (AOV)
#significant TFs:	
#GB48271	13c4g	FBgn0000210	Br-c	GB14070	-	broad-complex
#GB47052	3c15g	NA	LOC725508	GB18040	-	uncharacterized LOC725508
#GB51521	3c7g	FBgn0052105	LOC724778	GB19877	-	LIM homeobox transcription factor 1-beta-like
#GB49869	13c12g	FBgn0032904	Mtp	GB16245	-	microsomal triacylglycerol transfer protein
#GB53328	9c12g	FBgn0024887	kin17	GB16193	-	kin17 protein
#GB55837	3c8g	FBgn0032940	Mio	GB12214	-	Mlx interactor
#GB49969	15c19g	FBgn0039530	LOC409680	GB18478	-	tubby-related protein 4-like
#

	
	

	
	
	
# Old versus new genes? ------------------------		
trg = read.table(file="clipboard",header=T) #this is straigth from my PNAS paper (SD5)
test = trg 
qwd = rep("N", nrow(test)); 
qwd[test$GB %in% high.genes.perm] = "sig" 
test$sig = qwd
aggregate(test$sig, by = list(test$Taxa), table)	
#      Group.1  x.N x.sig
#1        Apis   86     2
#2     Apoidea  210     5
#3 Hymenoptera 1265    56
#4     Insecta 8241   445
#
	# Insecta versus non-insect (P value equals 0.0323) is significant.
	#More conserved genes than less conserved genes.
	
	
	
# Expression? ------------------------		
	#Use Johnson's Data (Johnson_TableS3.xlsx)
exp = read.table(file="clipboard",header=T)
test = exp 
qwd = rep("N", nrow(test)); 
qwd[test$GB %in% high.genes.perm] = "sig" 
test$sig = qwd	
aggregate(test$sig, by = list(test$Taxonomy), table)		
	#Similar to above, "Conserved" have more The two-tailed P value is less than 0.0001
#      Group.1  x.N x.sig
#1   Arthropod  208     6
#2         Bee   24     0
#3   Conserved 5150   264
#4 Hymenoptera  385     7
#5      Insect  549    14
#6      Orphan  191     0
#7      O-TRGs  780    27
#

#Not highly expressed (relative to all expression) in all tissues.

boxplot(log10(1+test$brn[test$brn>25])~test$sig[test$brn>25])
	#check insect and hymenoptera and O-TRG
test2 = test[which(test$brn>25),] 
boxplot(log10(1+test2$brn[test2$Taxonomy=="O-TRGs"])~test2$sig[test2$Taxonomy=="O-TRGs"])
#They are not over-expressed in any single Nurse tissue (Brain included).





# Overlap with Parker et al. data
	#http://www.genomebiology.com/2012/13/9/R81
	#I think the additional dataset labels might be incorrect-> ask LF
	# I'll read through the paper instead....
	#And I can't find some of these genes in AMEL....
#
#GID	local	GB
#Fas1	lobe	GB44844
#Fas1	lobe	GB44846
#Fas1	lobe	GB44847
#Fas1	lobe	GB44848
#Fas1	lobe	GB44849
#LamininA	lobe	NA17
#Ankerin2	lobe	GB50172
#Ankerin2	lobe	GB50368
#Ankerin2	lobe	GB50369
#Amph	nerve	GB43280
#








#SNP Characterization
	#What do the SNPs in genes I've got actually do?

	
#Run on AfrSNPs.raw.vcf, it contains all SNPs 
java -jar /usr/local/lib/snpEff2/snpEff.jar Amel -o txt --vcf /media/data1/forty3/drone/vcf_drone/DroneSelectionFinal.recode.vcf  -ss  > DRONESNP.snpeff.eff
sed '/SYNONYMOUS/ !d' Capensis.snpeff.eff > exons.eff
	#fixed ?/?
sed '/WARNING/ !d' exons.eff > warnexons.eff
sed '/WARNING/ d' exons.eff > exons1.eff






















# UTRs? ------------------------		
#Re-run SNPEFF with upstream/downstream stuff
	#Run on AfrSNPs.raw.vcf, it contains all SNPs 
	#java -jar /usr/local/lib/snpEff2/snpEff.jar Amel -o txt DroneSelectionFinal.recode.vcf    > DRONE.UP.snpeff.eff
	#sed '/UPSTREAM/ !d' DRONE.UP.snpeff.eff > UPSTREAM.eff
	#sed '/WARNING/ d' UPSTREAM.eff > UPSTREAM1.eff
upstream = read.table(file="/media/data1/forty3/drone/vcf_drone/UPSTREAM1.eff", header=F,colClasses=rep("character",13))
upstream = upstream[c(1,2,8,12)];upstream$SNP = paste(upstream$V1,upstream$V2,sep=":")
test = merge(upstream, sigFST, by = "SNP")
names(test)[4] = "GB"
#3268 genes with high FST upstream SNPs (5kb upstream)
#lenunique(test$GB[as.numeric(test$V12) <= 1000,])
	#1133<1000

		

# Which UTRs have more high FST than by chance alone? -------------------	
	#Generate expected distribution to see if number of high FST snps per gene is expected by chance alone.
sum.up = aggregate(upstream$V8, by = list(upstream$V8), length)
obs.up = aggregate(test$GB, by=list(test$GB),length)
sum.up = sum.up[sum.up$Group.1 %in% obs.up$Group.1,]
qs = PFST$q
exp = sum.up$x
obs = as.numeric(obs.up$x)
p=c()

for(i in 1:length(exp)){
	out = c()
	for(u in 1:100){
	x=sample(qs, exp[i])
	out = c(out,length(x[x<0.01]))
	}
	p = c(p,pnorm(obs[i], mean = mean(out), sd = sd(out),lower.tail=FALSE))
	print(i)
}
	
obs.up$exp = sum.up$x 
obs.up$p = p 
head(obs.up[obs.up$p<0.05,])	
UTR.Pvalues = obs.up
high.utr.perm = obs.up$Group.1[obs.up$p<0.05]	
	#625 significantly more than expected.
	
	
	
	
head(deps[deps %in% high.utr.perm ])
#"None

head(degs[degs$GB %in% high.utr.perm,])
#        Gene Chromosome      GB Reg
#11 LOC552229          1 GB46514 DWN
#34 LOC408414         13 GB47990 DWN
#37 LOC726040         13 GB40077 DWN
#74 LOC724749          7 GB42493  UP

test = data.frame(high.utr.perm);names(test)="GB"
FSTFBGN = GO.getter(test)
FSTFBGN = FSTFBGN[!duplicated(FSTFBGN$FGN),]
FSTFBGN = FSTFBGN[!is.na(FSTFBGN$FGN),]
write.list(FSTFBGN$FGN, file="GOFst.UTRs")


summary(outcomeBP) #Goode AOF
#       GOBPID      Pvalue OddsRatio    ExpCount Count Size
#1  GO:0003006 0.004886929 12.008929  0.37659574     3    9
#2  GO:0048610 0.004886929 12.008929  0.37659574     3    9
#3  GO:0030707 0.005027928 47.368421  0.12553191     2    3
#4  GO:0048468 0.007948047  4.718661  1.29716312     5   31
#5  GO:0030154 0.025461914  1.933860  9.79148936    16  234
#6  GO:0048869 0.041827350  1.794707 10.37730496    16  248
#7  GO:0051168 0.041843972       Inf  0.04184397     1    1
#8  GO:0007507 0.041843972       Inf  0.04184397     1    1
#9  GO:0007306 0.041843972       Inf  0.04184397     1    1
#10 GO:0072358 0.041843972       Inf  0.04184397     1    1
#11 GO:0072359 0.041843972       Inf  0.04184397     1    1
#12 GO:0022412 0.041843972       Inf  0.04184397     1    1
#13 GO:0007179 0.041843972       Inf  0.04184397     1    1
#14 GO:0006611 0.041843972       Inf  0.04184397     1    1
#15 GO:0045168 0.045688709  3.315361  1.38085106     4   33
#16 GO:0046331 0.045688709  3.315361  1.38085106     4   33
#                                                                  Term
#1                       developmental process involved in reproduction
#2                            cellular process involved in reproduction
#3                                    ovarian follicle cell development
#4                                                     cell development
#5                                                 cell differentiation
#6                                       cellular developmental process
#7                                                       nuclear export
#8                                                    heart development
#9                                            eggshell chorion assembly
#10                                   cardiovascular system development
#11                                      circulatory system development
#12 cellular process involved in reproduction in multicellular organism
#13          transforming growth factor beta receptor signaling pathway
#14                                         protein export from nucleus
#15                cell-cell signaling involved in cell fate commitment
#16                                                  lateral inhibition
#




#UTRs conserved?
exp = read.table(file="clipboard",header=T)
test = exp 
qwd = rep("N", nrow(test)); 
qwd[test$GB %in% high.utr.perm] = "sig" 
test$sig = qwd	
aggregate(test$sig, by = list(test$Taxonomy), table)
#      Group.1  x.N x.sig
#1   Arthropod  208     6
#2         Bee   24     0
#3   Conserved 5218   196
#4 Hymenoptera  385     7
#5      Insect  549    14
#6      Orphan  190     1
#7      O-TRGs  784    23
# The two-tailed P value equals 0.0050, again, mroe conserved genes. 



























































###########################################################################################










# Longer-term selection ------------------------
		#COMPARE SNP versus SNP Selected herein versus FST between.....

		
SC.pfst = PFST[c(1,4,5,6)]; names(SC.pfst)[4]="SCq"
M.pfst = read.table(file="M.counts",header=F)
M.pfst$q = qvalue0(M.pfst$V3)$q		
M.pfst$SNP = paste(M.pfst$V1, M.pfst$V2, sep=":")
test = merge(SC.pfst, M.pfst, by="SNP")
test$wt2 <- as.numeric(cut(test$SCq, 15))
test2 = aggregate(as.numeric(test$q), by =list(test$wt2), mean)
	#high as you go from 1-15
#M and SC SNP q values not correlated.
#C and SC SNP q values not correlated. 
#A and SC SNP q values not correlated.
	





#Come back to this later, there is somethign here.

# I repeated the same analysis for significant FST genes between pops (above for Sel vs Con). 
		#Are the saem genes and same SNPs high FST between these pops?
			#out of laziness, I'm going to use the same SNPs......
SC.pfst = PFST[c(1,4,5,6)]; names(SC.pfst)[4]="SCq"
M.pfst = read.table(file="A.counts",header=F)
M.pfst$q = qvalue0(M.pfst$V3)$q		
M.pfst$SNP = paste(M.pfst$V1, M.pfst$V2, sep=":")
M.pfst =  merge(M.pfst, SC.pfst, by="SNP")	
sigFST.M = M.pfst[M.pfst$q < 0.01,]
Genic.Fst.M = c()
chrom = intersect(gff$V1,sigFST.M$CHROM)
for(i in chrom){
	win.temp = gff[gff$V1==i,]
	deg.temp = sigFST.M[which(as.character(sigFST.M$CHROM)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$V3)), ">=") 
	blah1=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$V4)), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	Genic.Fst.M = rbind(temp,Genic.Fst.M)
	print(i)
}
names(Genic.Fst.M)[10] ="GB"
high.genes.M = as.character(unique(Genic.Fst.M$GB))	
sum(high.genes %in% high.genes.C)
#[1] 571 !!!! WHAT!~ 88% of high.genes is in this list......
sum(high.genes %in% high.genes.A) #high.genes.X contain the genes with high FST outliers in pop-level comparisons.
#573
sum(high.genes %in% high.genes.M) 
#559

	#check against my permuted SNPs sigFSTperm
high.genes.perm = sigFSTperm$GB[sigFSTperm$GB %in% nsyn$V8]
sum(high.genes.perm %in% high.genes.C)
#523
sum(high.genes.perm %in% high.genes.A) 
#522
sum(high.genes.perm %in% high.genes.M) 	
#511		
		
		
test = merge(Genic.Fst.M, nsyn.1, by="SNP")
high.genes.M = as.character(unique(test$GB))			
		
		
		
		
		
		
		
		
		
		
		
		
		



# How informative are the significant, permuted SNPs within control population?
testSNPs = Genic.Fst[which(Genic.Fst$GB %in% sigFSTperm$GB),]
write.list(testSNPs$SNP[!duplicated(testSNPs$SNP)],file="PermutedSIGFST")
	#Create a ped and map file of these SNPs
vcftools --vcf /media/data1/forty3/drone/vcf_drone/DroneSelectionFinal.recode.vcf --positions /media/data1/forty3/drone/FST/pFST/PermutedSIGFST --plink --out /media/data1/forty3/drone/FST/pFST/PERMUTEDFST --keep /media/data1/forty3/drone/vcf_drone/controlBees.txt




map=read.table(file="PERMUTEDFST.map")
dt = read.table("PERMUTEDFST.ped", colClasses=c(rep("character",16376)))
pheno=read.table(file="/media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt",header=F)
pheno=pheno[,-c(2)]
vec=merge(pheno,dt,by="V1")
vec=vec[,-c(3:7)]
names(vec)[2]="pheno"
nms=vec$V1
vec=vec[,-1]
phens=vec$pheno
vec=vec[,-1]



c1=vec[,seq(1,ncol(vec),2)]
c2=vec[,seq(2,ncol(vec),2)]

pfun <- function(x, y) paste(x, y, sep = "")
datNew <- vector("list", ncol(c1)/2)
for (i in 1:ncol(c1)) {
    datNew[[i]] <- pfun(c1[[i]], c2[[i]])
}

vec <- as.data.frame(datNew);names(vec)=seq(1:ncol(vec))
vec$pheno=phens
names(vec)=paste("SNP",names(vec),sep="")

lev = apply(vec, 2, function(x) length(table(x)))
vec = vec[,-c(which(lev=="1"))]

lm(SNPpheno~., data=vec)











####
#2) RF #need to optimize....
vec.rand = vec[,(sample(1:13188,1000))]
vec.rand$pheno = vec$SNPpheno

remo=ncol(vec.rand)-1
bestmtry = tuneRF(vec.rand[-remo], vec.rand$pheno,  ntreeTry=1000, stepFactor=0.5,improve=0.01, trace=TRUE, plot=FALSE, dobest=FALSE)
bt=data.frame(bestmtry)
mtry=as.numeric(bt[,1][which(bt[,2]==min(bt[,2]))])

RFfit <- randomForest(pheno ~ .,  data=vec.rand, importance=TRUE, ntree=10000,mtry=333, proximity=T,localImp=T)









#GLM to see how much these SNPs explain variation





library("rpart")


fit <- rpart(pheno~., data=vec, method="class", control=rpart.control(minsplit=2, minbucket=1, cp=0.01))
plot(fit);text(fit)

numvec=vec
id <- c(1:ncol(numvec))
numvec[,id] <- as.numeric(as.character(unlist(numvec[,id])))
fit <- rpart(pheno~., data=numvec, method="anova", control=rpart.control(minsplit=2, minbucket=1, cp=0.001))
x11();plot(fit);text(fit)
pfit<- prune(fit, cp=fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])
plot(pfit)
text(pfit)
summary(pfit)

	
















			
		
		
		
	
		
		
		
		
		
		
		
	
	