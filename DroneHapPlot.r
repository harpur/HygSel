###
# Plot 6.2 
###





source("/media/data1/forty3/drone/git/GenomeR/VCFFunctions.r")

sorter <- function(x){
	srts = sort(summary(as.factor(x)), decreasing=T)
	if(length(srts)<2){
		return(names(sort(summary(as.factor(x)), decreasing=T)[1]))
	}else{
		if(sum(srts[1] == srts[2])){
		return("NA")
		}else{
			return(names(sort(summary(as.factor(x)), decreasing=T)[1]))
		}
}
}

samps = read.table(file="/media/data1/forty3/drone/git/data/DroneSamps.txt")

system("vcftools --vcf out.indel.dp.q.miss.recode.vcf --bed VERYhigh.bed --recode --out HYGHIGHHAPS" )
i = "HYGHIGHHAPS.recode.vcf" 
vcf = Read.VCF()

system("vcftools --vcf /media/data1/forty3/brock/align/AllBees_SNPs.raw.vcf --bed VERYhigh.bed --recode --out ALL")
i = "ALL.recode.vcf" 
C.vcf = Read.VCF() 

system("vcftools --vcf /media/data1/forty3/brock/align/N.raw.vcf --bed VERYhigh.bed --recode --out NA")
i = "NA.recode.vcf" 
N.vcf = Read.VCF() 


#####################
#choose random SNPs
#####################

samps = read.table(file="/media/data1/forty3/drone/git/data/DroneSamps.txt")

system("cut -f 1,2 out.indel.dp.q.miss.recode.vcf > sites")
system("shuf -n 10000 sites > output")


system("vcftools --vcf out.indel.dp.q.miss.recode.vcf --positions output --recode --out HYGHIGHHAPS" )
i = "HYGHIGHHAPS.recode.vcf" 
vcf = Read.VCF()

system("vcftools --vcf /media/data1/forty3/brock/align/AllBees_SNPs.raw.vcf --positions output --recode --out ALL")
i = "ALL.recode.vcf" 
C.vcf = Read.VCF() 

system("vcftools --vcf /media/data1/forty3/brock/align/N.raw.vcf --positions output --recode --out NA")
i = "NA.recode.vcf" 
N.vcf = Read.VCF() 


############################







#Extract Alleles fro each Haploid Drone ---------------------
	#this is, for now, hard-coded
pos = sapply(vcf,function(x) return((x[2])))
pos = pos[-1]
chr = sapply(vcf,function(x) return((x[1])))[-1]
snp = paste(chr, pos, sep="_")
allele  = sapply(vcf,function(x) substr(unlist(x)[10:134],1 ,1 )) 
allele = allele[,-1]

maj = apply(allele[c(93:125),], 2, function(x) sorter(x)) #major allele in control population 
s.maj = apply(allele[c(1:92),], 2, function(x) sorter(x)) #major allele in selected population 
#maj = matrix(maj, nrow(allele), nc = ncol(allele), byrow = T)
#z = maj == allele
#allele = z


#Extract Alleles fro each C population ---------------------
	#this is, for now, hard-coded
C.pos=sapply(C.vcf,function(x) return((x[2])))
C.pos = C.pos[-1]
C.chr = sapply(C.vcf,function(x) return((x[1])))[-1]
C.snp = paste(C.chr, C.pos, sep="_")



#from WhichAlleleMin for c pop
one = sapply(C.vcf,function(x) length( grep("1/1",unlist(x)[c(10:18)])))[-1]
het  = sapply(C.vcf,function(x) length( grep("0/1",unlist(x)[c(10:18)])))[-1]
zer = sapply(C.vcf,function(x) length( grep("0/0",unlist(x)[c(10:18)])))[-1]
C.siz = (one + het + zer) 
C.maj = as.numeric(zer < one) #if it's the 0 allele, FALSE(0) is alt allele, TRUE(1)
C.min = as.numeric(!C.maj)

#from WhichAlleleMin for M pop
one = sapply(C.vcf,function(x) length( grep("1/1",unlist(x)[c(19:27)])))[-1]
het  = sapply(C.vcf,function(x) length( grep("0/1",unlist(x)[c(19:27)])))[-1]
zer = sapply(C.vcf,function(x) length( grep("0/0",unlist(x)[c(19:27)])))[-1]
M.siz = (one + het + zer) 
M.maj = as.numeric(zer < one) #if it's the 0 allele, FALSE(0) is alt allele, TRUE(1)
M.min = as.numeric(!M.maj)

#from WhichAlleleMin for NA pop
one = sapply(C.vcf,function(x) length( grep("1/1",unlist(x)[c(10:14)])))[-1]
het  = sapply(C.vcf,function(x) length( grep("0/1",unlist(x)[c(10:14)])))[-1]
zer = sapply(C.vcf,function(x) length( grep("0/0",unlist(x)[c(10:14)])))[-1]
N.siz = (one + het + zer) 
N.maj = as.numeric(zer < one) #if it's the 0 allele, FALSE(0) is alt allele, TRUE(1)
N.min = as.numeric(!M.maj)



#from WhichAlleleMin for c pop
one = sapply(C.vcf,function(x) length( grep("1/1",unlist(x)[c(38:48)])))[-1]
het  = sapply(C.vcf,function(x) length( grep("0/1",unlist(x)[c(38:48)])))[-1]
zer = sapply(C.vcf,function(x) length( grep("0/0",unlist(x)[c(38:48)])))[-1]
A.siz = (one + het + zer) 
A.maj = as.numeric(zer < one) #if it's the 0 allele, FALSE(0) is alt allele, TRUE(1)
A.min = as.numeric(!A.maj)


MC = data.frame(C.snp, C.maj, M.maj, A.maj, N.maj); names(MC)[1] = "snp"
MC$m.siz = M.siz
MC$c.siz = C.siz
MC$a.siz = A.siz


#MC = MC[which(MC$m.siz>0 & MC$c.siz>0 & MC$a.siz>0),]
MC = MC[which(MC$C.maj != MC$M.maj),]

#		#for absolute M alleles:
#		MC = MC[which(MC$C.maj == MC$A.maj & MC$A.maj != MC$M.maj),]


#Which alleles in CONT population are major alleles in C lineage?
CON = data.frame(snp, maj, s.maj)
test = merge(CON, MC, by = "snp")





#Take all sites in "allele" that I have in "test"	
shared.allele = allele[,which(snp %in% test$snp)] #so now test contains ONLY C lienage alleles. 
snp = snp[which(snp %in% test$snp)] 


majeC = test$C.maj #1 is C lineage allele in CON
majeC  = matrix(majeC , nrow(shared.allele), nc = ncol(shared.allele), byrow = T)
z = majeC == shared.allele
allele = z

allele[allele=="TRUE"] = "grey" #grey is where Cmaj = pop
allele[allele!="grey"] = "black" #black is Cmaj /= pop


#Show proportion of Cmin alleles fixed in each pop
scaff = as.character(gsub("_.*","",snp))
pos = as.numeric(gsub(".*_","",snp))


con = c()
sel = c()
scf = c()
for(i in unique(scaff)){
#this counts the number of black and grey counds, from above, in this case, grey is C lineage
	
	#con = c(table(allele[(c(93:125)),which(scaff==i)])/sum(table(allele[(c(93:125)),which(scaff==i)])),con) #C
	con = c(length(allele[c(93:125),which(scaff==i)][allele[c(93:125),which(scaff==i)]=="black"])/length(allele[c(93:125),which(scaff==i)]),con) #C
	
	#sel = c(table(allele[(c(1:92)),which(scaff==i)])/sum(table(allele[(c(1:92)),which(scaff==i)])), sel) #S
	sel = c(length(allele[c(1:92),which(scaff==i)][allele[c(1:92),which(scaff==i)]=="black"])/length(allele[c(1:92),which(scaff==i)]),sel)
	
	scf = c(i,scf)

}

#con = data.frame("C.nCmaj" = con[seq(1,length(con),2)], "C.Cmaj" = con[seq(2,length(con),2)] )
#sel = data.frame("nmaj" = sel[seq(1,length(sel),2)], "Cmaj" = sel[seq(2,length(sel),2)] )
#sel$scaff = scf
#sel$chr = gsub("[.].*","",sel$scaff)
#sel$C.Cmaj = con$C.Cmaj


sel = data.frame(cbind(con = 1-con, sel = 1-sel))
sel$scaff = scf
sel$chr = gsub("[.].*","",sel$scaff)



###
#M.def = aggregate((sel$sel-sel$con),by=list(sel$chr), mean) #average excess deficit of M allele in SEL
#M.def$se = aggregate((sel$sel-sel$con),by=list(sel$chr), function(x) sd(x)/sqrt(length(x)))$x #average excess deficit of M allele in SE
#M.def$chr = rep("S", nrow(M.def))

def = aggregate((sel$sel-sel$con),by=list(sel$chr), mean) #average excess deficit of M allele in SEL
def$se = aggregate((sel$sel-sel$con),by=list(sel$chr), function(x) sd(x)/sqrt(length(x)))$x #average excess deficit of M allele in SE
def$chr = rep("R", nrow(M.def))

M.def = rbind(def, M.def)


















pos = 1:length(snp)
#png(file="MlineageAlleles.png")
plot(pos,y=rep(1,length(pos)), col=allele[1,],pch=15, ylim = c(1,126), ylab='', yaxt='n', xlab = "Position (bp)")
for(i in 2:125){
	points(pos,y=rep(i,length(pos)), col=allele[i,],pch=15)
}
abline(h = 88, col="red", lty=2)













##ASSOC w/in control pop?

assoc = z#[(c(93:125)),]
nams = unlist(vcf[1])[c(10:134)]
row.names(allele) = nams
nams = gsub("-.*","",nams)
#nams = nams[c(93:125)]


geno = matrix(nc = ncol(z)+1, nr = length(unique(nams)))
for(i in unique(nams)){
	#i = "3780"
	vec = apply(z[c(which(nams == i)),],2, function(x) names(which.max(table(x)))) #DOMINANCE EFFECT!!
	geno[(which(unique(nams)==i)),] = c(i, vec)

	#vec = apply(z[c(which(nams == i)),],2, function(x) sum(x)) #geno EFFECT!!
	#vec[vec=="2"] =1 
	#geno[(which(unique(nams)==i)),] = c(i, vec)



}
#1 = M major alleles' 0 = C




geno = data.frame(geno); names(geno)[1] = "V2"
names(geno)[2:ncol(geno)] = paste("SNP",snp,sep="")
geno = merge(samps, geno, by ="V2")
geno.C = geno[geno$V3=="C",]

nums = apply(geno.C,2,function(x) length(unique(x)))
geno.C = geno.C[,which(nums>1)]


ps=c()
for(i in 4:ncol(geno.C)){
	p = summary(aov(geno.C$V4~unlist(geno.C[i])))[[1]][["Pr(>F)"]]
	ps = c( ps , p[1] )
}


sig.CM.snps = names(geno.C[c(4:ncol(geno.C))])[which(ps<0.05)]
sig.snps = gsub("SNP","",sig.snps)

NP11.18_1354137

summary(aov(geno.C$V4~geno.C$SNP6.2_510951))




##############

summary(aov(geno.C$V4~geno.C$SNP14.1_6398))
                    Df Sum Sq Mean Sq F value Pr(>F)
geno.C$SNP14.1_6398  2  705.1   352.5   4.257   0.05 *
Residuals            9  745.3    82.8
---
Signif. codes:  0 â***â 0.001 â**â 0.01 â*â 0.05 â.â 0.1 â â 1
> TukeyHSD(aov(geno.C$V4~geno.C$SNP14.1_6398))
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = geno.C$V4 ~ geno.C$SNP14.1_6398)

$`geno.C$SNP14.1_6398`
          diff       lwr         upr     p adj
1-0 -12.573179 -31.12832  5.98196562 0.1963551
3-0 -17.094365 -34.13835 -0.05037634 0.0493606
3-1  -4.521186 -23.92662 14.88424916 0.7967472



#create barplot df
melted = melt(geno.C, id.vars=c("SNP14.1_6398"),measure.vars=c("V4"))
means.sem = ddply(melted, c("SNP14.1_6398"), summarise,
	mean=mean(value,na.rm=T), 
	sem=sd(value,na.rm=T)/sqrt(length(value)))
means.sem = transform(means.sem, lower=mean-sem, upper=mean+sem)

names(means.sem)[1] = "Genotype"
means.sem$Genotype = c("C/C", "M/C", "M/M")



#plot 
p<-ggplot(means.sem, aes(x = factor(Genotype), y = mean,  fill = factor(Genotype))) +  
	geom_bar(position = "dodge",stat="identity",color = "black") + 
	coord_cartesian(ylim = c(0, 95)) + 
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
	x= "Genotype in Control Population",
	y = "Hygienic Performance(% Cells Removed)") 
	#scale_y_continuous(expand = c(0,0)) +
	#scale_x_discrete(expand = c(0.1,0))	+
	#geom_segment(aes(x = 1, y = .6, xend = 2, yend = .6), linetype=3) +
	#geom_segment(aes(x = 1, y = .6, xend = 1, yend = .59), linetype=3) +
	#geom_segment(aes(x = 2, y = .6, xend = 2, yend = .59), linetype=3) +
	#annotate("text", x = 1.5,  y = 0.62,label = "***" )
	
p











dev.off()
#ok, so if the abpve is correct--there is no difference in C lineage alleles in this particular site

#I re-ran with M alleles --big difference!


#plot with ordered hygiene values:
nams = unlist(vcf[1])[c(10:134)]
row.names(allele) = nams
nams = gsub("-.*","",nams)
nams =  data.frame(nams);names(nams)="V2"
nams = merge(nams, samps, by ="V2")

allele2 = (cbind(allele, nams))
allele2 = allele2[with(allele2, order(-nams$V4)), ]
allele2 = allele2[-c(719:722)]

###


#grey is now C alleles and black is M (major)
plot(pos,y=rep(1,length(pos)), col=as.character(allele[1,]),pch=15, ylim = c(1,126))
for(i in 2:125){
	points(pos,y=rep(i,length(pos)), col=allele[i,],pch=15)

}





pos = pos[1:43]

plot(pos,y=rep(1,length(pos)), col=allele[1,],pch=15, ylim = c(1,92))
for(i in 2:92){
	points(pos,y=rep(i,length(pos)), col=allele[i,],pch=15)

}



x11();plot(pos,y=rep(1,length(pos)), col=allele[93,],pch=15, ylim = c(92,126))
for(i in 94:125){
	points(pos,y=rep(i,length(pos)), col=allele[i,],pch=15)

}


plot(pos,y=rep(1,length(pos)), col=allele[1,],pch=15, ylim = c(1,126))
for(i in 2:125){
	points(pos,y=rep(i,length(pos)), col=allele[i,],pch=15)

}





#3173, 3507, 3521
#3822, 3824, 3780

x11();plot(pos,y=rep(1,length(pos)), col=allele[38,],pch=15, ylim = c(1,126))
for(i in c(39,40,55,54,53,77:80,102:104, 105:107,93:95, 120:122)){
	points(pos,y=rep(i,length(pos)), col=allele[i,],pch=15)

}

0	3145	S	81.18286409
1	3152	S	83.38084232
2	3154	S	97.66517635
3	3155	S	84.96376812
4	3156	S	90.18342391
5	3160	S	86.46111854
6	3162	S	88.28682674
7	3163	S	98.1544665
8	3164	S	75
9	3165	S	89.87549841
10	3169	S	74.12430759
11	3172	S	98.10606061

12	3173	S	100 38-40

13	3175	S	97.05957439
14	3498	S	97.79411765
15	3504	U	82.6
16	3505	U	100

17	3507	S	100 53 -55

18	3508	S	91.05833963
19	3510	S	92.3348064
20	3516	S	69.8326432
21	3517	S	98.52941176
22	3518	S	99.609375
23	3519	S	99.59677419
24	3520	S	98.82041406
25	3521	S	100 - 77-80
26	3522	S	97.11748634
27	3523	S	92.02256944
28	3525	S	94.57983193
29	3774	C	77.6
30	3780	C	58.59901089 -93
31	3800	C	75.92039801
32	3803	C	59.38177445
33	3822	C	50.85504106 - 102-104
34	3824	C	53.37600547 -105
35	3828	C	74.43722322
36	3841	C	73.97938351
37	3858	C	66.28032456
38	3862	C	67.86173667
39	3867	C	86.18573061 - 120:122
40	3870	C	83.00438596






 [1]"3145-1" "3145-2" "3145-3" "3152-1" "3152-2" "3152-3" "3154-1"
 [8] "3154-2" "3154-3" "3155-1" "3155-2" "3155-3" "3156-1" "3156-2" "3156-3"
 [16] "3160-1" "3160-2" "3160-3" "3162-1" "3162-2" "3162-3" "3163-1" "3163-2"
 [24] "3163-3" "3164-1" "3164-2" "3164-3" "3165-1" "3165-2" "3165-3" "3169-1"
 [32] "3169-2" "3169-3" "3169-4" "3172-1" "3172-2" "3172-3" "3173-1" "3173-2"
 [40] "3173-3" "3175-1" "3175-2" "3175-3" "3498-1" "3498-2" "3498-3" "3504-1"
 [48] "3504-2" "3504-3" "3505-1" "3505-2" "3505-3" "3507-1" "3507-2" "3507-3"
 [56] "3508-1" "3508-2" "3508-3" "3510-1" "3510-2" "3510-3" "3516-1" "3516-2"
 [64] "3516-3" "3517-1" "3517-2" "3517-3" "3518-1" "3518-2" "3518-3" "3519-1"
 [72] "3519-2" "3519-3" "3520-1" "3520-2" "3520-3" "3521-1" "3521-2" "3521-3"
 [80] "3521-4" "3522-1" "3522-2" "3522-3" "3523-1" "3523-2" "3523-3" "3525-1"
 [88] "3525-2" "3525-3" "3774-1" "3774-2" "3774-3" 
 
[93]"3780-1" "3780-2" "3780-3"
[96] "3800-1" "3800-2" "3800-3" "3803-1" "3803-2" "3803-3" "3822-1" "3822-2"
[104] "3822-3" "3824-1" "3824-2" "3824-3" "3828-1" "3828-2" "3828-3" "3841-1"
[112] "3841-2" "3841-3" "3858-1" "3858-3" "3862-1" "3862-2" "3862-3" "3862-4"
[120] "3867-1" "3867-2" "3867-3" "3870-1" "3870-2" "3870-3"

[c(93:125),]






