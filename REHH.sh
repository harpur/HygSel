

export WRKDIR=/media/data1/forty3/drone/working





#create data set -------	
vcftools --vcf Drone.Hap.recode.vcf --plink --out outSNP
Rscript /media/data1/forty3/drone/git/ScaffMaptoChr.r outSNP.map
mv outSNP.map outSNP_scaff.map
mv NA.chrom.map outSNP.map
plink --noweb --file outSNP --recode  --out AllSamplere  


#4. Trim out based on r
plink --noweb --file AllSamplere --indep 50 5 2  
plink --noweb --file AllSamplere --extract plink.prune.in --recode --out AllSamplereINDEP
 
#5. Output 1 file per Chromosome
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; 
do plink --noweb --file AllSamplere  --recode --out DroneSamps_$K --chr $K ; done
 


#Phase---------------
#cd /media/data1/afz/shapeit/
#here, I used rho 0.39 (from Wallberg) across the whole genome
for K in  16 2 3 4 5 6 7 8 9 10 11 12 13 14 15 1; \
	do ./shapeit --rho 0.39 -P  $WRKDIR/DroneSamps_$K -T 2 --window 0.5 -O  $WRKDIR/DroneSamps_$K.phased ; done
	
	
#Prep files for EHH
	#DroneSamps_5.phased.sample
	#.hap and .inp
#.inp
#F0100190 1 113642 1 2
#F0100220 1 244699 1 2
#F0100250 1 369419 1 2
#F0100270 1 447278 1 2
#F0100280 1 487654 1 2
#F0100290 1 524507 2 1

#haps
#1 2 1 1 2 1...
#2 2 2 1 1 1...
#3 2 1 1 1 1...

library("rehh")
source("/media/data1/forty3/brock/scripts/GOgetter.r")
source("/media/data1/forty3/brock/scripts/VarFunct.r")
source("/media/data1/forty3/brock/scripts/movingavg.r")
nms =read.table(file="nms.txt",header=T)



for(i in c(7:16, 1:4)){
	samp=read.table(file=paste("DroneSamps_",i,".phased.sample",sep=""),header=T)
	samp=samp[-1,];samp=samp[-3]
	samp=rbind(samp,samp);samp=samp[order(samp[,1]),]
	haps=read.table(file=paste("DroneSamps_",i,".phased.haps",sep=""),header=F)
	maps=haps[c(2,1,3)];
	haps=haps[-c(1,2,3,4,5)]
	names(haps)=paste(nms$ID,nms$pop,sep="_")
	haps=as.matrix(haps)
	FR = apply(haps[,c(59:82)],1,function(x) length(x[x=="0"])/length(x)) #
	
	#https://cran.r-project.org/web/packages/rehh/vignettes/rehh.pdf
	#By default alleles are assumed to be coded as 0 (missing data), 1 (ancestral allele) or 2 (derived allele).
	#Recoding of the alleles in this format, according to the SNP information data file (see 1.2) can be performed
	#with the recode.allele option of the function data2haplohh() (see 1.3).
	
	#Where FR > 0.5, all 0's become 1 and all 1's become 2 
	fr.haps=haps[which(FR >= 0.5),]
	fr.haps[which(fr.haps=="1" )] = "2"
	fr.haps[which(fr.haps=="0" )] = "1"
	haps[which(FR >= 0.5),] = fr.haps
	
	#Where FR < 0.5 all 0's become 2's and 1's become 1's
	fr.haps=haps[which(FR < 0.5),]
	fr.haps[which(fr.haps=="1")] = "1"
	fr.haps[which(fr.haps=="0")] = "2"
	haps[which(FR < 0.5),] = fr.haps

	#Major Allele in Control population will be "ancestral"
	maps$anc= rep("2", nrow(haps)) ;maps$rep=rep("1", nrow(haps))
	chaps = haps[,grep("_BM",colnames(haps))] #control pop
	shaps = haps[,grep("_SEL",colnames(haps))] #selected pops
		
	#Output dataframes
	thaps = t(chaps)
	thaps = as.matrix(data.frame(cbind(seq(1,nrow(thaps)),thaps)))
	write.list(thaps, file=paste("DroneSampsPH_",i,"_con.haps",sep=""))
	write.list(maps, file=paste("DroneSampsPH_",i,".inp",sep=""))
	
	thaps = t(shaps)
	thaps = as.matrix(data.frame(cbind(seq(1,nrow(thaps)),thaps)))
	write.list(thaps,file=paste("DroneSampsPH_",i,"_sel.haps",sep=""))
	print(i)
}
	
for(i in c(16:2)){
	#selected population
	data.sel<-data2haplohh(hap_file=paste("DroneSampsPH_",i,"_sel.haps",sep=""),map_file=paste("DroneSampsPH_",i,".inp",sep=""))
	#data.con<-data2haplohh(hap_file=paste("DroneSampsPH_",i,"_con.haps",sep=""),map_file=paste("DroneSampsPH_",i,".inp",sep=""))
	res.sel<-scan_hh(data.sel,threads=10)
	#res.con<-scan_hh(data.con,threads=10)
	ihs.cgu<-ihh2ihs(res.sel, freqbin = 0.1)
	ihs = data.frame(ihs.cgu$iHS)
	ihs$rs = row.names(ihs)
	names(ihs)[4] ="ihs.p"
	write.list(ihs,file="IHS.out",append=T)
	print(i)
}
	
	#wg.rsb<-ies2rsb(res.sel,res.con)
	#names(wg.rsb)[4] ="p" #negative, pop2 bigger
		
	#wg.xpehh<-ies2xpehh(res.sel,res.con)
	#names(wg.xpehh)[4] ="p"


	#ihs.cgu<-ihh2ihs(res.sel, freqbin = 0.18)
	

	
	
	

	
	
	
}



#page 22 of manual
hap=data2haplohh(hap_file="DroneSampsPH_1_sel.haps",map_file="DroneSampsPH_1.inp")
#hap_chr_5.con
#hap_chr_1.sel 




test = data.frame(ihs.cgu$res.ihs)
#con
data.con<-data2haplohh(hap_file="DroneSampsPH_1_con.haps",map_file="DroneSampsPH_1.inp")
res.con<-scan_hh(data.con)

plot(test$POSITION, test$Pvalue, cex=0.2,pch=19)

#wg.rsb<-ies2rsb(res.con,res.sel,"CON","SEL")
#head(wg.rsb$res.rsb)


#do for only selected pop.
#http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0050171


# In ihh2ihs(res.sel) :
#  Size of Allele Frequency Class: 0.925-0.95 <10: You should probably increase freqbin













