DroneSelectionFinal.recode.vcf
##File Creation
#1. Create PED and MAP files
cd /media/data1/forty3/drone/vcf_drone
vcftools --vcf DroneSelectionFinal.recode.vcf --plink --out AllSample

#2. Convert scaff to chrom
Rscript /media/data1/afz/git/ScaffMaptoChr.r AllSample.map

#3. Re-order
plink --noweb --file AllSample --recode  --out AllSamplere 


#4. Trim out based on r
plink --noweb --file AllSamplere --indep 50 5 2  
plink --noweb --file AllSamplere --extract plink.prune.in --recode --out AllSamplereINDEP
 
#5. Output 1 file per Chromosome
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; 
do plink --noweb --file AllSamplereINDEP  --recode --out DroneSamps_$K --chr $K ; done
 

#6. Run Shapeit
for K in  16 2 3 4 5 6 7 8 9 10 11 12 13 14 15 1; 
do ./shapeit -P /media/data1/forty3/drone/vcf_drone/DroneSamps_$K -T 2 -O/media/data1/forty3/drone/vcf_drone/DroneSamps_$K.phased ; done
 
cd /media/data1/afz/shapeit
./shapeit -P /media/data1/forty3/drone/vcf_drone/DroneSamps_5 -T 2 -O/media/data1/forty3/drone/vcf_drone/DroneSamps_5.phased 
 

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
nms =read.table(file="clipboard",header=T)



for(i in c(1:16)){
samp=read.table(file=paste("DroneSamps_",i,".phased.sample",sep=""),header=T)
samp=samp[-1,];samp=samp[-3]
samp=rbind(samp,samp);samp=samp[order(samp[,1]),]
haps=read.table(file=paste("DroneSamps_",i,".phased.haps",sep=""),header=F)
maps=haps[c(2,1,3)];
haps=haps[-c(1,2,3,4,5)]
names(haps)=paste(nms$ID_1,nms$grp,sep="_")
haps=as.matrix(haps)
FR = apply(haps[,c(59:82)],1,function(x) length(x[x=="0"])/length(x)) #
#Where FR > 0.5, all 0's become 2 and all 1's stay 1
fr.haps=haps[which(FR >= 0.5),]
fr.haps[which(fr.haps=="0")] = "2"
haps[which(FR >= 0.5),] = fr.haps
#Where FR<0.5 all 1's become 2's and 0's become 1's
fr.haps=haps[which(FR < 0.5),]
fr.haps[which(fr.haps=="1")] = "2"
fr.haps[which(fr.haps=="0")] = "1"
haps[which(FR < 0.5),] = fr.haps

#Major Allele in Control population will be "ancestral"
maps$anc= rep("2", nrow(haps)) ;maps$rep=rep("1", nrow(haps))
chaps = haps[,grep("_C",colnames(haps))] #control pop
shaps = haps[,grep("_S",colnames(haps))] #selected pops
#Output dataframes
thaps = t(chaps)
thaps = data.frame(cbind(seq(1,nrow(thaps)),thaps))
write.list(thaps, file=paste("DroneSampsPH_",i,"_con.haps",sep=""))
write.list(maps, file=paste("DroneSampsPH_",i,".inp",sep=""))
thaps = t(shaps)
thaps = data.frame(cbind(seq(1,nrow(thaps)),thaps))
write.list(thaps,file=paste("DroneSampsPH_",i,"_sel.haps",sep=""))


#selected population
data.sel<-data2haplohh(hap_file=paste("DroneSampsPH_",i,"_sel.haps",sep=""),map_file=paste("DroneSampsPH_",i,".inp",sep=""))
res.sel<-scan_hh(data.sel)
ihs.cgu<-ihh2ihs(res.sel)
test = data.frame(ihs.cgu$res.ihs)
test$SNP = row.names(test)
write.list(test[test$Pvalue<0.05,],file="IHS.out",append=T)
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




























haps=read.table(file="DroneSamps_6.phased.haps ")


6. Add in Allelic information into 
for K in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 ; \
do python  /home/sani/task1/Phased_Script_PP.py DroneSamps_16.phased.haps ; done


plink --noweb --file DroneSamps_16.phased.haps --recode  --out test --chr $K




