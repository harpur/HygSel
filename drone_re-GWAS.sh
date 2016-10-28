###
# test indels against pheno - redoing all the pipelin EXCEPT indel filtration
### 








gatk \
   -T VariantsToBinaryPed \
   -R /home/tanu/02_GATK_pipeline_testing/am45new.fasta  \
   -V drone_raw_INDEL.vcf \
   -m FSTDroneSelectionCandidatesALL.fam \
   -bed output.bed \
   -bim output.bim \
   --minGenotypeQuality 25 \
   -fam output.fam

bim = read.table(file="output_indel.bim",header=F, colClasses = rep("character", 6))
bim$V5 = rep("T", nrow(bim))
bim$V6 = rep("G", nrow(bim))
bim$V3 = rep("0", nrow(bim))
bim$V1 = paste("Group", bim$V1, sep="")

out = bim[-grep("Group1[7-8].", bim$V1),]
out2 = bim[grep("Group1[7-8].", bim$V1),]



write.table(out[c(2,1,4)],file="/media/data1/forty3/drone/ScaffConvert/CandidateMap",col.names=F,row.names=F,quote=F)
system("perl /media/data1/forty3/drone/ScaffConvert/scaffold_to_chr_BAH.pl /media/data1/forty3/drone/ScaffConvert/scaffolds_on_chr.txt /media/data1/forty3/drone/ScaffConvert/CandidateMap")

#map = read.table(file=fil,header=T)
conv = read.table(file="/media/data1/forty3/drone/ScaffConvert/CandidateMap_onChr.txt",header=F)
conv$V2 = gsub("chr","",conv$V2)

out$V1 = conv$V2
out$V4 = conv$V3
out2$V1 = rep("17", nrow(out2))

out = rbind(out, out2)


write.table(out, "output.bim",row.names=F,col.names=F,quote=F)


plink --noweb --bfile output --recode  --out out 

grep ^17 out.map | cut -f2 > chr17

plink --noweb --file out  --exclude chr17 --remove /media/data1/forty3/drone/vcf_drone/controlBees.txt --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --assoc --adjust --qq-plot --qt-means --out indels


ass = read.table("indels.qassoc.adjusted",header=T)


plot(-log10(ass$QQ), -log10(ass$UNADJ))
abline(0,1)

all = read.table(file="indels.qassoc",header=T)
plot(-log10(all$P), col = all$CHR, pch = 19, main = "INDEL")






#INDELS WORKS!!!!!! WTF!!! IDIOT!!!


vcftools --vcf /media/data1/forty3/drone/vcf_drone/DroneSelectionFinal.recode.vcf --plink --out outSNP


Rscript /media/data1/forty3/drone/git/ScaffMaptoChr.r outSNP.map

plink --noweb --file outSNP --recode  --out AllSamplere  

plink --noweb --file AllSamplere  --remove /media/data1/forty3/drone/vcf_drone/controlBees.txt --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --assoc --adjust --qq-plot


plot(-log10(ass$QQ), -log10(ass$UNADJ))
abline(0,1)




































###
# CHECK FST
###


plink --noweb --file out --remove /media/data1/forty3/drone/vcf_drone/controlBees.txt --freq
plink --noweb --file out --keep /media/data1/forty3/drone/vcf_drone/controlBees.txt --freq --out con


conf = read.table("con.frq",header=T)
self = read.table("plink.frq",header=T)


conf$sel = self$MAF


conf$diff = conf$MAF - conf$sel




#Output hetero SNPs --------------------
vcftools --vcf out.het.raw.vcf --hardy 
sed '/[/]0[/]/ d' out.hwe > test.out
sed '/[/]1[/]/ d' test.out > test1.out  
cut -f 1,2 test1.out > putativeCNV.list 
rm test.out
rm test1.out
grep -Fwf /data2/bombus/git/data/BombusCNV.list out.raw.vcf > cnv.vcf #for each species, I created one cnv.vcf.
sed '/^#/ !d' out.raw.vcf > VCFheader
cat VCFheader cnv.vcf > cnv1.vcf
rm cnv.vcf
mv cnv1.vcf cnv.vcf


cat out.het.raw.vcf | vcf-convert -v 4.1 > Ouput.vcf

sed -i '/^17[.]/ d' Ouput.vcf
sed -i '/^18[.]/ d' Ouput.vcf 


cat out.vcf | vcf-convert -v 4.1 > Ouput.vcf


# Filter Variants  -----------------------------------
	#gatk -R $REF/am45placed.fasta -T VariantFiltration \
	#	-V Ouput.vcf  \
	#	--filterExpression "QD < 5.0 || FS > 40.0 || MQ < 25.0 "   \
	#	--filterName "mpileFilters" \
	#	--mask out.raw.indels.vcf  \
	#	--maskExtension 10 \
	#	--maskName "InDel" -o  out.vcf 

gatk -R /home/amel45/AM45/am45new.fasta -T VariantFiltration \
	-V Output.vcf \
	--mask cnv.vcf \
	--maskExtension 5 \
	--maskName "cnv" \
	-o  out2.vcf

rm Ouput.vcf 
mv out2.vcf out.vcf 


cat out.vcf | vcf-convert -v 4.1 > Output.vcf
vcftools --vcf  out.vcf --recode --remove-filtered-all --max-alleles 2 --out out.cnv #will be called out.indel.recode.vcf



# Filter Depth and Quality outliers  -----------------------------------
Rscript VCFQualityDepthFilter.r "/media/data1/forty3/drone/bam_drone/out.cnv.recode.vcf" 


# Output Missingness  -----------------------------------
vcftools --vcf out.indel.dp.q.recode.vcf --missing 





R
site.miss = read.table(file="out.lmiss",header=T, colClass=c("character", "numeric","numeric","numeric","numeric","numeric")) 
ind.miss = read.table(file="out.imiss",header=T) #no individual had 5% missing loci. Good.
write.table(site.miss[c(1,2)][which(site.miss$F_MISS>0.05),], file="missing.loci", col.names=F,row.names=F,quote=F)

vcftools --vcf out.indel.dp.q.recode.vcf --exclude-positions missing.loci --recode --out out.indel.dp.q.miss
	#1328938 total SNPs

#cp out.indel.dp.q.miss.recode.vcf /mnt/nfs/data1/apis/forty3/drone/vcf_drone/





vcftools --vcf out.indel.dp.q.miss.recode.vcf --plink --out outSNP

Rscript /media/data1/forty3/drone/git/ScaffMaptoChr.r outSNP.map

mv outSNP.map outSNP_scaff.map
mv NA.chrom.map outSNP.map

plink --noweb --file outSNP --recode  --out AllSamplere  

plink --noweb --file   AllSamplere   --indep-pairwise 50 5 0.4 --maf 0.05
plink --noweb \
	--file AllSamplere   \
	--extract plink.prune.in \
	--recode \
	--out gwasINDEP2


plink --noweb --file gwasINDEP2  --remove /media/data1/forty3/drone/vcf_drone/controlBees.txt --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --assoc --adjust --qq-plot

ass = read.table("plink.qassoc.adjusted", header=T, nrow=100000)
plot(-log10(ass$QQ), -log10(ass$UNADJ))
abline(0,1)

ass = read.table("plink.qassoc", header=T)
plot(-log10(ass$P),col=ass$CHR, main = "SNPs", pch=19)







#original SP set?



vcftools --vcf DroneSelection.vcf --plink --out outDSNP

Rscript /media/data1/forty3/drone/git/ScaffMaptoChr.r outDSNP.map

mv outDSNP.map outDSNP_scaff.map
mv NA.chrom.map outDSNP.map

plink --noweb --file outDSNP --recode  --out AllSamplere  

#plink --noweb --file   AllSamplere   --indep-pairwise 50 5 0.4 --maf 0.05
#plink --noweb \
	#--file AllSamplere   \
	#--extract plink.prune.in \
	#--recode \
	#--out gwasINDEP2


plink --noweb --file AllSamplere  --remove /media/data1/forty3/drone/vcf_drone/controlBees.txt --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --assoc --adjust --qq-plot

ass = read.table("plink.qassoc.adjusted", header=T, nrow=100000)
plot(-log10(ass$QQ), -log10(ass$UNADJ))
abline(0,1)

ass = read.table("plink.qassoc", header=T)
plot(-log10(ass$P),col=ass$CHR, main = "SNPs_Trim", pch=19)



#Try the join statist?
















# Calling synonymous and non-synonymous sites -------------------------
	#note, v4.0 of SNPEFF no longer supports TXT format (I think)...annoying.
java -jar $SNPEFF/snpEff.jar Bimp \
	-o txt \
	out.indel.dp.q.recode.vcf \
	-no-downstream \
	-no-upstream > out.snpeff.eff
	
sed '/SYNONYMOUS/ !d' out.snpeff.eff > exons.eff #take only SNPS within exons from output file
sed '/WARNING_/ !d' exons.eff > warnexons.eff #take out warnings from output file
sed -i '/WARNING/ d' exons.eff 


# make VCF file of only SYN and NSYN SNPs -------------------------
cut -f 1,2 exons.eff  > synnsyn.list 
vcftools --vcf out.indel.dp.q.recode.vcf --positions synnsyn.list --recode --out out.nsyn

#create tabix index for merging ----------------------------------
bgzip out.nsyn.recode.vcf
tabix -p vcf out.nsyn.recode.vcf.gz



















vcftools --vcf Ouput.vcf  --plink --out outrawSNP


Rscript /media/data1/forty3/drone/git/ScaffMaptoChr.r outrawSNP.map

mv outrawSNP.map outrawSNP_scaff.map
mv NA.chrom.map outrawSNP.map

plink --noweb --file outrawSNP --recode  --out AllSamplere  

plink --noweb --file AllSamplere  --remove /media/data1/forty3/drone/vcf_drone/controlBees.txt --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --assoc --adjust --qq-plot --qt-means --out allSNPs


ass = read.table("allSNPs.qassoc.adjusted",header=T)


plot(-log10(ass$QQ), -log10(ass$UNADJ))
abline(0,1)

all = read.table(file="allSNPs.qassoc",header=T)
plot(-log10(all$P))

