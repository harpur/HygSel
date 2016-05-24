#trim SNPs

#Environment variables -------------------------------
export REF=/data2/reference


 	#REMOVE 17s!!!!
#Output hetero SNPs --------------------
vcftools --vcf out.het.raw.vcf --hardy 
sed -i '/^17[.]/ d' out.hwe 
sed -i '/^18[.]/ d' out.hwe 
sed -i '/[/]0[/]/ d' out.hwe 
sed -i '/[/]1[/]/ d' out.hwe  
sed -i '/[/]2[/]/ d' out.hwe 
sed -i '/[/]3[/]/ d' out.hwe 
sed -i '/[/]4[/]/ d' out.hwe 
cut -f 1,2 out.hwe  > putativeCNV.list 
grep -Fwf putativeCNV.list  out.raw.vcf > cnv.vcf #for each species, I created one cnv.vcf.
sed '/^#/ !d' out.raw.vcf > VCFheader
cat VCFheader cnv.vcf > cnv1.vcf
rm cnv.vcf
mv cnv1.vcf cnv.vcf

#putativeCNV.list merged and take all unique cases from each species to created BombusCNV.list


# Filter Variants  -----------------------------------
gatk -R $REF/am45placed.fasta -T VariantFiltration \
	-V out.raw.vcf \
	--filterExpression "QD < 5.0 || FS > 40.0 || MQ < 25.0 "   \
	--filterName "mpileFilters" \
	--mask out.raw.indels.vcf  \
	--maskExtension 10 \
	--maskName "InDel" -o  out.vcf 

gatk -R $REF/am45placed.fasta -T VariantFiltration \
	-V out.vcf \
	--mask cnv.vcf \
	--maskExtension 5 \
	--maskName "cnv" \
	-o  out2.vcf

rm out.vcf 
mv out2.vcf out.vcf 


cat out.vcf | vcf-convert -v 4.1 > Ouput.vcf
vcftools --vcf Ouput.vcf --recode --remove-filtered-all --max-alleles 2 --out out.indel #will be called out.indel.recode.vcf


# Filter Depth and Quality outliers  -----------------------------------
Rscript VCFQualityDepthFilter.r "out.indel.recode.vcf" 


# Output Missingness  -----------------------------------
vcftools --vcf out.indel.dp.q.recode.vcf --missing 








R
site.miss = read.table(file="out.lmiss",header=T, colClass=c("character", "numeric","numeric","numeric","numeric","numeric")) 
ind.miss = read.table(file="out.imiss",header=T) #no individual had 5% missing loci. Good.
write.table(site.miss[c(1,2)][which(site.miss$F_MISS>0.05),], file="missing.loci", col.names=F,row.names=F,quote=F)

vcftools --vcf out.indel.dp.q.recode.vcf --exclude-positions missing.loci --recode --out out.indel.dp.q.miss
	#1328938 total SNPs

cp out.indel.dp.q.miss.recode.vcf /mnt/nfs/data1/apis/forty3/drone/vcf_drone/




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




