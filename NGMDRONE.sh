###
#NGM Aligner for drone sequence
###




#Environment variables -------------------------------
export REF=/data2/reference
export APIS=/home/amel45/AM45


#Align ----------------------
	#create input file "fastq", a list of all fastq files

#create .bam
filename='fastq'
exec 4<$filename
echo Start
while read -u4 p ; do
	pidlist_sampe=""
	endfq="_R1.fastq"
	endfq2="_R2.fastq"
	fq1=$p$endfq
	fq2=$p$endfq2
	ngm -r $REF/am45new.fasta  -1 $fq1 -2 $fq2 -o $p.bam -t 31 -p -b
done

#create .sorted.bam
filename='fastq'
exec 4<$filename
echo Start
while read -u4 p ; do
	pidlist_sampe=""
	endfq=".bam "
	fq1=$p$endfq
	samtools sort -m 6G -@ 12 $fq1 $p.sorted
	samtools index $p.sorted.bam
done



#Remove dups and add read groups ----------------------
filename='samples'
exec 4<$filename
echo Start
while read -u4 p ; do
	picard MarkDuplicates I=$p.sorted.bam  O=$p.sorted.dp.bam METRICS_FILE=$p.Dups VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
	picard AddOrReplaceReadGroups INPUT= $p.sorted.dp.bam   OUTPUT=$p.sorted.rg.dp.bam RGID=$p RGPL=illumina RGLB=$p RGPU=run RGSM=$p VALIDATION_STRINGENCY=LENIENT
	samtools index $p.sorted.rg.dp.bam
done




# Call SNPs with GATK -----------------------------------
ls *.sorted.rg.dp.bam > bams.list  #assemble list of bams in species directory


#Call SNPs
gatk -R $REF/am45placed.fasta -T UnifiedGenotyper \
	-I bams.list  \
	-o out.raw.vcf  \
	-stand_call_conf 60.0 \
	-stand_emit_conf 40.0 \
	-dcov 200 \
	--min_base_quality_score 20  \
	-nt 15 -glm SNP  \
	-ploidy 1 &
 
 
#Call INDELs
gatk -R $REF/am45placed.fasta -T UnifiedGenotyper \
	-I bams.list  \
	-o out.raw.indels.vcf  \
	-stand_call_conf 60.0 \
	-stand_emit_conf 40.0 \
	-dcov 200 \
	--min_base_quality_score 20  \
	-nt 25 -glm indel  \
	-ploidy 1 &
 
  
#Call SNPs, assume diploid
nohup gatk -R $APIS/am45new.fasta -T UnifiedGenotyper \
	-I bams.list  \
	-o out.het.raw.vcf  \
	-stand_call_conf 60.0 \
	-stand_emit_conf 40.0 \
	-dcov 200 \
	--min_base_quality_score 20  \
	-nt 18 -glm SNP  \
	-ploidy 2  &
 
 
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



















