###
#NGM Aligner for drone sequence
###




#Environment variables -------------------------------
export REF=/data2/reference



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
gatk -R $APIS/am45placed.fasta  -T UnifiedGenotyper \
	-I bams.list  \
	-o out.het.raw.vcf  \
	-stand_call_conf 60.0 \
	-stand_emit_conf 40.0 \
	-dcov 200 \
	--min_base_quality_score 20  \
	-nt 20 -glm SNP  \
	-ploidy 2  
 
 
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
