# Workflow and Supplementary Material for Drone Selection Experiment






This project utilizes Illumina Sequence data from a selection experiment for hygienic behaviour. It selected bees for three generations using either a field assay for hygiene (called FAS population) or a metric of field assay and expression of marker proteins (called MAS population; see [Guarna et al. 2015](http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-014-1193-6)). These populations were maintained along with a mase-line, unselected population (BM). Each bee's population-origin is listed in DroneSamps.txt


To do:
Most of this is still (guiltily) hard-coded, so I need to go through and generalize it. 



##VCF Creation
I aligned 2 different data sets using NGM. First, all Drones individually and second, I merged drones into a single bam file where each merged file contained the ~3 drones sequenced per queen. Each fastq was trimmed with Trimmomatic v0.32 e.g:
<pre><code>java -jar /usr/share/java/trimmomatic-0.32.jar PE -threads 30 -phred33 -trimlog 3870-3.trimlog 3870-3_R1.fastq 3870-3_R2.fastq 3870-3_R1_TP.fastq 3870-3_R1_TU.fastq 3870-3_R2_TP.fastq 3870-3_R2_TU.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 </code></pre>
I next followed [GATK's Best Practices for Alignment](https://www.broadinstitute.org/gatk/guide/bp_step.php)  and aligned with NGM (NGMDrone.sh), removed duplicate reads, and re-aligned around indels. 
I hard-filtered sites within 10bp of indels and sites within 5bp of putative CNVs (sites that were called hetero in my haploid drones). This was performed in (trimdrone.sh). This script also trims QD < 5.0 || FS > 40.0 || MQ < 25.0 and removes SNPs with outlier Depth and Quality scores (VCFQualityDepthFilter.r)

The resulting VCF file is Drone.Hap.recode.vcf and is found in /vcf_drone


##SNP Functional Classification
I used SNPEFF to calssify mutations putative functional roles.
<pre><code>java -jar /usr/local/lib/snpEff2/snpEff.jar Amel -o txt Drone.Hap.recode.vcf -no-downstream -no-upstream  > HYG.snpeff.eff</code></pre>
<pre><code>java -jar /usr/local/lib/snpEff2/snpEff.jar Amel -o txt Drone.Hap.recode.vcf   > HYG-up_dwn.snpeff.eff</code></pre>

I repeated this for my high candate sites
<pre><code>java -jar /usr/local/lib/snpEff2/snpEff.jar Amel -o txt candidates.recode.vcf -no-downstream -no-upstream  > CAND.snpeff.eff</code></pre>
<pre><code>java -jar /usr/local/lib/snpEff2/snpEff.jar Amel -o txt candidates.recode.vcf   > CAND-up_dwn.snpeff.eff</code></pre>
	
These data are found in /vcf_drone

		
# Drone Analysis Pipeline

#Selection Analysis

##Hapflk 
I ran hapflk and summarized it using  hapflksummary.r. This script will concatenate all the hapflk outputs and also extract broad differentiated ranges (within 50Kb). The output is saved into "hapFLK50kb.RData"

For outlier analysis, I made use of the recent [hapflk software](https://forge-dga.jouy.inra.fr/projects/hapflk). This analysis was carried out using hapflk.sh, hapflksummary.r, and hapFLKPlot.r. Those scripts will run hapflk on a .vcf file, ouput phased .ped files and then run hapflk, summarizing the data in .RData files. I identified outliers as any site with Q > 0.01. 
		
When calculating Fst, I used [pFst and wcFst](https://github.com/jewmanchue/vcflib/wiki/Association-testing-with-GPAT) straight from final VCF files. 

<pre><code>
./pFst --target 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,17,18,19,20,21,22,23,24,25,26,27,28 \
	   --background 30,31,32,33,34,35,36,37,38,39,40 \
	   --deltaaf 0.05 \
	   --file /media/data1/forty3/drone/vcf_drone/Drone.Hap.recode.vcf \
	   --counts --type PL   > /media/data1/forty3/drone/FST/pFST/pFST.out
</code></pre>

<pre><code>
./wcFst --target 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28 \
	   --background 29,30,31,32,33,34,35,36,37,38,39,40 \
	   --type PL \
	   --file /media/data1/forty3/drone/vcf_drone/Drone.Hap.recode.vcf \
	    >/media/data1/forty3/drone/FST/pFST/wcFST.out 
</code></pre>

These data are found in /working


##Tajima's D
I ran this in 1000 bp windows using VCFTOOLS and only on the selected populations. The output is saved in S.Tajima.D. This is summarized with TDsummary.r and save in TD.RData. The script also performs a running median of 251 windows. 


##iHs
I calculated [iHS](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0040072) statistic across each chromosome within the selected population. Each allele, after phasing, was deemed ancestral or derived based on it's frequency in the baseline population. The script is REHH.sh and will run across chromosomes and compile the results for iHS into a single outfile (iHS.out). I used the [REHH package in R](https://cran.r-project.org/web/packages/rehh/vignettes/rehh.pdf)


##CSS
The composite measure takes in all the selection measures I used and combines them in a non-parametric way. I've implement CSS as per the [original article](http://bmcgenet.biomedcentral.com/articles/10.1186/1471-2156-15-34) using css.r, kindly provided by the first author. 

	
##HWE
Ran HWE test with vcftools, had to make some alterations.
		
		
###Population identity 
Here, I use the SNPs our group [previously identified](http://www.pnas.org/content/111/7/2614.abstract) in A, M, and C lineages and compare them to the selected populations to identify where the selected alleles originated from.

1. Create tabix index for merging
<pre><code>bgzip ALLSNP.recode.vcf
tabix -p vcf ALLSNP.recode.vcf.gz
bgzip Drone.Hap.recode.vcf
tabix -p vcf Drone.Hap.recode.vcf.gz
</code></pre>

2. Merge VCFs
<pre><code>vcf-merge Drone.Hap.recode.vcf.gz /media/data1/forty3/brock/align/ALLSNP.recode.vcf.gz | bgzip -c > HYG.vcf.gz
HYG.vcf.gz</code></pre>


####LAMP
I used LAMP.sh and LAMPanalysis.r across phased chromosomes for both selected and control populations together. I extracted this (LAMPanalysis.r) to find evidence for differential admixture at hygienic loci. I did the same thing with ADMIXTURE and by using A, M, and C major alleles. 

These data are found in /vcf_drone/lamp


####MAF and HWE
To look at allele frequencies within selected regions of the genome. Followed this up with anal.r
<pre><code>
vcftools --vcf DroneSelection.vcf --remove controlBees.txt --freq --max-alleles 2 --out S
vcftools --vcf Drone.Hap.recode.vcf --hardy --remove controlBees.txt --out S
</code></pre>


##DEGs and QTLs
There are lists of [DEGs](http://www.biomedcentral.com/1471-2164/16/500), [DEPs](http://www.biomedcentral.com/1471-2164/16/63), Some [GWAS](http://journals1.scholarsportal.info/pdf/1755098x/v12i0002/323_doa4sabihbmc.xml) and QTL papers ([1](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-294X.2010.04569.x/full) and [2](http://link.springer.com/article/10.1007/s00114-002-0371-6#page-1)) available for hygienic behaviour. I pulled these data in to see if I've got evidence of significant FST SNPs within them. I would expect some within QTLs, but not necessarily any in DEGs. 

#### QTLs and GWAS
I used Oxley's QTLs and mapped their location in AMEL4.5 by BLAST'ing the location of the nearest SNP within said QTL's peak LOD score. These are the best matches I could find using Solignac's marker set (scaffold:pos format).

|	Amel4	|	Start	|	End	|
|	-------------	|	:-------------:	|	-------------:	|
|	hyg1	|	2.19:1245006	|	2.19:86589	|
|	hyg2	|	5.14:527544	|	5.14:1558442	|
|	hyg3	|	16.2:42885	|	16.4:920721	|
NOTE: This has been further refined and updated. See GenomeR/QTLs dataset.

For the "GWAS" I pulled Spotter's putative, unpublished QTLs from Table 1:

|	LG	|	Start	|	End	|
|	---------	|	:---------:	|	---------	|
|	LG1	|	3039231	|	8453574	|
|	LG1	|	9418717	|	16819942	|
|	LG2	|	1	|	12503099	|
|	LG6	|	11206828	|	17739083	|
|	LG7	|	9515998	|	12848973	|
|	LG12	|	1	|	4003353	|
|	LG13	|	5247545	|	10266737	|
|	LG15	|	1	|	6643609	|
|	LG16	|	3196393	|	6242592	|

	
#### DEGs and DEPs
LF kindly provided significant DEPs. I pulled Boutin's Tables 4 and 5. 
With these lists I looked for genes with significant FST SNPs (HYGFSTAnalyses.r and boutinDEGs.r). I permuted SNPs across the genome for signiciance. 



<!--- 

vcftools --vcf Drone.Hap.recode.vcf --keep /media/data1/forty3/drone/git/data/FASBees.txt --window-pi 1000  --out FAS
vcftools --vcf Drone.Hap.recode.vcf --remove allbutMAS.txt --TajimaD 1000  --out MAS





##FST Analyses
(cd /media/data1/forty3/drone/FST/pFST/vcflib/bin)
Used [pFst and wcFst](https://github.com/jewmanchue/vcflib/wiki/Association-testing-with-GPAT) to estimate pairwise Fst and p-values between selected (pooled) and control:

<pre><code>
./pFst --target 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,17,18,19,20,21,22,23,24,25,26,27,28 \
	   --background 30,31,32,33,34,35,36,37,38,39,40 \
	   --deltaaf 0.05 \
	   --file /media/data1/forty3/drone/vcf_drone/Drone.Hap.recode.vcf \
	   --counts --type PL   > /media/data1/forty3/drone/FST/pFST/pFST.out
</code></pre>

<pre><code>
./wcFst --target 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28 \
	   --background 29,30,31,32,33,34,35,36,37,38,39,40 \
	   --type PL \
	   --file /media/data1/forty3/drone/vcf_drone/Drone.Hap.recode.vcf \
	    >/media/data1/forty3/drone/FST/pFST/wcFST.out 
</code></pre>

I then re-ran the analysis with only FAS against BM 

<pre><code>
./pFst --target 14,17,18,19,20,21,22,23,24,25,26,27,28 \
	   --background 30,31,32,33,34,35,36,37,38,39,40 \
	   --deltaaf 0.05 \
	   --file /media/data1/forty3/drone/vcf_drone/Drone.Hap.recode.vcf \
	   --counts --type PL   > /media/data1/forty3/drone/FST/pFST/pFST.FAS.out
</code></pre>

Run again, but only for MAS against BM

<pre><code>
./pFst --target 0,1,2,3,4,5,6,7,8,9,10,11,12,13 \
	   --background 30,31,32,33,34,35,36,37,38,39,40 \
	   --deltaaf 0.05 \
	   --file /media/data1/forty3/drone/vcf_drone/Drone.Hap.recode.vcf \
	   --counts --type PL   > /media/data1/forty3/drone/FST/pFST/pFST.MAS.out
</code></pre>



###Fst between Lineages and SEL/CON
To get an idea of where the alleles came from (kind of), I'm using FST between the major lineages
<pre><code>
vcftools --vcf HYG.vcf --weir-fst-pop /media/data1/forty3/drone/git/data/m.txt --weir-fst-pop /media/data1/forty3/drone/git/data/SEL.txt --maf 0.05 --out M_vs_SEL
vcftools --vcf HYG.vcf --weir-fst-pop /media/data1/forty3/drone/git/data/c.txt --weir-fst-pop /media/data1/forty3/drone/git/data/SEL.txt --maf 0.05  --out C_vs_SEL
vcftools --vcf HYG.vcf --weir-fst-pop /media/data1/forty3/drone/git/data/s.txt --weir-fst-pop /media/data1/forty3/drone/git/data/SEL.txt --maf 0.05  --out A_vs_SEL
vcftools --vcf HYG.vcf --weir-fst-pop /media/data1/forty3/drone/git/data/m.txt --weir-fst-pop /media/data1/forty3/drone/git/data/CONT.txt --maf 0.05  --out M_vs_CONT
vcftools --vcf HYG.vcf --weir-fst-pop /media/data1/forty3/drone/git/data/c.txt --weir-fst-pop /media/data1/forty3/drone/git/data/CONT.txt --maf 0.05  --out C_vs_CONT
vcftools --vcf HYG.vcf --weir-fst-pop /media/data1/forty3/drone/git/data/s.txt --weir-fst-pop /media/data1/forty3/drone/git/data/CONT.txt --maf 0.05  --out A_vs_CONT
vcftools --vcf HYG.vcf --weir-fst-pop /media/data1/forty3/drone/git/data/m.txt --weir-fst-pop /media/data1/forty3/drone/git/data/c.txt --maf 0.05  --out M_vs_C
vcftools --vcf HYG.vcf --weir-fst-pop /media/data1/forty3/drone/git/data/s.txt --weir-fst-pop /media/data1/forty3/drone/git/data/c.txt --maf 0.05  --out A_vs_C
vcftools --vcf HYG.vcf --weir-fst-pop /media/data1/forty3/drone/git/data/m.txt --weir-fst-pop /media/data1/forty3/drone/git/data/s.txt --maf 0.05  --out M_vs_A
vcftools --vcf HYG.vcf --keep /media/data1/forty3/drone/git/data/m.txt --TajimaD 1000 --maf 0.05  --out M
vcftools --vcf HYG.vcf --keep /media/data1/forty3/drone/git/data/c.txt --TajimaD 1000 --maf 0.05  --out C
vcftools --vcf /media/data1/forty3/brock/align/N.raw.vcf  --TajimaD 1000 --maf 0.05  --out N
</code></pre>


###Output High FST regions and plots
I munged the fst data using DroneAnalysis.r. This script takes in the outputs above, merges them, creates unique SNP IDs, and processed it into NCBI chromosomes. The latter is performed by a perl script developed by Amro Zayed and slightly modified by me (scaffold_to_chr.pl). Once prociessed into chromosomes, I run a creeping window average across the genome in 1 and 5 kb windows using the [Qanbari et al. 2012](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0049525) approach from my own [scripts](https://github.com/harpur/GenomeR). It outputs "ClusteredHighSNPsCreeper5kb", a list of high FST regions, and "RAWOUT.RData". 


##Extract high SNPs and characterize
"high.bed" contains all SNPs in AND codition for 1kb windows. Using this, I extracted the SNPs from my VCF file
<pre><code>
vcftools --vcf Drone.Hap.recode.vcf --bed high.bed --recode --out high
java -jar /usr/local/lib/snpEff2/snpEff.jar Amel -o txt high.recode.vcf -no-downstream -no-upstream  > HYG.high.snpeff.eff

</code></pre>

I characterized SNPs within high windows, high snps (FSTP >  4) relative to all SNPs within the genome. I used -ud 5000 

##Analysis
All is within DroneAnalysis.r
All analyses saved in "Hygiene.RData"




####Admixture Analysis

vcftools --vcf HYG.vcf --plink --remove /media/data1/forty3/drone/git/data/o.txt
plink --file out --noweb --thin 0.25 --maf 0.05 --make-bed 

for K in 1 2 3 4 5 ; \
do /home/brock/admixture/admixture  --cv=20 plink.bed $K -j12 | tee log${K}.out; done
grep -h CV log*.out 

for K in 4 ; \
do /home/brock/admixture/admixture  --cv=20 plink.bed $K -j12 | tee log${K}.out; done
mv plink.4.Q plink.4.Q.wholegenome

#trying now with O


vcftools --vcf HYG.vcf --bed highregions.bed --plink  --out HYGHIGHHAPS #--remove /media/data1/forty3/drone/git/data/o.txt
plink --file HYGHIGHHAPS --noweb  --make-bed 


for K in 4 ; \
do /home/brock/admixture/admixture  --cv=20 plink.bed $K -j7 | tee log${K}.out; done
mv plink.4.Q plink.4.Q.selectedsites



#testing with highest sites and NO A 

vcftools --vcf HYG.vcf --bed VERYhigh.bed --plink  --out VHYGHIGHHAPS --remove /media/data1/forty3/drone/git/data/o.txt --remove /media/data1/forty3/drone/git/data/s.txt
plink --file VHYGHIGHHAPS --noweb  --make-bed 


for K in 2 ; \
do /home/brock/admixture/admixture  --cv=20 plink.bed $K -j7 | tee log${K}.out; done
mv plink.2.Q plink.2.Q.Vselectedsites











###GO Analysis
To come



###TRN Analysis




###Gene Age



-->




























































<!--- 





#### GWAS 
I've stopped this analysis because I lack the power to perform GWAS associations
Using SNPs within the highest clusters, I ran 2 GWAS statistics to see if there is an association between genotype and phenotype. The samples we have are stratified, so I controlled for this by using [Permutation and a CMH test](http://pngu.mgh.harvard.edu/~purcell/plink/anal.shtml). The permutation used the quantitative phenotype and then CMH used a categorical case-control where case was defined as having >80% hygiene.

#####File Creation
1. Create PED and MAP files
<pre><code>cd /media/data1/forty3/drone/vcf_drone</code></pre>
<pre><code>vcftools --vcf HYG.vcf --plink --vcf HYG.vcf --bed VERYhigh.bed --out AllSample</code></pre>

2. Convert scaff to chrom
<pre><code>Rscript /media/data1/forty3/drone/git/ScaffMaptoChr.r AllSample.map</code></pre>

3. Re-order
<pre><code>plink --noweb --file AllSample --recode  --out AllSamplere</code></pre>


4. Trim out based on r
<pre><code>plink --noweb --file AllSamplere --indep-pairwise 50 5 0.5 </code></pre> 
<pre><code>plink --noweb \
	--file AllSamplere   \
	--extract plink.prune.in \
	--recode \
	--out AllSamplereINDEP
</code></pre> 



#####Association

<pre><code>plink --noweb --file AllSamplereINDEP --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHBcat.txt --mh2 --within SELCONcluster.txt --out  CLUSTEREDUNPCAT </code></pre>

<pre><code>plink --noweb --file  AllSamplereINDEP --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --assoc --adjust --perm --out PERMUN </code></pre>

















load(file="HYGRESULTS_FIGS.RDATA")


con = read.table(file="/media/data1/forty3/drone/vcf_drone/con.hwe",header=T, colClasses = c("character", "numeric", "character","character","numeric","numeric"))
sel = read.table(file="/media/data1/forty3/drone/vcf_drone/sel.hwe",header=T, colClasses = c("character", "numeric", "character","character","numeric","numeric"))

sel = con
sel$SNP = paste(sel$CHR, sel$POS, sep="_")
sel$GRP = paste("Group", sel$CHR, sep="")




# Pi vs FST ----------------------------------------
Pi.Fst = c()
chrom = intersect(regions$GRP, sel$GRP)
for(i in chrom){
	win.temp = regions[regions$GRP==i,]
	deg.temp = sel[which(as.character(sel$GRP)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$GPOS)), ">=") 
	blah1=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$end)), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	Pi.Fst = rbind(temp,Pi.Fst)
	print(i)
}


Pi.Fst[Pi.Fst$GRP=="Group6.2",]




6.2 479112              3/5/3   2.75/5.50/2.75 0.090909 1.000000
284019 6.2 480028              2/3/6   1.11/4.77/5.11 1.517551 0.232198
284020 6.2 480220              1/4/6   0.82/4.36/5.82 0.076389 1.000000
284021 6.2 480238   

CHR	POS	OBS(HOM1/HET/HOM2)	E(HOM1/HET/HOM2)	ChiSq	P
1.1	838	10/1/0	10.02/0.95/0.02	0.024943	1.000000
1.1	940	10/1/0	10.02/0.95/0.02	0.024943	1.000000



might update this for LROH.....


#proportion of homoz. in each region?



hom1 = as.numeric(gsub("[/]","",substr(Pi.Fst$OBS.HOM1.HET.HOM2.,1,2)))
hom2 = as.numeric(gsub("[/]","",substr(Pi.Fst$OBS.HOM1.HET.HOM2.,5, nchar(Pi.Fst$OBS.HOM1.HET.HOM2.))))

hom.prop = hom1 + hom2

#SEL:
length(hom.prop[hom.prop>20]) #(707;0.91)

#sel all genome:

sel = sel[which(!(sel$GRP %in% regions$GRP)),]

hom1 = as.numeric(gsub("[/]","",substr(sel$OBS.HOM1.HET.HOM2.,1,2)))
hom2 = as.numeric(gsub("[/]","",substr(sel$OBS.HOM1.HET.HOM2.,5, nchar(sel$OBS.HOM1.HET.HOM2.))))

hom.prop = hom1 + hom2

length(hom.prop[hom.prop>20]) #(707;0.91)


#CON
length(hom.prop[hom.prop>8]) #656


	Outcome 1	Outcome 2	     Total
Group 1	707	67	774
Group 2	656	118	774
Total	1363	185	1548
Fisher's exact test
  The two-tailed P value is less than 0.0001

  more homozygous loci in these regions.


substr(Pi.Fst$OBS.HOM1.HET.HOM2.,1, -3)



Pi.Fst = c()
chrom = intersect(regions$GRP, sel$GRP)
for(i in chrom){
	win.temp = regions[regions$GRP==i,]
	deg.temp = sel[which(as.character(sel$GRP)==as.character(i)),]
	blah=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$GPOS)), "<=") 
	blah1=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$end)), ">=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	Pi.Fst = rbind(temp,Pi.Fst)
	print(i)
}








##Output PLINK format
1. Convert to PLINK
<pre><code>vcftools --vcf Drone.Hap.recode.vcf --recode --plink</code></pre>

2. Convert scaff to chrom
<pre><code>Rscript /media/data1/forty3/drone/git/ScaffMaptoChr.r out.map drone</code></pre>

3. Re-order
<pre><code>plink --noweb --file drone.chrom --recode  --out drone.chrom</code></pre> 

4. Add Phenotype information
<pre><code>
cut -d " " -f1-6  drone.chrom.ped > first6
cut -d " " -f1-6 --complement drone.chrom.ped > drone.ped
paste -d " " first6 drone.ped > drone.test.ped 
rm drone.chrom.ped
rm drone.ped
mv drone.test.ped drone.chrom.ped
</code></pre> 






###compare to NA bees
<pre><code>
gatk -R /home/amel45/AM45/am45new.fasta -T UnifiedGenotyper \
	-I bams.list  \
	-o out.NA.raw.vcf  \
	-stand_call_conf 60.0 \
	-stand_emit_conf 40.0 \
	-dcov 200 \
	--min_base_quality_score 20  \
	-nt 16 -glm SNP  \
	-ploidy 2 &
	
vcftools --vcf out.NA.raw.vcf --max-alleles 2 --freq
 </code></pre>
 

I ran a quick-and-dirty SNP call of Ontarioan bees to see what proportion of my candidate SNPs are actually present in this poopulation as I'll be genotyping in them
 
 
 

##Nucleotide Diversity and HWE
vcftools --vcf Drone.Hap.recode.vcf --window-pi 5000 --remove controlBees.txt --out sel &
vcftools --vcf Drone.Hap.recode.vcf --window-pi 5000 --keep controlBees.txt --out con &
vcftools --vcf Drone.Hap.recode.vcf --window-pi 5000 --keep pop3.txt --out p3 & #sel pop
vcftools --vcf Drone.Hap.recode.vcf --window-pi 5000 --keep pop1.txt --out p1 & #sel pop


vcftools --vcf Drone.Hap.recode.vcf --hardy --remove controlBees.txt --out sel &
vcftools --vcf Drone.Hap.recode.vcf --hardy --keep controlBees.txt --out con &
vcftools --vcf Drone.Hap.recode.vcf --hardy --keep pop3.txt --out p3 & #sel pop
vcftools --vcf Drone.Hap.recode.vcf --hardy --keep pop1.txt --out p1 & #sel pop


vcftools --vcf DroneSelection.vcf --get-INFO MQ
vcftools --vcf DroneSelection.vcf --get-INFO MQ0 --out mq0

 
 
 vcftools --vcf Drone.Hap.recode.vcf --site-mean-depth --remove controlBees.txt --out sel
 
 
 --site-mean-depth
 
 
 ###Output HWE for each population
<pre><code>
vcftools --vcf Drone.Hap.recode.vcf --max-alleles 2 --freq
</code></pre>
 
 
 
 
 
 
 
 
 


##Admixture Analyses
I used the final VCF file, paired with SNPs called from Harpur et al. 2014 to run ADMIXTURE. 


* Get out common SNPs between AMC and AHB in a merged VCF file, create a .ped file. For this, I copied AMC.ped and AMC.map from my AFZ project. This contains SNPS for AMC with MAF 0.05. 

WARNING: H-scroll and block of code...sorry :)
<pre><code>
vcftools --vcf /media/data1/afz/VCF/AllAMCSNPs.recode.vcf --max-alleles 2 --plink  --out AMC
vcftools --vcf /media/data1/forty3/drone/vcf_drone/DroneSelectionFinal.recode.vcf --max-alleles 2 --plink  --out HYG
Rscript /media/data1/afz/git/intersectingMap.r AMC.map HYG.map 
vcftools --vcf /media/data1/forty3/drone/vcf_drone/DroneSelectionFinal.recode.vcf --positions Shared.map --recode --out HYG
vcftools --vcf /media/data1/afz/VCF/AllAMCSNPs.recode.vcf --positions Shared.map --recode --out AMC
gatk -T CombineVariants -R /home/amel45/AM45/am45new.fasta  --variant HYG.recode.vcf --variant AMC.recode.vcf -o HYGmergedAMC.vcf -genotypeMergeOptions REQUIRE_UNIQUE
vcftools --vcf HYGmergedAMC.vcf --max-alleles 2 --plink  --out AMCHYG
vcftools --vcf HYGmergedAMC.vcf --recode --positions HIGHFST.map --plink  --out AMCHYGHIGH
</code></pre>


* Use the .ped files to create .bim, .fam,  and .bed for ADMIXTURE. I did this for all significant FST SNPs and for a random selection of SNPs. Then, run ADMIXTURE for K=3.

<pre><code>
plink --file AMCHYGHIGH --noweb --make-bed --out AMCHYGHIGH
/home/brock/admixture/admixture  --cv=10 AMCHYGHIGH.bed 3 -j2 | tee log3.out

plink --file AMCHYGHIGH --noweb --make-bed --out AMCHYGHIGH
/home/brock/admixture/admixture  --cv=10 AMCHYGHIGH.bed 3 -j2 | tee log3.out
</code></pre>


for K in 1 2 3 4 5; \
do /home/brock/admixture/admixture  --cv=10 AMCHYGHIGH.bed $K -j20 | tee log${K}.out; done
grep -h CV log*.out


* I used this to look at differences in introgression between selected and control lines at significant SNPs


saved in HygieneHighFSTADMIXTURE.xlsx


##Association Analysis

###File Creation
1. Create PED and MAP files
<pre><code>cd /media/data1/forty3/drone/vcf_drone</code></pre>
<pre><code>vcftools --vcf DroneSelectionFinal.recode.vcf --plink --out AllSample</code></pre>

2. Convert scaff to chrom
<pre><code>Rscript /media/data1/afz/git/ScaffMaptoChr.r AllSample.map</code></pre>

3. Re-order
<pre><code>plink --noweb --file AllSample --recode  --out AllSamplere</code></pre> 

4. Trim out based on r
<pre><code>plink --noweb --file AllSamplere --indep 50 5 2 </code></pre> 
<pre><code>plink --noweb \
	--file AllSamplere   \
	--extract plink.prune.in \
	--recode \
	--out AllSamplereINDEP
</code></pre> 

5. Output 1 file per Chromosome
<pre><code>
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
do plink --noweb --file AllSamplereINDEP  --recode --out DroneSamps_$K --chr $K ; done
</code></pre> 


##If Phased Data wanted:
6. Run Shapeit
<pre><code>
for K in  16 2 3 4 5 6 7 8 9 10 11 12 13 14 15 1; \
do ./shapeit -P /media/data1/forty3/drone/vcf_drone/DroneSamps_$K -T 2 -O/media/data1/forty3/drone/vcf_drone/DroneSamps_$K.phased ; done
</code></pre> 

7. Add in Allelic information into 
<pre><code>
for K in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
do python  /media/data1/forty3/drone/git/Phased_Script_PP.py DroneSamps_$K.phased.haps ; done
</code></pre> 

8. Write out phased, tped files for association
<pre><code>
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
	do Rscript HapstoTPED.r /media/data1/forty3/drone/vcf_drone/DroneSamps_$K.phased.haps.out; done
</code></pre> 

8a. Frequency removal
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
	do plink --noweb --tfile DroneSamps_$K.phased.sample --maf 0.05  --recode --out DroneSamps_$K.phasedMAF ; done



###Associations with Cochran-Mantel-Haenszel (CMH) tests with unphased data
<pre><code>
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
	do (plink --noweb --file DroneSamps_$K.UNphasedMAF --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --mh2 --within SELCONcluster.txt --out  CLUSTEREDUNPHMAF2_$K ) & 
	done
</code></pre> 

###Associations with Permutation tests with unphased data
<pre><code>
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
	do (plink --noweb --file DroneSamps_$K.UNphasedMAF --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --assoc --perm --out PERMUNhap_$K ) & 
	done
</code></pre> 


###Compile results from above
Here, I focussed on phased outputs



























##DEGs and QTLs
There are lists of [DEGs](http://www.biomedcentral.com/1471-2164/16/500), [DEPs](http://www.biomedcentral.com/1471-2164/16/63), Some [GWAS](http://journals1.scholarsportal.info/pdf/1755098x/v12i0002/323_doa4sabihbmc.xml) and QTL papers ([1](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-294X.2010.04569.x/full) and [2](http://link.springer.com/article/10.1007/s00114-002-0371-6#page-1)) available for hygienic behaviour. I pulled these data in to see if I've got evidence of significant FST SNPs within them. I would expect some within QTLs, but not necessarily any in DEGs. 

#### QTLs and GWAS
I used Oxley's QTLs and mapped their location in AMEL4.5 by BLAST'ing the location of the nearest SNP within said QTL's peak LOD score. These are the best matches I could find using Solignac's marker set (scaffold:pos format).

|	Amel4	|	Start	|	End	|
|	-------------	|	:-------------:	|	-------------:	|
|	hyg1	|	2.19:1245006	|	2.19:86589	|
|	hyg2	|	5.14:527544	|	5.14:1558442	|
|	hyg3	|	16.2:42885	|	16.4:920721	|


For the "GWAS" I pulled Spotter's putative, unpublished QTLs from Table 1:

|	LG	|	Start	|	End	|
|	---------	|	:---------:	|	---------	|
|	LG1	|	3039231	|	8453574	|
|	LG1	|	9418717	|	16819942	|
|	LG2	|	1	|	12503099	|
|	LG6	|	11206828	|	17739083	|
|	LG7	|	9515998	|	12848973	|
|	LG12	|	1	|	4003353	|
|	LG13	|	5247545	|	10266737	|
|	LG15	|	1	|	6643609	|
|	LG16	|	3196393	|	6242592	|

	

#### DEGs and DEPs
LF kindly provided significant DEPs. I pulled Boutin's Tables 4 and 5. 
With these lists I looked for genes with significant FST SNPs (HYGFSTAnalyses.r and boutinDEGs.r). I permuted SNPs across the genome for signiciance. 



###GO Analysis
I've focussed solely on significant FST genes for the time. I used GOstats (BostatsBEE.r) with a gene universe composed of fly orthologs to honey bee genes. For the first test, I used any-old gene with a significant SNP. For the second analysis, I used only genes with NSYN SNPs that had significantly more significant SNPs (that's fun to say) than expected by chance. My permutation procedure can be seen in HYGFSTAnalyses.r (~line 182). 


###TL;DR
Check in HYGFSTAnalyses.r






###Selection over longer time frames?
I'm going to pull out the list of significant genes and see if they have evidence of selection over longer times within Apis (Gamma), and between populations (Pi, TD, Fst).



###I'D LIKE TO TRY EXEHH http://hgdp.uchicago.edu/Software/



###Plotting Data
All FST plot scripts can be found as .r files

* FST Histogram, FST by chromosome (and stupid legend) = FSTHygienPlot.r

* Boxplot of admixture in FST SNPs = HighFSTAdmixPlot.r

* I ran REVIGO using the web app (no R version :( ). Uploaded REVIGO_AHighFST.r and REVIGO_AllNSYNHighFST.r for REVIGO plots of all sig FST and all sig NSYN Fst. 

* I may re-run Rcircos, see Rcircos_Drone.r


--->













