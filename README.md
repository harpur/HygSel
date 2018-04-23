# Workflow and Supplementary Material for Drone Selection Experiment






This project utilizes Illumina Sequence data from a selection experiment for hygienic behaviour. It selected bees for three generations using either a field assay for hygiene (called FAS population) or a metric of field assay and expression of marker proteins (called MAS population; see [Guarna et al. 2015](http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-014-1193-6)). These populations were maintained along with a base-line, unselected population (BM). Each bee's population-origin is listed in DroneSamps.txt



## VCF Creation
I aligned 2 different data sets using NGM. First, all Drones individually and second, I merged drones into a single bam file where each merged file contained the ~3 drones sequenced per queen. Each fastq was trimmed with Trimmomatic v0.32 e.g:
<pre><code>java -jar /usr/share/java/trimmomatic-0.32.jar PE -threads 30 -phred33 -trimlog 3870-3.trimlog 3870-3_R1.fastq 3870-3_R2.fastq 3870-3_R1_TP.fastq 3870-3_R1_TU.fastq 3870-3_R2_TP.fastq 3870-3_R2_TU.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 </code></pre>
I next followed [GATK's Best Practices for Alignment](https://www.broadinstitute.org/gatk/guide/bp_step.php)  and aligned with NGM (NGMDrone.sh), removed duplicate reads, and re-aligned around indels. 
I hard-filtered sites within 10bp of indels and sites within 5bp of putative CNVs (sites that were called hetero in my haploid drones). This was performed in (trimdrone.sh). This script also trims QD < 5.0 || FS > 40.0 || MQ < 25.0 and removes SNPs with outlier Depth and Quality scores (VCFQualityDepthFilter.r)

<!---
The resulting VCF file is Drone.Hap.recode.vcf and is found in /vcf_drone
-->

## SNP Functional Classification
I used SNPEFF to calssify mutations putative functional roles.
<pre><code>java -jar /usr/local/lib/snpEff2/snpEff.jar Amel -o txt Drone.Hap.recode.vcf -no-downstream -no-upstream  > HYG.snpeff.eff</code></pre>
<!---
<pre><code>java -jar /usr/local/lib/snpEff2/snpEff.jar Amel -o txt Drone.Hap.recode.vcf   > HYG-up_dwn.snpeff.eff</code></pre>
I repeated this for my high candate sites
<pre><code>java -jar /usr/local/lib/snpEff2/snpEff.jar Amel -o txt candidates.recode.vcf -no-downstream -no-upstream  > CAND.snpeff.eff</code></pre>
-->
<!---	
These data are found in /vcf_drone
-->
		
# Drone Analysis Pipeline

# Selection Analysis

## Hapflk 
I ran hapflk and summarized it using  hapflksummary.r. This script will concatenate all the hapflk outputs and also extract broad differentiated ranges (within 50Kb). 

<!---	
The output is saved into "hapFLK50kb.RData"
-->

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
<!---	
These data are found in /working
-->

## Tajima's D
I ran this in 1000 bp windows using VCFTOOLS and only on the selected populations. The output is saved in S.Tajima.D. This is summarized with TDsummary.r and save in TD.RData. The script also performs a running median of 251 windows. 


## iHs
I calculated [iHS](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0040072) statistic across each chromosome within the selected population. Each allele, after phasing, was deemed ancestral or derived based on it's frequency in the baseline population. The script is REHH.sh and will run across chromosomes and compile the results for iHS into a single outfile (iHS.out). I used the [REHH package in R](https://cran.r-project.org/web/packages/rehh/vignettes/rehh.pdf)


## CSS
The composite measure takes in all the selection measures I used and combines them in a non-parametric way. I've implement CSS as per the [original article](http://bmcgenet.biomedcentral.com/articles/10.1186/1471-2156-15-34) using css.r, kindly provided by the first author. 

	
### Population identity 
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


#### LAMP
I used LAMP.sh and LAMPanalysis.r across phased chromosomes for both selected and control populations together. I extracted this (LAMPanalysis.r) to find evidence for differential admixture at hygienic loci. I did the same thing with ADMIXTURE and by using A, M, and C major alleles. 
<!---	
These data are found in /vcf_drone/lamp
-->

#### MAF and HWE
To look at allele frequencies within selected regions of the genome. Followed this up with anal.r
<pre><code>
vcftools --vcf DroneSelection.vcf --remove controlBees.txt --freq --max-alleles 2 --out S
vcftools --vcf Drone.Hap.recode.vcf --hardy --remove controlBees.txt --out S
</code></pre>


## QTLs
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





