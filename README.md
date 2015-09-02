# Workflow and Supplementary Material for Drone Selection Experiment









##VCF Creation:
* VCFCreation_DroneSelection.txt

* DroneVCFcreation.txt




##FST Analyses
###After Alignment, SNP calling, and trimming procedures (VCF Creation above)
* Fst between selected (N=2 pops or N=1) and unselected using DroneCalculateFSTfromVCF.sh (copy and paste)

* Run  MovingAveFSTfromPLINK.r <PLINK Fst> <WInsize> to create moving averages of FST along chromosomes (converted with updated perl script within R) saves an .RData file for next steps

* Get out high FST windows with HighWindowsFST.r <RData> <cutoff> (saved as dataframe called HighSNPs)

* If 2 independent FST analyses where used, get overlapping windows.

e.g.

<pre><code>Rscript MovingAveFSTfromPLINK.r 
	/media/data1/forty3/drone/vcf_drone/pop2_vs_sel.weir.fst 10000</code></pre>
<pre><code>Rscript HighWindowsFST.r \
	/media/data1/forty3/drone/vcf_drone/pop2_vs_sel.RDATA 0.95</code></pre>

save.image(file="/media/data1/forty3/drone/FST/SelvsCon/SelectedvsControlFSTPvals.RData")

* I also ran [pFST](https://github.com/jewmanchue/vcflib/wiki/Association-testing-with-GPAT) using pFST.sh and added this to the previous dataframe of FST values. Then ran FDR using Storey's Q method. Significant FST was based on q<0.01 (or, -log10(P)>2.61)


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

* I used this to look at differences in introgression between selected and control lines at significant SNPs

<!---
saved in HygieneHighFSTADMIXTURE.xlsx
-->

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



###Plotting Data
All FST plot scripts can be found as .r files

* FST Histogram, FST by chromosome (and stupid legend) = FSTHygienPlot.r

* Boxplot of admixture in FST SNPs = HighFSTAdmixPlot.r

* I ran REVIGO using the web app (no R version :( ). Uploaded REVIGO_AHighFST.r and REVIGO_AllNSYNHighFST.r for REVIGO plots of all sig FST and all sig NSYN Fst. 

* I may re-run Rcircos, see Rcircos_Drone.r
















