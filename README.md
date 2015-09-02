# Workflow and Supplementary Material for Drone Selection Experiment









##VCF Creation:
- VCFCreation_DroneSelection.txt
- DroneVCFcreation.txt




##FST Analyses
###After Alignment, SNP calling, and trimming procedures (VCF Creation above)
1. Fst between selected (N=2 pops or N=1) and unselected using DroneCalculateFSTfromVCF.sh (copy and paste)

2. Run  MovingAveFSTfromPLINK.r <PLINK Fst> <WInsize> to create moving averages of FST along chromosomes (converted with updated perl script within R) saves an .RData file for next steps

3. Get out high FST windows with HighWindowsFST.r <RData> <cutoff> (saved as dataframe called HighSNPs)

4. If 2 independent FST analyses where used, get overlapping windows.

e.g.

<pre><code>Rscript MovingAveFSTfromPLINK.r 
	/media/data1/forty3/drone/vcf_drone/pop2_vs_sel.weir.fst 10000</code></pre>
<pre><code>Rscript HighWindowsFST.r \
	/media/data1/forty3/drone/vcf_drone/pop2_vs_sel.RDATA 0.95</code></pre>

save.image(file="/media/data1/forty3/drone/FST/SelvsCon/SelectedvsControlFSTPvals.RData")

5. I also ran [pFST](https://github.com/jewmanchue/vcflib/wiki/Association-testing-with-GPAT) using pFST.sh and added this to the previous dataframe of FST values. Then ran FDR using Storey's Q method. Significant FST was based on q<0.01 (or, -log10(P)>2.61)


##Admixture Analyses
I used the final VCF file, paired with SNPs called from Harpur et al. 2014 to run ADMIXTURE. 


1. Get out common SNPs between AMC and AHB in a merged VCF file, create a .ped file. For this, I copied AMC.ped and AMC.map from my AFZ project. This contains SNPS for AMC with MAF 0.05. 

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


2. Use the .ped files to create .bim, .fam,  and .bed for ADMIXTURE. I did this for all significant FST SNPs and for a random selection of SNPs. Then, run ADMIXTURE for K=3.

<pre><code>
plink --file AMCHYGHIGH --noweb --make-bed --out AMCHYGHIGH
/home/brock/admixture/admixture  --cv=10 AMCHYGHIGH.bed 3 -j2 | tee log3.out

plink --file AMCHYGHIGH --noweb --make-bed --out AMCHYGHIGH
/home/brock/admixture/admixture  --cv=10 AMCHYGHIGH.bed 3 -j2 | tee log3.out
</code></pre>

<!---
#saved in HygieneHighFSTADMIXTURE.xlsx
-->


###Plotting Data
All FST plot scripts can be found as .r files

1. FST Histogram, FST by chromosome (and stupid legend) = FSTHygienPlot.r

2. Boxplot of admixture in FST SNPs = HighFSTAdmixPlot.r

3. I ran REVIGO using the web app (no R version :( ). Uploaded REVIGO_AHighFST.r and REVIGO_AllNSYNHighFST.r for REVIGO plots of all sig FST and all sig NSYN Fst. 

N. I may re-run Rcircos, see Rcircos_Drone.r















