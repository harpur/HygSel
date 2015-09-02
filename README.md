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


###Plotting Data
All FST plot scripts can be found as .r files

1. FST Histogram, FST by chromosome (and stupid legend) = FSTHygienPlot.r

2. Boxplot of admixture in FST SNPs = HighFSTAdmixPlot.r

3. I ran REVIGO using the web app (no R version :( ). Uploaded REVIGO_AHighFST.r and REVIGO_AllNSYNHighFST.r for REVIGO plots of all sig FST and all sig NSYN Fst. 
N. I may re-run Rcircos, see Rcircos_Drone.r



