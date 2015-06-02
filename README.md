# Workflow and Supplementary Material for Drone Selection Experiment









##VCF Creation:
- VCFCreation_DroneSelection.txt
- DroneVCFcreation.txt




##FST Analyses
###After Alignment, SNP calling, and trimming procedures (VCF Creation above)
1. Fst between selected (n=2 pops or N=1) and unselected using DroneCalculateFSTfromVCF.sh (copy and paste)
2. Run  MovingAveFSTfromPLINK.r <PLINK Fst> <WInsize> to create moving averages of FST along chromosomes (converted with updated perl script within R) saves an .RData file for next steps
3. Get out high FST windows with HighWindowsFST.r <RData> <cutoff> (saved as dataframe called HighSNPs)
4. If 2 independent FST analyses where used, get overlapping windows.

e.g.

<pre><code>Rscript MovingAveFSTfromPLINK.r 
	/media/data1/forty3/drone/vcf_drone/pop2_vs_sel.weir.fst 10000</code></pre>
<pre><code>Rscript HighWindowsFST.r 
	/media/data1/forty3/drone/vcf_drone/pop2_vs_sel.RDATA 0.95</code></pre>




##Phase
1. Create PED and MAP files
<pre><code>cd /media/data1/forty3/drone/vcf_drone</code></pre>
<pre><code>vcftools --vcf DroneSelectionFinal.recode.vcf --plink --out AllSample</code></pre>

2. Convert scaff to chrom
<pre><code>Rscript /media/data1/afz/git/ScaffMaptoChr.r AllSample.map</code></pre>

3. Re-order
<pre><code>plink --noweb --file AllSample --recode  --out AllSamplere</code></pre> 


4. Output 1 file per Chromosome
<pre><code>
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
do plink --noweb --file AllSamplere  --recode --out DroneSamps_$K --chr $K ; done
</code></pre> 









<!---

		5. Run Shapeit
		for K in  16 2 3 4 5 6 7 8 9 10 11 12 13 14 15 1; \
		do ./shapeit -P /media/data1/forty3/drone/vcf_drone/DroneSamps_$K -T 2 -O/media/data1/forty3/drone/vcf_drone/DroneSamps_$K.phased ; done




						haps=read.table(file="DroneSamps_6.phased.haps ")


						6. Add in Allelic information into 
						for K in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 ; \
						do python  /home/sani/task1/Phased_Script_PP.py DroneSamps_16.phased.haps ; done


						plink --noweb --file DroneSamps_16.phased.haps --recode  --out test --chr $K


						R
						source("/media/data1/forty3/brock/scripts/GOgetter.r")
						source("/media/data1/forty3/brock/scripts/VarFunct.r")
						source("/media/data1/forty3/brock/scripts/movingavg.r")


						samp=read.table(file="DroneSamps_6.phased.sample",header=T)
						samp=samp[-1,];samp=samp[-3]
						samp=rbind(samp,samp);samp=samp[order(samp[,1]),]
						haps=read.table(file="DroneSamps_6.phased.haps",header=F)
						haps$zer=rep("0", nrow(haps))
						maps=haps[c(1,2,ncol(haps), 3)]
						write.list(maps, file="DroneSampsPH_6.map")
						haps=haps[-c(1,2,ncol(haps), 3,4,5)]
						haps1=as.matrix(haps)
						haps1[haps1=="0"]=23;haps1[haps1=="1"]=2;haps1[haps1=="23"]=1;
						
			
						thaps=cbind(maps, haps1)
						write.list(thaps, file="DroneSampsPH_6.tped")
						write.list(samp, file="DroneSampsPH_6.tfam")




#should trim these of R<0.2
plink --noweb --tfile DroneSampsPH_6  --indep 50 5 2 
plink --noweb --tfile DroneSampsPH_6  --extract plink.prune.in --make-bed --out DroneSampsPH_6P
plink --bfile DroneSampsPH_6P  --hap-window 3 --hap-assoc --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb

#should consider changing the model

plink --tfile DroneSampsPH_16  --assoc --qt-means --adjust --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb

plink --file DroneSamps_6  --assoc --qt-means --adjust --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --out unphase


#testing SKAT
plink --tfile DroneSamps_16   --make-bed --out CHR16 --noweb		
plink --file DroneSamps_16   --recode12 --out CHR16 --noweb	

R
require("SKAT")
geno=read.table(file="CHR16.ped",header=F)






./famhap19 CHR16FAMMAP name2




					
						
				
## Associations
1. 




### Old, remove later
1. DronePLINKset.sh identifies high FST across the genome between selected and control populations, then 
2. Pass FST to DroneFST.r and identify high FST regions within the genome- outputs all SNPs within those regions (CandidateSNPs.snp)
3. DronePLINKset.sh uses a modified CandidateSNPs.snp (Candidates98_ALL.set ) to run the PLINK set test within each high FST region
4. DronePLINKset.sh uses  CandidateSNPs.snp to dientify independant SNPs within each region (called CandidateSET.xxx)
5. CandidateSET.xxx is then used for Recursive Partitioning and Regression Tree (RpartAssociations.r) 

















-->