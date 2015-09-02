
<!---
There be dragons here:


##File Creation
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


4a. output LD in  populations (Haven't done here)
<pre><code>
plink --noweb --file AllSamplereINDEP 
	--ld-window-kb 1000  
	--remove controlBees.txt 
	--out DroneSelLD
</code></pre> 

<pre><code>
plink --noweb --file AllSamplereINDEP 
	--ld-window-kb 1000  
	--keep controlBees.txt  
	--out DroneCONLD
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





	
##If UnPhased Data wanted:	
6. Frequency removal
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
	do plink --noweb --file DroneSamps_$K --maf 0.05  --recode --out DroneSamps_$K.UNphasedMAF ; done








##Association Analysis
1.  Cochran-Mantel-Haenszel (CMH) tests with unphased data
<pre><code>
plink --noweb --file  DroneSamps_6.UNphasedMAF --pheno /media/data1/forty3/drone/vcf_drone/CMHpheno.txt --mh2 --within SELCONcluster.txt --out CLUSTEREDUNPHMAF2_6 
</code></pre> 

2. Haplotype Association with disease staus
<pre><code>
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; 
	do plink --noweb 
		--file  DroneSamps_$K.UNphasedMAF 
		--pheno /media/data1/forty3/drone/vcf_drone/CMHpheno.txt 
		--hap-window 3  
		--hap-assoc 
		--out Dhap$K ; done
</code></pre> 
	
3. Haplotype Association with quantitiative and unphased
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
	do plink --noweb --file DroneSamps_$K.phasedMAF --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --hap-window 3 --hap-assoc --out hap_$K ; done
	

	
	
4. Quantitative permutation with phased
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
	do plink --noweb --file DroneSamps_$K.phasedMAF --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --assoc --perm --out PERMhap_$K ; done

	
4. Quantitative permutation with unphased
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
	do (plink --noweb --file DroneSamps_$K.UNphasedMAF --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --assoc --perm --out PERMUNhap_$K ) & 
	done
	
	


5. Rpart across all
I ran testSFS.r on each chromosome	and outputs to RPARTAssocSNPs




16.2:431139 comes up every time,....
plink --noweb --file DroneSamps_16.phasedMAF --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --proxy-assoc 16.2:431139 
http://pngu.mgh.harvard.edu/~purcell/plink/proxy.shtml

	
10. See Analysis.R and Rcircos_Drone.r



-->




<!---
##########
9. Write out phased, tped files for association
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
	do plink --noweb --file DroneSamps_$K --pheno /media/data1/forty3/drone/vcf_drone/CMHpheno.txt --mh --within SELCONcluster.txt --out CLUSTEREDnoShape_$K ; done








	
#>75 and with clusters
plink --noweb --tfile DroneSamps_6.phased.sample --pheno /media/data1/forty3/drone/vcf_drone/CMHpheno.txt --mh 
	
plink --noweb --tfile DroneSamps_6.phased.sample --pheno /media/data1/forty3/drone/vcf_drone/CMHpheno.txt --mh --within SELCONcluster.txt --out CLUSTERED
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	Rscript HapstoTPED.r /media/data1/forty3/drone/vcf_drone/DroneSamps_16.phased.haps.out; done


plink --noweb --tfile DroneSamps_6.phased.sample --homozyg	


plink --noweb --tfile DroneSamps_6.phased.sample --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --remove controlBees.txt --hap-window 3 --hap-assoc --adjust --out DroneSampsOutSELPHASE_6	
	
plink --noweb --tfile DroneSamps_6.phased.sample  --homozyg --out HOMODef --homozyg-snp 50 --homozyg-kb 500
	

plink --noweb --tfile DroneSamps_16.phased.sample --recode --make-bed
plink --noweb bfile plink --assoc 


plink --noweb --tfile DroneSamps_6.phased.sample --pheno /media/data1/forty3/drone/vcf_drone/CMHpheno.txt --mh


#>75 and with clusters
plink --noweb --tfile DroneSamps_6.phased.sample --pheno /media/data1/forty3/drone/vcf_drone/CMHpheno.txt --mh 
	
plink --noweb --tfile DroneSamps_6.phased.sample --pheno /media/data1/forty3/drone/vcf_drone/CMHpheno.txt --mh --within SELCONcluster.txt --out CLUSTERED

	
plink --noweb --bfile plink --hap-window 3  --hap-assoc --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --out DroneSampsOutPHASE_16	
	
	

9. Haplotype associations
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
	do plink --noweb --tfile DroneSamps_$K.phased.sample --hap-window 3 --hap-assoc --adjust --out DroneSampsOutPHASE_$K; done


plink --noweb --tfile DroneSamps_16.phased.sample --hap-window 3 --hap-assoc --adjust --out DroneSampsOutPHASE_16	

plink --noweb --tfile DroneSamps_6.phased.sample --recode --make-bed
	
	
	
plink --file AllSamplereINDEP  --hap-window 3 --hap-assoc --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --adjust --noweb



for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
	do plink --noweb --file DroneSamps_$K  --homozyg --out DroneHET_$K; done



plink --file AllSamplereINDEP  --homozyg --noweb



for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
	do plink --noweb --file DroneSamps_$K --keep controlBees.txt --hap-window 3 --hap-assoc --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --adjust --out DroneSampsOutCON_$K; done

















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