###
# Running hapFLK across hygiene samples
###

# HAPFLK Tutorial:
#These three files are needed to run hapflk:
	#hapmap3-lct.ped : the ped file (plink format) with a particularity. First column of the file (FID) contains the population identifier of individuals (rather than a family ID).
	#hapmap3-lct.map : the map file
	#kinship.txt: the kinship matrix estimated from the Reynolds distances

	
#I'm running this using control bees against selected pops	
	

#create data set -------	
vcftools --vcf DroneSelectionFinal.recode.vcf --plink --out outSNP
Rscript ScaffMaptoChr.r outSNP.map
mv outSNP.map outSNP_scaff.map
mv NA.chrom.map outSNP.map
plink --noweb --file outSNP --recode  --out AllSamplere  



#output PED files for each chromosome -----
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
do plink --noweb --file AllSamplere  --recode --out DroneSamps_$K --chr $K ; done
	#plink --noweb --file AllSamplere  --recode --out DroneSamps_5 --chr 5

#Phase---------------
#here, I used rho 0.39 (from Wallberg) across the whole genome
for K in  16 2 3 4 5 6 7 8 9 10 11 12 13 14 15 1; \
	do ./shapeit --rho 0.39 -P DroneSamps_$K -T 2 -O DroneSamps_$K.phased ; done
	
for K in  16 2 3 4 5 6 7 8 9 10 11 12 13 14 15 1; \
do
./shapeit -convert \
    --input-haps DroneSamps_$K.phased \
    --output-vcf gwas_$K.phased.vcf
vcftools --vcf gwas_$K.phased.vcf --plink --out gwas_$K --max-alleles 2; done



#define pops-------------------
	#cp /media/data1/forty3/drone/vcf_drone/gwas_[0-9]* /media/data1/forty3/drone/PLINK/hapflk-1.3.0/
for K in  16 2 3 4 5 6 7 8 9 10 11 12 13 14 15 1; \
do
	awk '!($1="")' gwas_$K.ped > gwas2.ped
	paste fist2 gwas2.ped > gwas_$K.ped; done


#Run HAPFLK ------------------	
for K in  16 2 3 4 5 6 7 8 9 10 11 12 13 14 15 1; \
do
	hapflk --file gwas_$K --keep-outgroup --outgroup C --phased -K 20 --nfit=20 --ncpu=17
	mv hapflk.hapflk hapflk_$K.hapflk; done


#summarize hapflk results------------
Rscript <HERE?> 16 /media/data1/forty3/drone/PLINK/hapflk-1.3.0/chr/


	
	













