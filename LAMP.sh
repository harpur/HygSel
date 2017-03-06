###
# LAMP and Admixture with high vs low regions
###


#Runs LAMP against hygiene data set vs AMC pops
	#cd /media/data1/forty3/drone/vcf_drone/lamp
	#AllBees_SNPs.raw.vcf  - AMC SNPs
	#DroneSelectionFinal.recode.vcf - Drone SNPs
	
	
	
	
	
#extract SNPs ------------------------------	
vcftools --vcf AllBees_SNPs.raw.vcf  --max-alleles 2 --recode --out ALLMAF
vcftools --vcf DroneSelectionFinal.recode.vcf --maf 0.05 --max-alleles 2 --max-missing-count 5 --recode --out droneMAF
cut -f 1,2  droneMAF.recode.vcf > SNPs
sed -i '/#/d' SNPs

#Calculate FST, extract --------------------
vcftools --vcf ALLMAF.recode.vcf --max-alleles 2 --positions SNPs  --weir-fst-pop c.txt --weir-fst-pop s.txt --out CA
vcftools --vcf ALLMAF.recode.vcf  --max-alleles 2 --positions SNPs --weir-fst-pop c.txt --weir-fst-pop m.txt --out CM
vcftools --vcf ALLMAF.recode.vcf --max-alleles 2 --positions SNPs  --weir-fst-pop m.txt --weir-fst-pop s.txt --out MA
Rscript extractFST.r

#extract sites from high FST for all bees, create MAP and PED -----------------------
vcftools --vcf ALLMAF.recode.vcf --positions aims --recode --out chrAIMS


bgzip chrAIMS.recode.vcf
bgzip droneMAF.recode.vcf
tabix -p vcf chrAIMS.recode.vcf.gz
tabix -p vcf droneMAF.recode.vcf.gz
vcf-merge chrAIMS.recode.vcf.gz droneMAF.recode.vcf.gz > lamp.samp.vcf


vcftools --vcf lamp.samp.vcf --maf 0.05 --max-alleles 2 --thin 500 --max-missing-count 5 --recode --out  lamp.samp
vcftools --vcf lamp.samp.recode.vcf --plink --keep c.txt  --out C
vcftools --vcf lamp.samp.recode.vcf --plink  --keep s.txt  --out A
vcftools --vcf lamp.samp.recode.vcf --plink  --keep m.txt  --out M 
cat C.ped M.ped A.ped > AMC.ped 
mv C.map AMC.map





#Convert scaff to chrom --------------------------
Rscript /media/data1/forty3/drone/git/ScaffMaptoChr.r AMC.map
mv AMC.map AMCSCAFF.map
mv NA.chrom.map AMC.map
plink --noweb --file AMC --recode  --out AMCr


#Output 1 file per Chromosome ---------------------
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
do plink --noweb --file AMCr --recode --out AMC_$K --chr $K ; done


#plink --noweb --file AMCr --recode --out AMC_1 --chr 1 




#Phase all chrs -----------------------------------
	#cd /media/data1/afz/shapeit/
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
do ./shapeit --rho 0.39 -P /media/data1/forty3/drone/vcf_drone/lamp/AMC_$K -T 2 --window 0.5 --output-max /media/data1/forty3/drone/vcf_drone/lamp/AMC_$K ; done

#./shapeit --rho 0.39 -P /media/data1/forty3/drone/vcf_drone/lamp/AMC_1 -T 2 --window 0.5 --output-max /media/data1/forty3/drone/vcf_drone/lamp/AMC_1



cd /media/data1/forty3/drone/vcf_drone/lamp
#create LAMP input files for ancestors -----------------------------------	
	#ancestral haplotype files!
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
do 
	cut -d ' ' -f6-23 AMC_$K.haps > C_$K
	cut -d ' ' -f24-41 AMC_$K.haps > M_$K
	cut -d ' ' -f42-63 AMC_$K.haps > A_$K; done

for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
do 
	bash turnfile.sh C_$K
	bash turnfile.sh M_$K
	bash turnfile.sh A_$K
	sed -i 's/ //g' C_${K}_r 
	sed -i 's/ //g' M_${K}_r 
	sed -i 's/ //g' A_${K}_r ; done



		
#pos files!
	#for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
	#do 
	#	cut -f 4 AMC_$K.map > pos_$K; done
	#for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
	#do 
	#	grep ^$K[.] aims > pos_$K; done
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
do 
cut -d' ' -f3 AMC_${K}.haps > pos_${K} ; done
#sed -i 's/:/\t/g' pos_${K}; done
	
 
	#cut -d' ' -f 2 AMC_1.haps > pos_1
	#sed -i 's/:/\t/g' pos_1





#Create Sample Inputs ------------
	#for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
	#do 
	#	cut -f 2 AMC_${K}.map > focalSNPs_$K
	#	sed -i 's/:/\t/g' focalSNPs_$K; done


for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
do 
	vcftools --vcf lamp.samp.recode.vcf --recode --keep drones.txt --positions pos_${K}  --out samps_${K}
	Rscript VCFtoLAMP.r samps_${K}.recode.vcf; done
	
	
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
do 	
	cut -d ' ' -f2- genos.samps_${K}.recode.vcf > samp${K}.gen
	sed -i 's/ //g' samp${K}.gen; done
	

		
mv samp[0-9]*.gen /media/data1/afz/lamp/LAMPLD-v1.1
mv *r /media/data1/afz/lamp/LAMPLD-v1.1
mv pos* /media/data1/afz/lamp/LAMPLD-v1.1



cd /media/data1/afz/lamp/LAMPLD-v1.1



for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
do 
	perl run_LAMPLD.pl pos_$K C_${K}_r M_${K}_r A_${K}_r samp${K}.gen chr${K} #0,1,2
	perl convertLAMPLDout_BAH.pl chr${K} chr${K}.out; done	

	
#summarize results ------------------
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
do 
	Rscript LAMPanalysis.r $K; done	
	
	
	
	
	
###########################################################################################################	
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
do 
	Rscript VCFtoLAMP.r samps_${K}.recode.vcf; done
	
	
	
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
do 
	cut -d ' ' -f2- genos.samps_${K}.recode.vcf > samp${K}.gen
	sed -i 's/ //g' samp${K}.gen; done	
	


	
#Run LAMP---------------------------------	
	
mv samp[0-9]*.gen /media/data1/afz/lamp/LAMPLD-v1.1
mv *r /media/data1/afz/lamp/LAMPLD-v1.1
mv pos* /media/data1/afz/lamp/LAMPLD-v1.1

cd /media/data1/afz/lamp/LAMPLD-v1.1



for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
do 
	perl run_LAMPLD.pl pos_$K C_${K}_r M_${K}_r A_${K}_r samp${K}.gen chr${K} #0,1,2
	perl convertLAMPLDout_BAH.pl chr${K} chr${K}.out; done	



#summarize results ------------------
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; \
do 
	Rscript LAMPanalysis.r $K; done	

	
	

	
	
	
	
	for K in 1 ; \
do 
	perl run_LAMPLD.pl pos_$K C_${K}_r M_${K}_r A_${K}_r samp${K}.gen chr${K} #0,1,2
	perl convertLAMPLDout_BAH.pl chr${K} chr${K}.out; done	
	
	
	