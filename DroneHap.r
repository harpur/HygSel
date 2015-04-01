 ###
 #Drone Selection - Haploblock identification (HaploView)

 
#I want to find haploblocks in each population to identify which regions where selected. 
	#Selected versus unselected haploblocks

cd /media/data1/forty3/drone/HaploView


java -jar Haploview.jar -n -pedfile DroneSelection.ped -map DroneSelection.map -chromosome 6 -blockoutput ALL -memory 10000 -skipcheck -startpos 0.01 -endpos 5
	

	
	
	
	
	
	
	
	
#Alternatively, LD blocks in VCFtools:
vcftools --maf 0.05 --ld-window-bp 5000 --vcf DroneSelectionFinal.recode.vcf --keep Selpop.txt --out Selected 
vcftools --maf 0.05 --ld-window-bp 5000 --vcf DroneSelectionFinal.recode.vcf --remove Selpop.txt --out Control 


#Or in PLINK:
plink --file DroneSelection --ld-window-kb 15 --noweb --allow-no-sex --keep controlBees.txt --out Control
plink --file DroneSelection --ld-window-kb 15 --noweb --allow-no-sex --remove controlBees.txt --out Selected
plink --file DroneSelection --ld-window-kb 15 --noweb --allow-no-sex --remove FASBees.txt --out FASSelected
	
	
	
R
sel = read.table(file = "Selected.ld",header = T)
con = read.table(file = "Control.ld",header = T)

#for CHR6, 1,1700000
con = con[con$CHR_A =="6",]
con = con[con$BP_A <= 1700000,]
sel = sel[sel$CHR_A = ="6",]
sel = sel[sel$BP_A <= 1700000,]

dis = abs(sel$BP_A - sel$BP_B)#/1000
brks = hist(dis, breaks=50, plot=F)$breaks
grp = cut(dis, breaks=brks)
r2meansSEL = tapply(sel$R2, grp, FUN = mean)
plot(r2meansSEL,
	ylab = "Average R2",
	xlab = "Binned Distance",
	pch = 19)
dis = abs(con$BP_A - con$BP_B)#/1000
brks = hist(dis, breaks=50, plot=F)$breaks
grp = cut(dis, breaks=brks)
r2meansCON = tapply(con$R2, grp, FUN = mean)
points(r2means, col="blue", pch=19)




#	To get the the LD decay plot generate a textfile with two columns - the first is the distance between the two SNPs (i.e. BP_B - BP_A, so if the first snp is at position 1000 on the chromosome, and the second is at position 2150, the distance is 1150 ) and the second is the corresponding R2 value (yoru last column). So you get a file with R2 values for SNPs certain distances apart. distance R2 5 0.2 5 0.3 67 0.2 67 0.4 67 0.5

#Then for each distance apart, calculate an average R2 (you need to generate a script to do this e.g. in python/perl) i.e. distance averageR2 5 0.25 67 0.3667

#Then plot that file to get the LD decay.

#Depending on how many points you have, you may want to using a sliding window average script for the plot. In the example, I think this is what they did - in terms of 1kb intervals, i.e. they did a sliding window - with a window size of 1kb and step size of 1kb, and calculated average r2.
	
	
	
	
	
	
	
	
	
	
	
	
	
#File format for Haploblock	
	
source("/media/data1/forty3/brock/scripts/VCFFunctions.r")	
i = "/media/data1/forty3/drone/vcf_drone/DroneSelectionFinal.recode.vcf"
x = Read.VCF()