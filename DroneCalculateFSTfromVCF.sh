




#####
#Fst Analyses with VCFTools
	#Calculate FST  all selected VS control (2 vs 1+3=Selpop.txt)
vcftools --maf 0.05 --vcf DroneSelectionFinal.recode.vcf --weir-fst-pop pop2.txt --weir-fst-pop Selpop.txt --out pop2_vs_sel 



#and each Separately:
vcftools --maf 0.05 --vcf DroneSelectionFinal.recode.vcf --weir-fst-pop pop2.txt --weir-fst-pop pop1.txt --out pop2_vs_pop1 
vcftools --maf 0.05 --vcf DroneSelectionFinal.recode.vcf --weir-fst-pop pop2.txt --weir-fst-pop pop3.txt --out pop2_vs_pop3