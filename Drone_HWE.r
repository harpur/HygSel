###
# HWE Filter sites
###











args <- commandArgs(trailingOnly = TRUE)
fil = as.character(args[1]) #file name 

#Output site HWE with vcftools --------
system(paste("vcftools --vcf ",fil,"--hardy", sep=" "))


#Load and mung, get significant sites ---------------------
hwe=read.table(file="out.hwe",header=T,colClasses=c("character","numeric","character","character","numeric", "numeric"))
hwe = hwe[which(hwe$P<0.05),]

hwe$obs = gsub("^[0-9]*/", "", hwe$OBS.HOM1.HET.HOM2.)
hwe$obs = as.numeric(gsub("/[0-9]*$", "", hwe$obs))
hwe$ex = gsub("/[0-9]*.[0-9]*$", "", hwe$E.HOM1.HET.HOM2.)
hwe$ex = as.numeric(gsub("^[0-9]*.[0-9]*/", "", hwe$ex))
hwe$OE = hwe$obs - hwe$ex
 table(hwe$CHR[hwe$OE>0])
