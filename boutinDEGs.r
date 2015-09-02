###
#Boutin DEGs overlap with my analyses
###



# Boutin's Paper  -------------------
	# http://www.biomedcentral.com/1471-2164/16/500/
	# I pulled the DEG list and saved it as BoutinDEGs.xlsx

	
# Load datasets  -------------------	

#Fst Data
load(file = "/media/data1/forty3/drone/vcf_drone/pop2_vs_sel.RDATA")

#GFF files
gff = read.table(file = "/media/data1/forty3/brock/scripts/GFF/NCBIGFF.txt")
gff=gff[gff$V5=="gene",] #only need genes

#DEG list
	#this is only protein-coding genes in placed scaffolds of the genome
degs = read.table(file="/media/data1/forty3/drone/git/boutinDEGs.txt",header=T)
gff.DEG=gff[gff$V2 %in% degs$GB,] #gff for only these degs


# Are these genes within high FST windows?  -------------------	


DEG.Fst = c()
chrom = intersect(gff.DEG$V1, HighWindows$chr)
for(i in chrom){
	win.temp = HighWindows[HighWindows$chr==i,]
	deg.temp = gff.DEG[gff.DEG$V1==i,]
	blah=outer(as.numeric(deg.temp$V3), as.numeric(as.character(win.temp$Rwin))+50000, "<=") 
	blah1=outer(as.numeric(deg.temp$V4), as.numeric(as.character(win.temp$Fwin))-50000, ">=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	DEG.Fst = rbind(temp, DEG.Fst)
	print(i)
}






