###
#Calculate Fst by Gene between populations
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

# FST for each Gene  -------------------	
Genic.Fst = c()
chrom = intersect(gff$V1,fst23$CHROM)
for(i in chrom){
	win.temp = gff[gff$V1==i,]
	deg.temp = fst23[fst23$CHROM==i,]
	blah=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$V3)), ">=") 
	blah1=outer(as.numeric(deg.temp$POS), as.numeric(as.character(win.temp$V4)), "<=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	temp = deg.temp[blah[,1],]
	temp = cbind(temp, win.temp[blah[,2],])
	Genic.Fst = rbind(temp,Genic.Fst)
	print(i)
}
Genic.Fst.mean = with(Genic.Fst,aggregate(as.numeric(FST), by=list(V2),mean))
Genic.Fst.len = with(Genic.Fst,aggregate(as.numeric(FST), by=list(V2),length))
Genic.Fst.mean$len = Genic.Fst.len$x 
Genic.Fst.len = NULL



# Are DEGs higher FST?  -------------------	

#DEG list
	#this is only protein-coding genes in placed scaffolds of the genome
degs = read.table(file="/media/data1/forty3/drone/git/boutinDEGs.txt",header=T)
gff.DEG=gff[gff$V2 %in% degs$GB,] #gff for only these degs

`
test = Genic.Fst.mean[Genic.Fst.mean$len>10,]

boxplot(test$x[!(test$Group.1 %in% gff.DEG$V2)], test$x[test$Group.1 %in% gff.DEG$V2],notch=T)







