###
# Plot a gene region
###



#example here is RAD17, I guess





#Arguments --------------------------
GB = "GB45774"

#Load functions and packages -----------------------------
source("/media/data1/forty3/drone/git/GenomeR/VCFFunctions.r")
sorter <- function(x){
	srts = sort(summary(as.factor(x)), decreasing=T)
	if(length(srts)<2){
		return(names(sort(summary(as.factor(x)), decreasing=T)[1]))
	}else{
		if(sum(srts[1] == srts[2])){
		return("NA")
		}else{
			return(names(sort(summary(as.factor(x)), decreasing=T)[1]))
		}
}
}


#Load Dataframes ----------------------------------------------
Functional.Sites = read.table(file="HYG.snpeff.eff", header=F,sep="\t", quote="")
Gene.Regions = read.table(file="/home/amel45/AM45/GFF3_2.bed")
samps = read.table(file="/media/data1/forty3/drone/git/data/DroneSamps.txt")

#Extract Gene Region ---------------------
	#this is, for now, hard-coded
Gene.Regions = Gene.Regions[Gene.Regions$V4==GB,]
Functional.Sites = Functional.Sites[Functional.Sites$V10==GB,]
#table(as.character(Functional.Sites$V16))
#Functional.Sites = Functional.Sites[-grep("INTRAGENIC:*", Functional.Sites$V16),]
st = min(Gene.Regions$V2)
en = max(Gene.Regions$V3)
chr = as.character(Gene.Regions$V1[1])
system(paste("vcftools --vcf out.indel.dp.q.miss.recode.vcf --from-bp", st, "--to-bp", en, "--chr", chr, "--recode --out HYGHIGHHAPS",sep=" ") ) 
i = "HYGHIGHHAPS.recode.vcf" 
vcf = Read.VCF()


#Extract Alleles for each Haploid Drone ---------------------
	#this is, for now, hard-coded
pos = sapply(vcf,function(x) return((x[2])))
pos = pos[-1]
chr = sapply(vcf,function(x) return((x[1])))[-1]
#snp = paste(chr, pos, sep="_")
allele  = sapply(vcf,function(x) substr(unlist(x)[10:134],1 ,1 )) 
allele = allele[,-1]
snps = data.frame(cbind(chr, pos))



#Prepare Plot Parameters ------------------
ysize = nrow(allele)+25
mid = (ysize+(ysize-20))/2
allele[allele=="0"] = "grey"
allele[allele!="grey"] = "blue"


allele[pos %in% Functional.Sites$V2] = 




#Plot -------------------------------------
plot(NULL, yaxt = "n", xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
	xlim = c(st-100, en+100), 
	ylim = c(0, ysize), 
	xlab = NA, 
	ylab = NA, 
	las = 1, 
	pch = 19)
segments(st-100,mid,en+100,mid) 	#plot intron line
rect(Gene.Regions$V2, ysize-20, Gene.Regions$V3, ysize, col="grey")	#plot exons

#png(file="MlineageAlleles.png")
points(pos,y=rep(1,length(pos)), col=allele[1,],pch=15,cex=0.3)
for(i in 2:125){
	points(pos,y=rep(i,length(pos)), col=allele[i,],pch=15,cex=0.3)
}

abline(h = 88, col="red", lty=2)







    ysize <- ncol(hap1.N)+gap
    par(bg=bg, fg=fg, col.lab=fg, col.axis=fg, cex=0.3, cex.sub=0.4, cex.main=0.4, cex.axis=0.3, cex.lab=0.4, lwd=0.5)
    par(mfrow=c(1,1), oma = c(3,3,1,1), mar = c(0.5,0.5,0.5,0.5))

    
    # Output window to console
	
	
	




























