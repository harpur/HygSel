#!/usr/bin/Rscript
###
# Identify High FST windows within .RData


#Requires RData file with FST created from MovingAveFSTfromPLINK.r 














###
#Arguments
args = commandArgs(trailingOnly = TRUE)
RData = args[1]
cutoff = args[2] #0.95%

###
#Functions
source("/media/data1/forty3/brock/scripts/GOgetter.r")
source("/media/data1/forty3/brock/scripts/VarFunct.r")
source("/media/data1/forty3/brock/scripts/movingavg.r")
#####


###
#Script Begins
load(file=RData)
highFstWindow=quantile(unlist(Fst$x),as.numeric(cutoff))
HighFstPositions=as.matrix(apply(Fst, 1, function(x) unlist(as.numeric(unlist(x[3]))[which((as.numeric(unlist(x[2]))>=highFstWindow))])))
#SNPs within high FST areas
	#NOTE: Row names in HighFstPositions are INCORRECT. shoudl be 1,2,3,etc.
row.names(HighFstPositions)=c(1:16)
matr=c()
for(i in 1:16){
	#i=3
	test=as.vector(unlist(HighFstPositions[i,1]))#Two regions in CHR3
	if(length(test)==0){print(row.names(HighFstPositions)[i])}else{
		test2=c(test,test[length(test)]);test2=test2[-1]
		x=test2-test
		highs =  which(x>5000)
		Fwin = c(test[1],test2[highs])
		Rwin = c(test[highs], +test[x==0])
		#abline(v=c(Rwin),col="red")
		#abline(v=c(Fwin),col="blue")
		matr=rbind(matr,cbind(Fwin,Rwin,"chr"=rep(row.names(HighFstPositions)[i], length(Fwin))))		
}
}



#Matr now contains the start and end of all high FST windows at >=cutoff% data
HighWindows = data.frame(matr)
HighSNPs=c()
chr=(unique(HighWindows$chr))
for (i in chr){
	test=fst23[fst23$CHROM==i,]
	low=as.numeric(as.character(unlist(HighWindows$Fwin[which(HighWindows$chr==i)])))
	high=as.numeric(as.character(unlist(HighWindows$Rwin[which(HighWindows$chr==i)])))
	blah=outer(as.numeric(unlist(low)), as.numeric(unlist(test$POS)), "<=") 
	blah1=outer(as.numeric(unlist(high)), as.numeric(unlist(test$POS)), ">=") 
	blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T)) #The gene region will be the colum variable
	blah=(test[blah[,2],])
	blah=blah[!duplicated(blah),]
	HighSNPs=rbind(HighSNPs,blah)
	print(i)
}

HighSNPs$test = paste(HighSNPs$CHROM, HighSNPs$POS, sep=":")
fst.raw$V2=gsub("chr","",fst.raw$V2)
fst.raw$test = paste(fst.raw$V2, fst.raw$V3, sep=":")

test=merge(HighSNPs, fst.raw, by="test")
test$CHROM.y=gsub("Group","",test$CHROM.y)
HighSNPs = test[c(2,3,4,11)] ;  names(HighSNPs)=c("CHROM","POS","FST","SNP")


save.image(file=RData)






