####################################################
# I want to try a SET TEST
#######
#file="ControlHighFST.map"
#http://pngu.mgh.harvard.edu/~purcell/plink/anal.shtml#set

plink --file DroneSelection --extract CandidateSNPs.snp --out ControlHighFST --noweb --make-bed --keep controlBees.txt 
plink --bfile ControlHighFST  --recode --out ControlHighFST --noweb

plink --file ControlHighFST --set-test --set highFST.set --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --allow-no-sex --assoc --mperm 1000 --out CONTROLSET --set-p 0.01 --set-max 10 --qt-means


plink --file DroneSelection --extract CandidateSNPs.snp --out SelHighFST --noweb --make-bed --remove controlBees.txt 
plink --bfile SelHighFST  --recode --out SelHighFST --noweb

plink --file SelHighFST --set-test --set highFST.set --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --allow-no-sex --assoc --mperm 1000 --out SelSET --set-p 0.01 --set-max 10 --qt-means




##OK

R
load(file="/media/data1/forty3/drone/FST/SelvsCon/SelectedvsControlFST.RData")


map=read.table(file="ControlHighFST.map",header=F)
mns=read.table("CONTROLSET.qassoc.means",header=T)


CONset=read.table("CONTROLSET.qassoc.set.mperm",header=T)
CONset=unlist(strsplit(as.vector((CONset$SNPS)),"[|]"));CONset=CONset[which(CONset!="NA")]
CONpos=map[map$V2 %in% CONset,]
chroms=unique(CONpos$V1)
# 3  5  6  7  9 11 12 14 15 16

SELset=read.table("SelSET.qassoc.set.mperm",header=T)
SELset=unlist(strsplit(as.vector((SELset$SNPS)),"[|]"));SELset=SELset[which(SELset!="NA")]
SELpos=map[map$V2 %in% SELset,]