




# Functions -------------------	
source(file="/media/data1/forty3/brock/scripts/VarFunct.r")

# Rpart Data -------------------	
#All associated SNPs are saved in RPARTAssocSNPs -> CHR 1 couldn't run, too many SNPs

# CMH Data  -------------------	
#All associated SNPs are saved in CLUSTEREDMAF_*.cmh and were unphased

cmh.assoc = c()
fils = list.files(pattern="^CLUSTEREDUNPHMAF_")
for(fil in fils){
	#fil = fils[13]
	assoc = read.table(file=fil,header=T)
	assoc = assoc[complete.cases(assoc),]
	assoc = assoc[FDRcontrol(assoc$P, FDR=0.05),]
	cmh.assoc = rbind(assoc, cmh.assoc)
}
#ran 0.1 and 0.05


# Unphased haplotype  Data  -------------------	
#All associated SNPs are saved in Dhap*.assoc.hap and were unphased

hap.assoc = c()
fils = list.files(pattern="^Dhap*")
for(fil in fils){
	#fil = fils[13]
	assoc = read.table(file=fil,header=T)
	omni.assoc = assoc[assoc$HAPLOTYPE=="OMNIBUS",]
	assoc = assoc[assoc$HAPLOTYPE!="OMNIBUS",]
	omni.assoc = omni.assoc[FDRcontrol(omni.assoc$P, FDR=0.1),]
	assoc = assoc[assoc$LOCUS %in% omni.assoc$LOCUS,]
	assoc = assoc[FDRcontrol(assoc$P, FDR=0.1),]
	hap.assoc = rbind(hap.assoc, assoc)
}


# Unphased Permutations  Data  -------------------	
#All associated SNPs are saved in PERMUNhap_*.qassoc and were unphased




perm.assoc = c()
fils = list.files(pattern="^PERMUNhap_*")
fils = fils[grep("*.perm", fils)]

for(fil in fils){
	#fil = fils[13]
	assoc = read.table(file=fil,header=T)
	assoc = assoc[complete.cases(assoc),]
	assoc = assoc[FDRcontrol(assoc$EMP1, FDR=0.05),]
	perm.assoc = rbind(perm.assoc,assoc)
}

permq.assoc = c()
fils = list.files(pattern="^PERMUNhap_*")
fils = fils[grep("*.qassoc$", fils)]
for(fil in fils){
	#fil = fils[13]
	assoc = read.table(file=fil,header=T)
	permq.assoc = rbind(permq.assoc,assoc)
}

permq.assoc = permq.assoc[c(1,2,3)]
perm.assoc = merge(perm.assoc, permq.assoc, by="SNP")





save.image(file="AssociatedSNPTests.RData")

######################################
#Phased - THis looks like a good way to go.
	#Phased data with overlap of P valeues (tried 0.01 and 0.001)

# CMH Data  -------------------	
#All associated SNPs are saved in CLUSTEREDMAF_*.cmh and were phased




cmh.assoc = c()
fils = list.files(pattern="^CLUSTEREDMAF_")
for(fil in fils){
	#fil = fils[13]
	assoc = read.table(file=fil,header=T)
	assoc = assoc[complete.cases(assoc),]
	#assoc = assoc[FDRcontrol(assoc$P, FDR=0.1),]
	assoc = assoc[assoc$P<0.01,]
	cmh.assoc = rbind(assoc, cmh.assoc)
}


# Phased haplotype  Data  -------------------	
#All associated SNPs are saved in Dhap*.assoc.hap and were unphased

hap.assoc = c()
fils = list.files(pattern="^hap*")
for(fil in fils){
	#fil = fils[13]
	assoc = read.table(file=fil,header=T)
	#assoc = assoc[FDRcontrol(assoc$P, FDR=0.1),]
	assoc = assoc[assoc$P<0.01,]
	hap.assoc = rbind(hap.assoc, assoc)
}


# Phased Permutations  Data  -------------------	
#All associated SNPs are saved in PERMhap_*.qassoc and were phased

perm.assoc = c()
fils = list.files(pattern="^PERMhap_*")
fils = fils[grep("*.perm", fils)]

for(fil in fils){
	#fil = fils[13]
	assoc = read.table(file=fil,header=T)
	assoc = assoc[complete.cases(assoc),]
	#assoc = assoc[FDRcontrol(assoc$EMP1, FDR=0.1),]
	assoc = assoc[assoc$EMP1<0.01,]
	perm.assoc = rbind(perm.assoc,assoc)
}

permq.assoc = c()
fils = list.files(pattern="^PERMhap_*")
fils = fils[grep("*.qassoc$", fils)]
for(fil in fils){
	#fil = fils[13]
	assoc = read.table(file=fil,header=T)
	permq.assoc = rbind(permq.assoc,assoc)
}

permq.assoc = permq.assoc[c(1,2,3)]
perm.assoc = merge(perm.assoc, permq.assoc, by="SNP")



#load RPart associations
rpart.assoc=read.table(file="RPARTAssocSNPs",header=F); names(rpart.assoc)=c("CHR", "SNP", "v3", "BP", "IMP")











#merging them together -> how much overlap in SIG (P<0.01) snps?
snp.df = data.frame("SNP"=unlist(strsplit(unlist(as.vector(hap.assoc$SNPS)), "[|]")))
test = merge(cmh.assoc, perm.assoc, by ="SNP")
nrow(test[test$SNP %in% snp.df$SNP,])
test = test[test$SNP %in% snp.df$SNP,]

table(test$CHR.x)



test1=test[test$CHR!="1",]
test1[test1$SNP %in% rpart.assoc$SNP,]


13.5:506219  13  2709902  T 0.1098  C  7.338 0.0067500 0.02683 1.6100
16.8:396292  16  6320817  A 0.1098  G  7.338 0.0067500 0.11760 1.0220
3.16:726396   3 11807633  T 0.1829  C  8.164 0.0042740 0.11560 0.8191
9.12:226487  




#args = commandArgs(trailingOnly = TRUE)
pedfile = "DroneSamps_16.ped"
mapfile = "DroneSamps_16.map"



map=read.table(file=mapfile)
dt = read.table(pedfile, colClasses=c(rep("character",nrow(map)+6)))
pheno=read.table(file="DronePhenoHB.txt",header=F)
pheno=pheno[,-c(2)]
vec=merge(pheno,dt,by="V1")
pop=read.table(file="SELCON.txt",header=F);pop=pop[,-c(2)]
vec=vec[,-c(3:7)]
names(vec)[2]="pheno"
nms=vec$V1
vec=vec[,-1]
phens=vec$pheno
vec=vec[,-1]

print("loaded") #I will do this again with haplotypes and NOT genotypes later (i.e. ranomdly select columns)

c1=vec[,seq(1,ncol(vec),2)]
c2=vec[,seq(2,ncol(vec),2)]

pfun <- function(x, y) paste(x, y, sep = "")
datNew <- vector("list", ncol(c1)/2)
for (i in 1:ncol(c1)) {
    datNew[[i]] <- pfun(c1[[i]], c2[[i]])
}

vec <- as.data.frame(datNew);names(vec)=seq(1:ncol(vec))
vec$pheno=phens
vec$fam=pop$V3
names(vec)=paste("SNP",names(vec),sep="")








which(map$V2=="16.8:396292")
boxplot(vec$SNPpheno~vec$SNP10864*vec$SNPfam)















plot(test$BP.x[test$CHR=="12"], 
	-log10(test$EMP1[test$CHR=="12"]),
	 col ="red")
	
	
par(new=TRUE)
plot(movingAverage(cmh.assoc$BP, n=500), 
	movingAverage(-log10(cmh.assoc$P),n=500),
	ylim =c(0,1),
	xlim=c(0,30000000))





















#let's test CHR6:
source("/media/data1/forty3/brock/scripts/movingavg.r")

cmh.assoc = c()
fils = list.files(pattern="^CLUSTEREDMAF_")
for(fil in fils){
	#fil = fils[12]
	assoc = read.table(file=fil,header=T)
	assoc = assoc[complete.cases(assoc),]
	assoc = assoc[FDRcontrol(assoc$P, FDR=0.05),]
	cmh.assoc = rbind(assoc, cmh.assoc)
}
#ran 0.1 and 0.05



# Phased Permutations  Data  -------------------	
#All associated SNPs are saved in PERMhap_*.qassoc and were phased

perm.assoc = c()
fils = list.files(pattern="^PERMhap_*")
fils = fils[grep("*.perm", fils)]

for(fil in fils){
	#fil = fils[12]
	assoc = read.table(file=fil,header=T)
	assoc = assoc[complete.cases(assoc),]
	assoc = assoc[FDRcontrol(assoc$EMP1, FDR=0.05),]
	perm.assoc = rbind(perm.assoc,assoc)
}

permq.assoc = c()
fils = list.files(pattern="^PERMhap_*")
fils = fils[grep("*.qassoc$", fils)]
for(fil in fils){
	#fil = fils[12]
	assoc = read.table(file=fil,header=T)
	permq.assoc = rbind(permq.assoc,assoc)
}

permq.assoc = permq.assoc[c(1,2,3)]
perm.assoc = merge(perm.assoc, permq.assoc, by="SNP")

#

perm.assoc=perm.assoc[with(perm.assoc, order(BP)),]
cmh.assoc=cmh.assoc[with(cmh.assoc, order(BP)),]


plot(movingAverage(perm.assoc$BP, n=500), 
	movingAverage(-log10(perm.assoc$EMP1),n=500),
	ylim =c(0,1),
	xlim=c(0,30000000), col ="red")
par(new=TRUE)
plot(movingAverage(cmh.assoc$BP, n=500), 
	movingAverage(-log10(cmh.assoc$P),n=500),
	ylim =c(0,1),
	xlim=c(0,30000000))








