




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



