#pFST for FST p-values



#https://github.com/jewmanchue/vcflib/wiki

./pFst --target 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,18,19,20,21,22,23,24,25,26,27,28,29 --background 16,17,30,31,32,33,34,35,36,37,38,39,40,41 --deltaaf 0.0 --file /media/data1/forty3/drone/vcf_drone/DroneSelectionFinal.recode.vcf --counts --type PL   > /media/data1/forty3/drone/FST/pFST/pFstDRONE.counts

#NOT APPLICABLE ANALYSES:
./pFst --target 1,2,3,4,5,6,7,8,9  --background 10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29  --deltaaf 0.0 --file /media/data1/afz/VCF/AllAMCSNPs.recode.vcf --counts --type PL   > /media/data1/forty3/drone/FST/pFST/C.counts

./pFst --target 19,20,21,22,23,24,25,26,27,28,29  --background 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18  --deltaaf 0.0 --file /media/data1/afz/VCF/AllAMCSNPs.recode.vcf --counts --type PL   > /media/data1/forty3/drone/FST/pFST/A.counts

./pFst --target 10,11,12,13,14,15,16,17,18  --background 1,2,3,4,5,6,7,8,9,19,20,21,22,23,24,25,26,27,28,29 --deltaaf 0.0 --file /media/data1/afz/VCF/AllAMCSNPs.recode.vcf --counts --type PL   > /media/data1/forty3/drone/FST/pFST/M.counts
