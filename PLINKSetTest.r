####################################################
# PLINK Set Test 
#######
#file="ControlHighFST.map"
#http://pngu.mgh.harvard.edu/~purcell/plink/anal.shtml#set

plink --file DroneSelection --extract CandidateSNPs.snp --out ControlHighFST --noweb --make-bed --keep controlBees.txt 
plink --bfile ControlHighFST  --recode --out ControlHighFST --noweb

plink --file ControlHighFST --set-test --set highFST.set --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --allow-no-sex --assoc --mperm 1000 --out CONTROLSET --set-p 0.01 --set-max 10 --qt-means


plink --file DroneSelection --extract CandidateSNPs.snp --out SelHighFST --noweb --make-bed --remove controlBees.txt 
plink --bfile SelHighFST  --recode --out SelHighFST --noweb

plink --file SelHighFST --set-test --set highFST.set --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --allow-no-sex --assoc --mperm 1000 --out SelSET --set-p 0.01 --set-max 10 --qt-means
