###
# Drone PLINK Set Tests
###



#Runs FST with VCFtools between 2 populations, then outputs SNPs in high FST windows and looks for associated SNPs within those areas















#####
#Fst Analyses with VCFTools
	#Calculate FST  all selected VS control (2 vs 1+3=Selpop.txt)
vcftools --maf 0.05 --vcf DroneSelectionFinal.recode.vcf --weir-fst-pop pop2.txt --weir-fst-pop Selpop.txt --out pop2_vs_sel 

#####
#Analyses with R for getting out high regions:

#Run DroneFST.r
	#Outputs CandidateSNPs98.snps, I updated that by hand to Candidates98_ALL.set for set tests



#####
#Association Analyses
	
	
	
#Select and Control together
	
	#For only Hygiene:
	plink --file DroneSelection --set-test --set Candidates98_ALL.set --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --linear --mperm 10000 --out AllBees98Candidates --set-p 0.01 --set-max 5
	plink --file DroneSelection --set-test --set Candidates98_ALL.set --mpheno 1 --pheno /media/data1/forty3/drone/vcf_drone/DronePhenoHB.txt --noweb --linear --mperm 10000 --out AllBees98Candidates001 --set-p 0.001 --set-max 5

		# AllBees98Candidates001 and  AllBees98Candidates are P<0.001 and P<0.01 for 5 maximum SNPs after 10000 permutations 
		# Both files use all bees for the hygiene phenotype. 
		
		
		#For all expression candidates:
	
	for i in {1..9}
	do
	   plink --file DroneSelection --set-test --set Candidates98_ALL.set --mpheno $i --pheno /media/data1/forty3/drone/vcf_drone/HBexpression.txt --noweb --linear --mperm 10000 --out ALL_PHENO$i --set-p 0.001 --set-max 5
	done
	
	
	
	

	
	
####
#List of Independant SNPs (taken from RpartCandidateTest.r)
	#For Rpart analyses.

###	
#Get out LINKED SNPs within Sets:	
	#plink --file ControlHighFST  --ld-snp-list highFST.list  --noweb --allow-no-sex --out CONTROL
	#plink --file SelHighFST  --ld-snp-list highFST.list  --noweb --allow-no-sex --out SELECT
	#plink --file DroneSelection --indep-pairwise 50 5 0.5 --extract highFST.list --noweb --allow-no-sex  --keep controlBees.txt --out CONTROL
	#plink --file DroneSelection --indep-pairwise 50 5 0.5 --extract highFST.list --noweb --allow-no-sex  --remove controlBees.txt --out SELECT

plink --file DroneSelection --indep-pairwise 50 5 0.5 --extract highFST.list --noweb --allow-no-sex --out ALL
plink --file DroneSelection --extract ALL.prune.in --out CandidateSET --noweb --make-bed 
plink --bfile CandidateSET  --recode --out  CandidateSET --noweb
	
	#this will output a set of PLINK called CandidateSET.xxx
	
	
	

	
#Plotting these data:	
	#source(DronePlot1_ChromosomeFST.r)
		#outputs FSTSelectTest.Chromosome*_*_updated.pdf
	
	
	