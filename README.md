# Workflow and Supplementary Material for Drone Selection Experiment









##VCF Creation:
- VCFCreation_DroneSelection.txt
- DroneVCFcreation.txt




##Association Analyses
1. DronePLINKset.sh identifies high FST across the genome between selected and control populations, then 
2. Pass FST to DroneFST.r and identify high FST regions within the genome- outputs all SNPs within those regions (CandidateSNPs.snp)
3. DronePLINKset.sh uses a modified CandidateSNPs.snp (Candidates98_ALL.set ) to run the PLINK set test within each high FST region
4. DronePLINKset.sh uses  CandidateSNPs.snp to dientify independant SNPs within each region (called CandidateSET.xxx)
5. CandidateSET.xxx is then used for Recursive Partitioning and Regression Tree (RpartAssociations.r) 
