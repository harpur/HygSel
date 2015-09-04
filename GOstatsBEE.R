###
# GO Analysis for bees
###

#Set up right now for Apis, but can be changed to w/e is needed.
#uses drosophila orthologs and FBGN.
#BAH 25-Aug-15


# Important background --------------------
#This is where most/all of the code below was assembled from. Read if anything is unclear.
#http://www.bioconductor.org/packages/release/bioc/vignettes/GOstats/inst/doc/GOstatsForUnsupportedOrganisms.pdf
#https://www.bioconductor.org/help/workflows/annotation-data/
#http://www.bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf

# Libraries and functions ------------------
source("http://bioconductor.org/biocLite.R")
biocLite("GOstats")
library("GOstats")
library("AnnotationForge")
library("GSEABase")
#available.dbschemas() #lists available data packages
library("drosophila2.db") #Download a database from those that are available
library("KEGG.db")

#load database for GO analyses and prepare GO data.frame and GeneSetCollection ---------------
k = (keys(drosophila2.db,keytype="GO")) # I read in our dros. orthologs
#orths = unique(unlist(read.table(file="clipboard")$V1))
goframeData = select(drosophila2.db, keys=k, cols=c("EVIDENCE", "FLYBASE"), keytype="GO")
goframeData = goframeData[goframeData$FLYBASE %in% orths, ]
goframeData = goframeData[!duplicated(goframeData$FLYBASE), ]
names(goframeData) =c("frame.go_id", "frame.Evidence","frame.ontology", "frame.gene_id")
goframeData$frame.ontology = NULL
#goframeData.BP = goframeData[goframeData$frame.ontology =="BP",]
#goframeData.BP$frame.ontology = NULL
goFrame = GOFrame(goframeData,organism="Apis mellifera")
goAllFrame = GOAllFrame(goFrame) #this is dropping my geneIDs.
gsc = GeneSetCollection(goAllFrame, setType = GOCollection())
	

# Prepare gene universe and your test set --------------------------------
universe = goframeData$frame.gene_id #define your gene universe
genes = sample(universe,567,replace=FALSE)#define your genes of interest
#genes = as.character(unique(unlist(read.table(file="clipboard")$V1)))

#  Perform tests --------------------------------
#Kindly provded by Sasha (thanks, man)
outcomeBP <- hyperGTest(GSEAGOHyperGParams(name = "lrt", 
			geneSetCollection=gsc,
			geneIds = as.character(genes),
			universeGeneIds = as.character(universe),
			ontology = "BP",
			pvalueCutoff = 0.05,
			conditional = FALSE,
			testDirection = "over"))
outcomeMF <- hyperGTest(GSEAGOHyperGParams(name = "lrt", 
			geneSetCollection=gsc,
			geneIds = as.character(genes),
			universeGeneIds = as.character(universe),
			ontology = "MF",
			pvalueCutoff = 0.05,
			conditional = FALSE,
			testDirection = "over"))
outcomeCC <- hyperGTest(GSEAGOHyperGParams(name = "lrt", 
			geneSetCollection=gsc,
			geneIds = as.character(genes),
			universeGeneIds = as.character(universe),
			ontology = "CC",
			pvalueCutoff = 0.05,
			conditional = FALSE,
			testDirection = "over"))
#head(summary(outcomeBP))
#head(summary(outcomeMF))			
#head(summary(outcomeCC))
		



		
#  KEGG Pathways --------------------------------			
	#		
	#k = (keys(drosophila2.db,keytype="FLYBASE")) # I read in our dros. orthologs
	##orths = unique(unlist(read.table(file="clipboard")$V1))
	#keggframeData = select(drosophila2.db, keys=k, cols=c("PATH"), keytype="FLYBASE")
	#keggframeData = keggframeData[keggframeData$FLYBASE %in% orths,]
	#keggframeData = keggframeData[!is.na(keggframeData$PATH),]
	#keggFrame=KEGGFrame(keggframeData,organism="Apis mellifera")
	#
	#
	#
	#goframeData = goframeData[goframeData$FLYBASE %in% orths, ]
	#goframeData = goframeData[!duplicated(goframeData$FLYBASE), ]
	#names(goframeData) =c("frame.go_id", "frame.Evidence","frame.ontology", "frame.gene_id")
	#goframeData$frame.ontology = NULL
	##goframeData.BP = goframeData[goframeData$frame.ontology =="BP",]
	##goframeData.BP$frame.ontology = NULL
	#goFrame = GOFrame(goframeData,organism="Apis mellifera")
	#goAllFrame = GOAllFrame(goFrame) #this is dropping my geneIDs.
	#gsc = GeneSetCollection(goAllFrame, setType = GOCollection())
	#
	#
	#		
	#			
	#> frame = toTable(org.Hs.egPATH)
	#> keggframeData = data.frame(frame$path_id, frame$gene_id)
	#> head(keggframeData)
	#
	#			
				
			
			
# I'll try REVIGO for a plot:





