# ---
# title: "ssGSEA"
# author: "Shelley MacNeil"
# date: "July 25, 2016"
# purpose: run ssGSEA
# ---- 

# import packages
library(GSEABase)
library(GSVA)
#library(GSVAdata)
library(limma)
#library(gage)
#data(c2BroadSets)
library(Biobase)
library(genefilter)

data= commandArgs()[7]
print (data)
outfile_ssGSVA_enrich_C2 = commandArgs()[8]
outfile_ssGSVA_genesets_C2 = commandArgs()[9]
outfile_ssGSVA_enrich_C6 = commandArgs()[10]
outfile_ssGSVA_genesets_C2 = commandArgs()[11]

#classFilePath = commandArgs()[8]
#print(classFilePath)

#specify paths and file names
C2gmtFile <- "~/Documents/PhDProjects/GMT_Files/c2.cp.v4.0.symbols.gmt"
C6gmtFile <- "~/Documents/PhDProjects/GMT_Files/c6.all.v5.1.symbols.gmt"

#read in the R data 
data= as.matrix(read.table(data, stringsAsFactors = FALSE, header=FALSE,sep='\t', quote="", row.names = 1))
dim(data)
data=data[apply(data[,1:26]==0,1,mean) < 0.85,]


#create the ssGSEA function
ssGSEA <- function (inputData, gmtFile, outFile_enrichment, outFile_geneSets){
  
  # Within GSEABase, function converts gmt files to an object that has those gene sets (convert gene set collection class, because
  # the c2BroadSets gene set collection class uses gene IDs instead of gene symbols); use the getGMT
  gmtSet <- getGmt(gmtFile)
  
  # Run GSVA in ssGSEA format and output to file
  gsva <- gsva(as.matrix(inputData), gmtSet, min.sz=10, max.sz=500, verbose=TRUE, rnaseq=TRUE, method="ssgsea")
  View(gsva)
  
  #write the ssGSEA results
  write.table(gsva, outFile_enrichment, col.names=NA, quote=FALSE, sep="\t")
  
  #make the model matrix
  design <- model.matrix(~ factor(classFile[,2]))
  design
  colnames(design) <- c("ALL", "CtrlvsTreated")
  
  # run the DE analysis
  fit <- lmFit(gsva, design)
  View(fit)
  fit  <- eBayes(fit)
  geneSets <- topTable(fit, coef="CtrlvsTreated", number=Inf, adjust.method = "BH", sort.by = "p" )
  
  #write the file 
  write.table(geneSets, outFile_geneSets, col.names=NA, quote=FALSE, sep="\t")
  return(geneSets)
  
}

#Run the ssGSEA function
# inputData, classFile, gmtFile, outFile_enrichment, outFile_geneSets
ssGSEA_C2= ssGSEA(data, C2gmtFile, outfile_ssGSVA_C2, outfile_geneSets_C2)
print("Done with C2")
ssGSEA_C6= ssGSEA(data, C6gmtFile, outfile_ssGSVA_C6, outfile_geneSets_C6)
print("Done with C6")
print("Complete")

#output files 
#outfile_ssGSVA_C2= "~/Documents/PhDProjects/Lee_LSD_Biomarker/results/ssGSEA/ssGSEA_enrichment_TCGA_BRCA_C2.txt"
#outfile_ssGSVA_C6= "~/Documents/PhDProjects/Lee_LSD_Biomarker/results/ssGSEA/ssGSEA_enrichment_TCGA_BRCA_C6.txt"
#outfile_geneSets_C2="~/Documents/PhDProjects/Lee_LSD_Biomarker/results/ssGSEA/ssGSEA_geneSets_TCGA_BRCA_C2.txt"
#outfile_geneSets_C6="~/Documents/PhDProjects/Lee_LSD_Biomarker/results/ssGSEA/ssGSEA_geneSets_TCGA_BRCA_C6.txt"
#inFile 
#inFile="~/Documents/PhDProjects/Multipathway_Modeling/Data/PANCAN24_BRCA_1119_TPMlog2.txt"








