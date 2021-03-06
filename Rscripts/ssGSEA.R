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


#create the ssGSEA function
ssGSEA <- function (inputData, classFile, gmtFile, outFile_enrichment, outFile_geneSets){

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
return(geneSets)
#write the file 
write.table(geneSets, outFile_geneSets, col.names=NA, quote=FALSE, sep="\t")

}


#specify paths and file names
C2gmtFile <- "~/Documents/PhDProjects/Data/GMT_Files/c2.cp.v4.0.symbols.gmt"
#C6gmtFile <- "~/Documents/PhDProjects/Data/GMT_Files/"

#output files 
outfile_ssGSVA_Lee = "~/Documents/PhDProjects/Lee_LSD_Biomarker/results/ssGSEA/ssGSEA_enrichment_Lee.txt"
outfile_ssGSVA_Lee 
outfile_geneSets_Lee="~/Documents/PhDProjects/Lee_LSD_Biomarker/results/ssGSEA/ssGSEA_geneSets_Lee.txt"
outfile_geneSets_Lee
outfileFile_ssGSVA_LSD <- "~/Documents/PhDProjects/Lee_LSD_Biomarker/results/ssGSEA/ssGSEA_enrichment_LSD.txt"
outfileFile_ssGSVA_LSD
outfileFile_geneSets_LSD <- "~/Documents/PhDProjects/Lee_LSD_Biomarker/results/ssGSEA/ssGSEA_geneSets_LSD.txt"
outfileFile_geneSets_LSD


#load in the R data from Combine_Single_RNAseq_Files.Rmd 
Lee_LSD_Data = load("data/Lee_LSD_RNAseqData.RData")
colnames(Lee_Ctrl_Both)
colnames(LSD_Ctrl_Both)

#make Lee class file
class_Lee=as.matrix(c(rep("control", 8),rep("Lee", 8) ))
class_Lee=cbind(colnames(Lee_Ctrl_Both), class_Lee) 
class_Lee

#make LSD Class File
class_LSD=as.matrix(c(rep("control", 8),rep("LSD", 8) ))
class_LSD=cbind(colnames(LSD_Ctrl_Both), class_LSD) 
class_LSD

#Run the ssGSEA function

# inputData, classFile, gmtFile, outFile_enrichment, outFile_geneSets

#Lee
Lee_ssGSEA_C2= ssGSEA(Lee_Ctrl_Both, class_Lee, C2gmtFile, outfile_ssGSVA_Lee, outfile_geneSets_Lee)

View(Lee_ssGSEA_C2)

LSD_ssGSEA_C2 = ssGSEA(LSD_Ctrl_Both, class_LSD,C2gmtFile,outputFile_ssGSVA_LSD, outputFile_geneSets_LSD  )




  