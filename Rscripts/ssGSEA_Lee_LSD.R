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

#write the file 
write.table(geneSets, outFile_geneSets, col.names=NA, quote=FALSE, sep="\t")
return(geneSets)

}

#specify paths and file names
C2gmtFile <- "~/Documents/PhDProjects/GMT_Files/c2.cp.v4.0.symbols.gmt"
C6gmtFile <- "~/Documents/PhDProjects/GMT_Files/c6.all.v5.1.symbols.gmt"

#output files 
outfile_ssGSVA_Lee_C2= "~/Documents/PhDProjects/Lee_LSD_Biomarker/results/ssGSEA/ssGSEA_enrichment_Lee_C2.txt"
outfile_geneSets_Lee_C2="~/Documents/PhDProjects/Lee_LSD_Biomarker/results/ssGSEA/ssGSEA_geneSets_Lee_C2.txt"
outfile_ssGSVA_Lee_C6= "~/Documents/PhDProjects/Lee_LSD_Biomarker/results/ssGSEA/ssGSEA_enrichment_Lee_C6.txt"
outfile_geneSets_Lee_C6="~/Documents/PhDProjects/Lee_LSD_Biomarker/results/ssGSEA/ssGSEA_geneSets_Lee_C6.txt"

outfile_ssGSVA_LSD_C6 <- "~/Documents/PhDProjects/Lee_LSD_Biomarker/results/ssGSEA/ssGSEA_enrichment_LSD_C6.txt"
outfile_geneSets_LSD_C6 <- "~/Documents/PhDProjects/Lee_LSD_Biomarker/results/ssGSEA/ssGSEA_geneSets_LSD_C6.txt"
outfile_ssGSVA_LSD_C2 <- "~/Documents/PhDProjects/Lee_LSD_Biomarker/results/ssGSEA/ssGSEA_enrichment_LSD_C2.txt"
outfile_geneSets_LSD_C2 <- "~/Documents/PhDProjects/Lee_LSD_Biomarker/results/ssGSEA/ssGSEA_geneSets_LSD_C2.txt"

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

Lee_ssGSEA_C2= ssGSEA(Lee_Ctrl_Both, class_Lee, C2gmtFile, outfile_ssGSVA_Lee_C2, outfile_geneSets_Lee_C2)
Lee_ssGSEA_C6= ssGSEA(Lee_Ctrl_Both, class_Lee, C6gmtFile, outfile_ssGSVA_Lee_C6, outfile_geneSets_Lee_C6)

# LSD
LSD_ssGSEA_C2 = ssGSEA(LSD_Ctrl_Both, class_LSD,C2gmtFile,outfile_ssGSVA_LSD_C2, outfile_geneSets_LSD_C2)
LSD_ssGSEA_C6 = ssGSEA(LSD_Ctrl_Both, class_LSD,C6gmtFile,outfile_ssGSVA_LSD_C6, outfile_geneSets_LSD_C6)







  