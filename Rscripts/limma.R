# ---
# title: "ssGSEA"
# author: "Shelley MacNeil"
# date: "July 29, 2016"
# purpose: run a differential expression analysis with limma 
# ---- 

#load in libraries 
library(limma)

#specficy file names 
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

# name the outfiles 
limma_lee_outfile="~/Documents/PhDProjects/Lee_LSD_Biomarker/results/limma/limma_Lee.txt"
limma_lsd_outfile="~/Documents/PhDProjects/Lee_LSD_Biomarker/results/limma/limma_LSD.txt"


# Make limma function
limma <- function (inputData, classFile, outFile){
  
#make the model matrix
design <- model.matrix(~ factor(classFile[,2]))
design
colnames(design) <- c("ALL", "CtrlvsTreated")  

# run the DE analysis
fit <- lmFit(inputData, design)
fit  <- eBayes(fit)
DEresults <- topTable(fit, coef="CtrlvsTreated", number=Inf, adjust.method = "BH", sort.by = "p" )
  
#write the file  
write.table(DEresults, outFile, col.names=NA, quote=FALSE, sep="\t")

DEgenes_01 <- DEresults[ which(DEresults$adj.P.Val < 0.01),]
DEgenes_05 <- DEresults[ which(DEresults$adj.P.Val < 0.05),]
DEgenes_50 <- topTable(fit, coef="CtrlvsTreated", number=50, adjust.method = "BH", sort.by = "p" )


print("Number of genes with FDR < 0.01")
print(dim(DEgenes_01))
print("Number of genes with FDR < 0.05")
print(dim(DEgenes_05))

print("Complete")
return(DEgenes_50)
}


View(Lee_Ctrl_Both)

dim(Lee_Ctrl_Both)
limma_lee_05 = limma(Lee_Ctrl_Both, class_Lee, limma_lee_outfile)
View(limma_lee_05)

limma_lee_01 = limma(Lee_Ctrl_Both, class_Lee, limma_lee_outfile)
View(limma_lee_01)

# "Number of genes with FDR < 0.01"  --- 154   6
# "Number of genes with FDR < 0.05"  --- 439   6

limma_lsd_05 = limma(LSD_Ctrl_Both, class_LSD, limma_lee_outfile)
View(limma_lsd)
limma_lsd_50 = limma(LSD_Ctrl_Both, class_LSD, limma_lee_outfile)

#"Number of genes with FDR < 0.01"   ----  3  6
# "Number of genes with FDR < 0.05"  ---- 14  6


# save the limma results 
save.image(file="~/Documents/PhDProjects/Lee_LSD_Biomarker/data/Lee_LSD_limma.RData")




