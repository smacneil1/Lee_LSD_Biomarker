# ---
# title: "ssGSEA"
# author: "Shelley MacNeil"
# date: "July 29, 2016"
# purpose: correlate binreg results with ICBP/TCGA data
# ---- 

source("~/Documents/PhDProjects/GFRN_signatures/Key_ASSIGN_functions_balancedsig.R")
source("~/Documents/PhDProjects/PI3K_Signature_Validation/scripts/CommonRFunctions.R")

# read in ICBP drug and protein data
icbp_drug_response = "~/Documents/PhDProjects/Multipathway_Modeling/data/ICBP_drugs.txt"
ICBP_drugs<-read.delim(icbp_drug_response, header=1, sep='\t',row.names=1, stringsAsFactors = FALSE)
ICBP_drugs=ICBP_drugs[,11:100]
icbp_protein = "~/Documents/PhDProjects/Multipathway_Modeling/data/ICBP_proteomics.txt"
ICBP_protein<-read.delim(icbp_protein, header=1, sep='\t',row.names=1, stringsAsFactors = FALSE)

# read in the Lee files 
setwd("~/Documents/PhDProjects/Lee_LSD_Biomarker/results/binReg/Lee/")
filenames=system("ls", intern=TRUE)
filenames

for(i in 1:length(filenames))
{
   f <- read.delim(filenames[i], quote="", skip= 33, row.names = 2, header=F)
   colnames(f)
   rownames(f)
   dim(f)
   f[,3]
   f= f[5]
   f
   #rownames(f)=f[,1]
   colnames(f)=c(gsub('\\..*$','',filenames[i]))
  if(i!=1){
    lee<-cbind(lee,f)
  }
  else {
    lee<-f
  }
}

View(lee)

# read in the LSD files 
setwd("~/Documents/PhDProjects/Lee_LSD_Biomarker/results/binReg/LSD/")

# read in one file at at time, and cbind all columns
filenames=system("ls", intern=TRUE)
filenames

for(i in 1:length(filenames))
{
  f <- read.delim(filenames[i], quote="", skip= 33, row.names = 2, header=F)
  f[,3]
  f= f[5]
  colnames(f)=c(gsub('\\..*$','',filenames[i]))
  if(i!=1){
    lsd<-cbind(lsd,f)
  }
  else {
    lsd<-f
  }
}
View(lsd)

# correlate with the ICBP drug data
icbp_drug_lsd=merge_drop(lsd,ICBP_drugs)
icbp_drug_lee=merge_drop(lee,ICBP_drugs)
drug_cor_lsd = cor(icbp_drug_lsd[,1:7], icbp_drug_lsd[,8:97], use="pairwise",method="spearman") 
drug_cor_lee = cor(icbp_drug_lee[,1:7], icbp_drug_lee[,8:97], use="pairwise",method="spearman") 
write.table(drug_cor_lsd, "~/Documents/PhDProjects/Lee_LSD_Biomarker/results/binReg/Cors/ICBP_drug_cors_LSD.txt", col.names=NA,sep='\t',quote=F)
write.table(drug_cor_lee, "~/Documents/PhDProjects/Lee_LSD_Biomarker/results/binReg/Cors/ICBP_drug_cors_LEE.txt", col.names=NA,sep='\t',quote=F)


# correalte with the ICBP protein data
icbp_protein_prob_lsd=merge_drop(lsd, ICBP_protein)
icbp_protein_prob_lee=merge_drop(lee, ICBP_protein)
prot_cor_lsd = cor(icbp_protein_prob_lsd[,1:7], icbp_protein_prob_lsd[,8:77], method="spearman") 
prot_cor_lee = cor(icbp_protein_prob_lee[,1:7], icbp_protein_prob_lee[,8:77] ,method="spearman") 
write.table(prot_cor_lsd, "~/Documents/PhDProjects/Lee_LSD_Biomarker/results/binReg/Cors/ICBP_prot_cors_LSD.txt", col.names=NA,sep='\t',quote=F)
write.table(prot_cor_lee, "~/Documents/PhDProjects/Lee_LSD_Biomarker/results/binReg/Cors/ICBP_prot_cors_LEE.txt", col.names=NA,sep='\t',quote=F)

# END

-------------------------
# finish this function later

correlate_binReg_ICBP = function(filedir, name, outFileDrugCors) {

setwd(filedir)
filenames=system("ls", intern=TRUE)
filenames

for(i in 1:length(filenames)) {
  f <- read.delim(filenames[i], quote="", skip= 33, row.names = 2, header=F)
  colnames(f)
  rownames(f)
  dim(f)
  f[,3]
  f= f[5]
  f
  #rownames(f)=f[,1]
  colnames(f)=c(gsub('\\..*$','',filenames[i]))
  if(i!=1){
    name<-cbind(name,f)
  }
  else {
    name<-f
  }
}

icbp_drug_prob=merge_drop(name,ICBP_drugs)
icbp_drug_cor = cor(icbp_drug_prob[,1:7], icbp_drug_prob[,8:97], use="pairwise",method="spearman") 
write.table(icbp_drug_cor, outFileDrugCors, col.names=NA,sep='\t',quote=F)

icbp_protein_prob=merge_drop(name,TCGA.BRCA.RBN)



}


