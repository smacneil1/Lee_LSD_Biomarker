# ---
# title: "ssGSEA"
# author: "Shelley MacNeil"
# date: "August 1, 2016"
# purpose: See if EMT pathways in TCGA correlate with growth factor pathways
# ---- 


# load the libraries 

# read in the data
TCGA_ssGSEA_path="~/Documents/PhDProjects/Lee_LSD_Biomarker/results/TCGA_ssGSEA/ssGSEA_enrich_TCGA_BRCA_all.txt"
TCGA_ssGSEA=as.matrix(read.table(TCGA_ssGSEA_path, row.names =1, header= FALSE, skip=1))
TCGA_ssGSEA_t=t(TCGA_ssGSEA)

EMT_RTK= TCGA_ssGSEA_t[]

colnames(TCGA_ssGSEA_t)
row.names(TCGA_ssGSEA_t)
head(TCGA_ssGSEA_t)
dim(TCGA_ssGSEA_t)

SAM_RTK = as.matrix(TCGA_ssGSEA_t[ , grep("SAM_RTKS", colnames(TCGA_ssGSEA_t))])

Anastassiou_EMT = as.matrix(TCGA_ssGSEA_t[ , grep("ANASTASSIOU_CANCER_MESENCHYMAL_TRANSITION_SIGNATURE", colnames(TCGA_ssGSEA_t))])

JECHLINGER_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_UP=as.matrix(TCGA_ssGSEA_t[ , grep("JECHLINGER_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_UP", colnames(TCGA_ssGSEA_t))])
JECHLINGER_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_DN=as.matrix(TCGA_ssGSEA_t[ , grep("JECHLINGER_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_DN", colnames(TCGA_ssGSEA_t))])

Wang_Lsd_down=as.matrix(TCGA_ssGSEA_t[ , grep("WANG_LSD1_TARGETS_DN", colnames(TCGA_ssGSEA_t))])
Wang_Lsd_up=as.matrix(TCGA_ssGSEA_t[ , grep("WANG_LSD1_TARGETS_UP", colnames(TCGA_ssGSEA_t))])

FOSTER_KDM1A_TARGETS_DN=as.matrix(TCGA_ssGSEA_t[ , grep("FOSTER_KDM1A_TARGETS_DN", colnames(TCGA_ssGSEA_t))])
FOSTER_KDM1A_TARGETS_UP=as.matrix(TCGA_ssGSEA_t[ , grep("FOSTER_KDM1A_TARGETS_UP", colnames(TCGA_ssGSEA_t))])

EPIDERMAL_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY=as.matrix(TCGA_ssGSEA_t[ , grep("FOSTER_KDM1A_TARGETS_UP", colnames(TCGA_ssGSEA_t))])
TRANSMEMBRANE_RECEPTOR_PROTEIN_TYROSINE_KINASE_SIGNALING_PATHWAY=as.matrix(TCGA_ssGSEA_t[ , grep("TRANSMEMBRANE_RECEPTOR_PROTEIN_TYROSINE_KINASE_SIGNALING_PATHWAY", colnames(TCGA_ssGSEA_t))])

HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION=as.matrix(TCGA_ssGSEA_t[ , grep("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", colnames(TCGA_ssGSEA_t))])
HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION

# Correlate with eachother
# EMT and RTK
cor = cor(Anastassiou_EMT,SAM_RTK, method = "spearman")
ttest= t.test(Anastassiou_EMT, SAM_RTK)
reg1 <- lm(Anastassiou_EMT~SAM_RTK)
plot(Anastassiou_EMT~SAM_RTK, main = "Pathway Correlation in TCGA", xlab="SAM_RTK", ylab="Anastassiou_EMT")
abline(reg1)
title(sub=paste("P-value = ", format(ttest$p.value, digits=5)), adj=1, line=4, font=1)
title(sub=paste("Spearman = ", format(cor, digits=4)), adj=1, line=-11, font=1)

# create the plotting function
plot_cors = function(pathway1, pathway2, xlab, ylab) {
cor = cor(pathway1, pathway2, method = "spearman")
ttest= t.test(pathway1, pathway2)
print(ttest)
reg1 <- lm(pathway1~pathway2)
plot(pathway1~pathway2, main = "Pathway Correlation in TCGA", xlab=xlab , ylab=ylab)
abline(reg1)
#boxplot(pathway1, pathway2, main = "Enrichment Scores TCGA", xlab=xlab , ylab=ylab)
title(sub=paste("P-value (t.test) = ", format(ttest$p.value, digits=3)), adj=1, line=4, font=1)
title(sub=paste("R = ", format(cor, digits=3)), adj=1, line=-2, font=1) 
}

# put the second one first

# Sams's RTK and numerous EMT sigs

pdf("~/Documents/PhDProjects/Lee_LSD_Biomarker/results/TCGA_ssGSEA/TCGA_BRCA_ssGSEA_cors.pdf")

# Sams's RTK and numerous EMT sigs
plot_cors(HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, SAM_RTK,"SAM_RTK", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" )
plot_cors(Anastassiou_EMT, SAM_RTK, "SAM_RTK", "Anastassiou_EMT")
# not a high correlation between EMT up and RTK
plot_cors(JECHLINGER_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_UP,SAM_RTK,"SAM_RTK", "JECHLINGER_EMT_UP" )
#high coorelation between EMT down and RTK 
plot_cors(JECHLINGER_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_DN,SAM_RTK,"SAM_RTK", "JECHLINGER_EMT_DOWN" )

plot_cors(EPIDERMAL_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY, HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "EGFR_SIGNALING_PATHWAY" )
plot_cors(EPIDERMAL_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY, Anastassiou_EMT,"Anastassiou_EMT", "EGFR_SIGNALING_PATHWAY" )
plot_cors(EPIDERMAL_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY, JECHLINGER_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_UP,"JECHLINGER_EMT_UP", "EGFR_SIGNALING_PATHWAY" )
plot_cors(EPIDERMAL_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY, JECHLINGER_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_DN,"JECHLINGER_EMT_DN", "EGFR_SIGNALING_PATHWAY" )


plot_cors(TRANSMEMBRANE_RECEPTOR_PROTEIN_TYROSINE_KINASE_SIGNALING_PATHWAY, HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "RECEPTOR_PROTEIN_TYROSINE_KINASE_SIGNALING" )
plot_cors(TRANSMEMBRANE_RECEPTOR_PROTEIN_TYROSINE_KINASE_SIGNALING_PATHWAY, Anastassiou_EMT,"Anastassiou_EMT", "RECEPTOR_PROTEIN_TYROSINE_KINASE_SIGNALING" )
plot_cors(TRANSMEMBRANE_RECEPTOR_PROTEIN_TYROSINE_KINASE_SIGNALING_PATHWAY, JECHLINGER_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_UP,"JECHLINGER_EMT_UP", "RECEPTOR_PROTEIN_TYROSINE_KINASE_SIGNALING" )
plot_cors(TRANSMEMBRANE_RECEPTOR_PROTEIN_TYROSINE_KINASE_SIGNALING_PATHWAY, JECHLINGER_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_DN,"JECHLINGER_EMT_DN", "RECEPTOR_PROTEIN_TYROSINE_KINASE_SIGNALING" )


# EMT and KDM1A
# When KDM1A down negative correlates with EMT 
# when KDM1A up coreelates with EMT 
plot_cors(FOSTER_KDM1A_TARGETS_UP, Anastassiou_EMT,"Anastassiou_EMT", "FOSTER_KDM1A_TARGETS_UP" )
plot_cors(FOSTER_KDM1A_TARGETS_DN, Anastassiou_EMT,"Anastassiou_EMT", "FOSTER_KDM1A_TARGETS_DN" )
plot_cors(FOSTER_KDM1A_TARGETS_UP,  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION," HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "FOSTER_KDM1A_TARGETS_UP" )
plot_cors(FOSTER_KDM1A_TARGETS_DN,  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION," HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "FOSTER_KDM1A_TARGETS_DN" )
# EMT and LSD

# When LSD is down EMT is down

#lower correlated of when LSD is down EMT  # not a string correlation between LSD up and EMT
plot_cors(Wang_Lsd_down, Anastassiou_EMT,"Anastassiou_EMT", "WANG_LSD1_TARGETS_DN" )
plot_cors(Wang_Lsd_up, Anastassiou_EMT,"Anastassiou_EMT", "WANG_LSD1_TARGETS_UP" )

plot_cors(Wang_Lsd_up,HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "WANG_LSD1_TARGETS_UP" )
plot_cors(Wang_Lsd_down,HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "WANG_LSD1_TARGETS_down" )


plot_cors(Wang_Lsd_down, JECHLINGER_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_DN ,"JECHLINGER_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_DN", "Wang_Lsd_down" )
# no correlation between LSD down and EMT up
plot_cors(Wang_Lsd_down, JECHLINGER_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_UP ,"JECHLINGER_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_UP", "Wang_Lsd_down" )
# when LSD is up EMT is down
plot_cors(Wang_Lsd_up, JECHLINGER_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_DN ,"JECHLINGER_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_DN", "Wang_Lsd_up" )
# When lsd id up, EMT is up.
plot_cors(Wang_Lsd_up, JECHLINGER_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_UP ,"JECHLINGER_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_UP", "Wang_Lsd_up" )

dev.off()

