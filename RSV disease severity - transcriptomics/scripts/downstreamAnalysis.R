
title: "Downstream analysis of RESCEU case-control transcriptomics data"
Project: "RESCEU - case control cohort transcriptomics"
author: "Deniz Öner (ID Biomarker)"

#Create a function to install packages
install_packages <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  
  if(length(new_packages)) {
    install.packages(new_packages, dependencies = TRUE)
  } else {
    print("All packages are already installed.")
  }
}

# List of libraries to install
libraries <- c("reshape2", "tibble", "ggplot2", "superb", "reshape", 
               "GO.db", "GOstats", "ggrepel", "limma", "Biobase", 
               "esetVis", "dplyr", "sva", "tidyr", "rstudioapi")

# Call the function to install the packages
install_packages(libraries)

#required R packages
library("reshape2")
library("reshape")
library("tibble")
library("ggplot2")
library("superb")
library("GO.db")
library("GOstats")
library("ggrepel")
library("ggplot2")
library("limma")
library("Biobase")
library("esetVis")
library("dplyr")
library("sva")
library("tidyr")
library("rstudioapi")

# set the working directory
setwd(dirname(getActiveDocumentContext()$path))

#load eset data
load("../data/eset.Rdata")

#Only RSV visit and healthy
eset <- eset_comb

eset_Subset_rsvv <- eset[,eset$Visit %in% c("Healthy", "RSV_visit")]
dim(eset_Subset_rsvv)
#322

eset_Subset_rsvv <- eset_Subset_rsvv[,complete.cases(eset_Subset_rsvv$DS_Invasiveventilation)]
dim(eset_Subset_rsvv)
#268

colorPalette <- c("#D55E00", "#E69F00", "#009E73", "#56B4E9")

###
print(esetSpectralMap(eset = eset_Subset_rsvv, 
                      title = "Respiratory support",
                      colorVar = "DS_Invasiveventilation", color = colorPalette,
                      topSamples = 0, topGenes = 5, topGenesVar = "SYMBOL", cloudGenes = TRUE)) +
  theme(plot.title = element_text(size=12, lineheight = .9, family = "Verdana", colour = "black")) +
  theme(legend.position="none")

print(esetSpectralMap(eset = eset_Subset_rsvv, 
                      title = "Sex",
                      colorVar = "Sex", color = colorPalette,
                      topSamples = 0, topGenes = 5, topGenesVar = "SYMBOL", cloudGenes = TRUE)) +
  theme(plot.title = element_text(size=12, lineheight = .9, family = "Verdana", colour = "black"))

print(esetSpectralMap(eset = eset_Subset_rsvv, 
                      title = "Age",
                      colorVar = "Age_groups", color = colorPalette,
                      topSamples = 0, topGenes = 5, topGenesVar = "SYMBOL", cloudGenes = TRUE)) +
  theme(plot.title = element_text(size=12, lineheight = .9, family = "Verdana", colour = "black"))


print(esetSpectralMap(eset = eset_Subset_rsvv, 
                      title = "Resvinet severity groups",
                      colorVar = "Resvinet severity groups", color = colorPalette,
                      topSamples = 0, topGenes = 5, topGenesVar = "SYMBOL", cloudGenes = TRUE)) +
  theme(plot.title = element_text(size=12, lineheight = .9, family = "Verdana", colour = "black"))

pData$DS_resp_support <- factor(pData$DS_resp_support)

#Differential expression

esetLimma <- eset[,!is.na(eset$Sex) & !is.na(eset$baseline_age_at_visit)]

gps <- paste(pData(esetLimma)[,"Visit"],pData(esetLimma)[,"DS_Invasiveventilation"],sep=".")

myDesign <- model.matrix(~ 0 + factor(gps) + Age_groups + Sex, data = pData(esetLimma))
colnames(myDesign) <- sub("factor\\(gps\\)", "", colnames(myDesign))

# Build model
corfit <- duplicateCorrelation(esetLimma, myDesign, block = pData(esetLimma)$Subject_ID)
fit    <- lmFit(esetLimma, myDesign, correlation = corfit$consensus)

contMatrix <- makeContrasts(
  HOSPITALISATION_Outpatients.vs.Healthy.Controls =  RSV_visit.Non_hospitalized - Healthy.Healthy,
  
  HOSPITALISATION_Inpatients.vs.Healthy.Controls=     (RSV_visit.Hospitalized_wo_IV + RSV_visit.Invasive_ventilation)/2 - Healthy.Healthy,
  
  HOSPITALISATION_Ventilated.infants.vs.Outpatients =  RSV_visit.Invasive_ventilation - RSV_visit.Non_hospitalized,
  
  HOSPITALISATION_Inpatients.vs.Outpatients=     (RSV_visit.Hospitalized_wo_IV + RSV_visit.Invasive_ventilation)/2 - RSV_visit.Non_hospitalized,
  
  HOSPITALISATION_HospitilizedwithoutMV.vs.healthy = RSV_visit.Hospitalized_wo_IV - Healthy.Healthy,
  
  HOSPITALISATION_HospitilizedwitMV.vs.HospitilizedwithoutMV = RSV_visit.Invasive_ventilation - RSV_visit.Hospitalized_wo_IV,
  
  levels = myDesign)

lmFit <- eBayes(contrasts.fit(fit, contMatrix))
coefs <- colnames(lmFit) # coefficients of interest

group_means <- apply(fit$coefficients, 2, mean)
group_stddev <- apply(fit$coefficients, 2, sd)

# Outpatients vs. Healthy controls
Allgenestable <- topTable(lmFit, coef=1, number=22593)

Allgenestable <- Allgenestable %>% mutate(significant0.05 = adj.P.Val < 0.05)
Allgenestable <- Allgenestable %>% mutate(significant0.05_up = (adj.P.Val < 0.05 & logFC > 0))
Allgenestable <- Allgenestable %>% mutate(significant0.05_down = (adj.P.Val < 0.05 & logFC < 0))

P1 <- 
  ggplot(Allgenestable, aes(x=logFC , y=-log10(adj.P.Val))) + 
  geom_point(aes(color=significant0.05), size=1, alpha=0.8) + 
  scale_colour_manual(values=c("#999999", "#0072B2")) +
  #geom_vline(xintercept=1, colour= "#56B4E9", linetype="dashed", size=1) + 
  #geom_vline(xintercept = -1,  colour= "#56B4E9", linetype="dashed", size=1) +
  geom_hline(yintercept=-log10(0.05), colour= "#56B4E9", linetype="dashed", size=1) +
  #ggtitle("Outpatients vs. Healthy controls") +
  theme(plot.title = element_text(size=14, lineheight = .9, family = "Verdana", colour = "black"),
        axis.title = element_text(size=14, lineheight = .9, family = "Verdana", colour = "black"),
        axis.text = element_text(size=12, colour = "black")) +
  theme(legend.position="none") + xlab("logFC") + ylab("-log(adjusted p value)")


top <- Allgenestable %>% filter(Allgenestable$adj.P.Val < 0.05) 

P1+geom_text_repel(data=head(top, 10), aes(label=SYMBOL)) 

#GO analysis
down <- Allgenestable %>% filter(significant0.05_down == TRUE) %>% dplyr::select(ENTREZID)
selectGenesDown <- unlist(down)
up <- Allgenestable %>% filter(significant0.05_up == TRUE) %>% dplyr::select(ENTREZID)
selectGenesUp <- unlist(up)

universeGenes = unique(Allgenestable$ENTREZID)
cutOff=0.05
#BiocManager::install("org.Hs.eg.db")

upParams = new("GOHyperGParams",
               geneIds=selectGenesUp,
               universeGeneIds =universeGenes,
               annotation="org.Hs.eg.db",
               ontology="BP",
               pvalueCutoff=cutOff,
               conditional=FALSE,
               testDirection="over")

downParams = new("GOHyperGParams",
                 geneIds=selectGenesDown,
                 universeGeneIds =universeGenes,
                 annotation="org.Hs.eg.db",
                 ontology="BP",
                 pvalueCutoff=cutOff,
                 conditional=FALSE,
                 testDirection="over")


upBP = hyperGTest(upParams)
summary(upBP)[1:10,]
up_sum <- summary(upBP)[1:5,]
up_sum$Expression <- "up"

downBP = hyperGTest(downParams)
summary(downBP)[1:10,]
down_sum <- summary(downBP)[1:5,]
down_sum$Expression <- "down"

topGO = rbind(up_sum, down_sum)

ggplot(topGO, aes(x=-log(Pvalue), y=reorder(Term, Pvalue), size = Count)) +
  geom_point(aes(color=Expression), alpha=0.7) +
  scale_colour_manual(values=c("#D55E00", "#009E73")) +
  scale_size(range=c(1,10)) +
  theme_bw() +
  ggtitle("GO Terms Biological Pathways Outpatients vs. healhy infants") +
  theme(plot.title = element_text(size=14, lineheight = .9, family = "Verdana", colour = "black"),
        axis.title = element_text(size=14, lineheight = .9, family = "Verdana", colour = "black"),
        axis.text = element_text(size=12, colour = "black")) +
  xlab("-log(p value)") + ylab("GO Terms BP") # + theme(legend.position="none")

# Inpatients vs. Healthy controls
Allgenestable <- topTable(lmFit, coef=2, number=22593)

Allgenestable <- Allgenestable %>% mutate(significant0.05 = adj.P.Val < 0.05)
Allgenestable <- Allgenestable %>% mutate(significant0.05_up = (adj.P.Val < 0.05 & logFC > 0))
Allgenestable <- Allgenestable %>% mutate(significant0.05_down = (adj.P.Val < 0.05 & logFC < 0))

P1 <- 
  ggplot(Allgenestable, aes(x=logFC , y=-log10(adj.P.Val))) + 
  geom_point(aes(color=significant0.05), size=1, alpha=0.8) + 
  scale_colour_manual(values=c("#999999", "#0072B2")) +
  #geom_vline(xintercept=1, colour= "#56B4E9", linetype="dashed", size=1) + 
  #geom_vline(xintercept = -1,  colour= "#56B4E9", linetype="dashed", size=1) +
  geom_hline(yintercept=-log10(0.05), colour= "#56B4E9", linetype="dashed", size=1) +
  #ggtitle("Inpatients vs. Healthy controls") +
  theme(plot.title = element_text(size=14, lineheight = .9, family = "Verdana", colour = "black"),
        axis.title = element_text(size=14, lineheight = .9, family = "Verdana", colour = "black"),
        axis.text = element_text(size=12, colour = "black")) +
  theme(legend.position="none") + xlab("logFC") + ylab("-log(adjusted p value)")


top <- Allgenestable %>% filter(Allgenestable$adj.P.Val < 0.05)

P1+geom_text_repel(data=head(top, 10), aes(label=SYMBOL)) 

#GO analysis
down <- Allgenestable %>% filter(significant0.05_down == TRUE) %>% dplyr::select(ENTREZID)
selectGenesDown <- unlist(down)
up <- Allgenestable %>% filter(significant0.05_up == TRUE) %>% dplyr::select(ENTREZID)
selectGenesUp <- unlist(up)

universeGenes = unique(Allgenestable$ENTREZID)
cutOff=0.05

upParams = new("GOHyperGParams",
               geneIds=selectGenesUp,
               universeGeneIds =universeGenes,
               annotation="org.Hs.eg.db",
               ontology="BP",
               pvalueCutoff=cutOff,
               conditional=FALSE,
               testDirection="over")

downParams = new("GOHyperGParams",
                 geneIds=selectGenesDown,
                 universeGeneIds =universeGenes,
                 annotation="org.Hs.eg.db",
                 ontology="BP",
                 pvalueCutoff=cutOff,
                 conditional=FALSE,
                 testDirection="over")


upBP = hyperGTest(upParams)
summary(upBP)[1:10,]
up_sum <- summary(upBP)[1:10,]
up_sum$Expression <- "up"


downBP = hyperGTest(downParams)
summary(downBP)[1:10,]
down_sum <- summary(downBP)[1:10,]
down_sum$Expression <- "down"

topGO = rbind(up_sum, down_sum)

ggplot(topGO, aes(x=-log(Pvalue), y=reorder(Term, Pvalue), size = Count)) +
  geom_point(aes(color=Expression), alpha=0.7) +
  scale_colour_manual(values=c("#D55E00", "#009E73")) +
  scale_size(range=c(1,10)) +
  theme_bw() +
  ggtitle("GO Terms Biological Pathways Inpatients vs. Healthy") +
  theme(plot.title = element_blank(),
        axis.title = element_text(size=14, lineheight = .9, family = "Verdana", colour = "black"),
        axis.text = element_text(size=12, colour = "black")) +
  xlab("-log(p value)") + ylab("GO Terms BP") # + theme(legend.position="none")


#Ventilated infants vs. Outpatients
Allgenestable <- topTable(lmFit, coef=3, number=22593)

Allgenestable <- Allgenestable %>% mutate(significant0.05 = adj.P.Val < 0.05)
Allgenestable <- Allgenestable %>% mutate(significant0.05_up = (adj.P.Val < 0.05 & logFC > 0))
Allgenestable <- Allgenestable %>% mutate(significant0.05_down = (adj.P.Val < 0.05 & logFC < 0))

P1 <- 
  ggplot(Allgenestable, aes(x=logFC , y=-log10(adj.P.Val))) + 
  geom_point(aes(color=significant0.05), size=1, alpha=0.8) + 
  scale_colour_manual(values=c("#999999", "#0072B2")) +
  #geom_vline(xintercept=1, colour= "#56B4E9", linetype="dashed", size=1) + 
  #geom_vline(xintercept = -1,  colour= "#56B4E9", linetype="dashed", size=1) +
  geom_hline(yintercept=-log10(0.05), colour= "#56B4E9", linetype="dashed", size=1) +
  #ggtitle("Ventilated infants vs. Outpatients") +
  theme(plot.title = element_text(size=14, lineheight = .9, family = "Verdana", colour = "black"),
        axis.title = element_text(size=14, lineheight = .9, family = "Verdana", colour = "black"),
        axis.text = element_text(size=12, colour = "black")) +
  theme(legend.position="none") + xlab("logFC") + ylab("-log(adjusted p value)")


top <- Allgenestable %>% filter(Allgenestable$adj.P.Val < 0.05) # & (Allgenestable$logFC > 1 | Allgenestable$logFC < -1))
P1+geom_text_repel(data=head(top, 10), aes(label=SYMBOL)) 

#GO analysis
down <- Allgenestable %>% filter(significant0.05_down == TRUE) %>% dplyr::select(ENTREZID)
selectGenesDown <- unlist(down)
up <- Allgenestable %>% filter(significant0.05_up == TRUE) %>% dplyr::select(ENTREZID)
selectGenesUp <- unlist(up)

universeGenes = unique(Allgenestable$ENTREZID)
cutOff=0.05

upParams = new("GOHyperGParams",
               geneIds=selectGenesUp,
               universeGeneIds =universeGenes,
               annotation="org.Hs.eg.db",
               ontology="BP",
               pvalueCutoff=cutOff,
               conditional=FALSE,
               testDirection="over")

downParams = new("GOHyperGParams",
                 geneIds=selectGenesDown,
                 universeGeneIds =universeGenes,
                 annotation="org.Hs.eg.db",
                 ontology="BP",
                 pvalueCutoff=cutOff,
                 conditional=FALSE,
                 testDirection="over")


upBP = hyperGTest(upParams)
summary(upBP)[1:10,]
up_sum <- summary(upBP)[1:5,]
up_sum$Expression <- "up"


downBP = hyperGTest(downParams)
summary(downBP)[1:10,]
down_sum <- summary(downBP)[1:5,]
down_sum$Expression <- "down"

topGO = rbind(up_sum, down_sum)

ggplot(topGO, aes(x=-log(Pvalue), y=reorder(Term, Pvalue), size = Count)) +
  geom_point(aes(color=Expression), alpha=0.7) +
  scale_colour_manual(values=c("#D55E00", "#009E73")) +
  scale_size(range=c(1,10)) +
  theme_bw() +
  ggtitle("GO Terms Biological Pathways \nVentilated infants vs. Outpatients") +
  theme(plot.title = element_text(size=14, lineheight = .9, family = "Verdana", colour = "black"),
        axis.title = element_text(size=14, lineheight = .9, family = "Verdana", colour = "black"),
        axis.text = element_text(size=12, colour = "black")) +
  xlab("-log(p value)") + ylab("GO Terms BP") # + theme(legend.position="none")


# HospitilizedwithoutMV.vs.healthy

Allgenestable <- topTable(lmFit, coef=5, number=22593)

Allgenestable <- Allgenestable %>% mutate(significant0.05 = adj.P.Val < 0.05)
Allgenestable <- Allgenestable %>% mutate(significant0.05_up = (adj.P.Val < 0.05 & logFC > 0))
Allgenestable <- Allgenestable %>% mutate(significant0.05_down = (adj.P.Val < 0.05 & logFC < 0))

P1 <- 
  ggplot(Allgenestable, aes(x=logFC , y=-log10(adj.P.Val))) + 
  geom_point(aes(color=significant0.05), size=1, alpha=0.8) + 
  scale_colour_manual(values=c("#999999", "#0072B2")) +
  #geom_vline(xintercept=1, colour= "#56B4E9", linetype="dashed", size=1) + 
  #geom_vline(xintercept = -1,  colour= "#56B4E9", linetype="dashed", size=1) +
  geom_hline(yintercept=-log10(0.05), colour= "#56B4E9", linetype="dashed", size=1) +
  #ggtitle(" HospitilizedwithoutMV.vs.nothospitalized") +
  theme(plot.title = element_text(size=14, lineheight = .9, family = "Verdana", colour = "black"),
        axis.title = element_text(size=14, lineheight = .9, family = "Verdana", colour = "black"),
        axis.text = element_text(size=12, colour = "black")) +
  theme(legend.position="none") + xlab("logFC") + ylab("-log(adjusted p value)")


top <- Allgenestable %>% filter(Allgenestable$adj.P.Val < 0.05) # & (Allgenestable$logFC > 1 | Allgenestable$logFC < -1))
P1+geom_text_repel(data=head(top, 10), aes(label=SYMBOL)) 

#GO analysis
down <- Allgenestable %>% filter(significant0.05_down == TRUE) %>% dplyr::select(ENTREZID)
selectGenesDown <- unlist(down)
up <- Allgenestable %>% filter(significant0.05_up == TRUE) %>% dplyr::select(ENTREZID)
selectGenesUp <- unlist(up)

universeGenes = unique(Allgenestable$ENTREZID)
cutOff=0.05

upParams = new("GOHyperGParams",
               geneIds=selectGenesUp,
               universeGeneIds =universeGenes,
               annotation="org.Hs.eg.db",
               ontology="BP",
               pvalueCutoff=cutOff,
               conditional=FALSE,
               testDirection="over")

downParams = new("GOHyperGParams",
                 geneIds=selectGenesDown,
                 universeGeneIds =universeGenes,
                 annotation="org.Hs.eg.db",
                 ontology="BP",
                 pvalueCutoff=cutOff,
                 conditional=FALSE,
                 testDirection="over")


upBP = hyperGTest(upParams)
summary(upBP)[1:10,]
up_sum <- summary(upBP)[1:5,]
up_sum$Expression <- "up"


downBP = hyperGTest(downParams)
summary(downBP)[1:10,]
down_sum <- summary(downBP)[1:5,]
down_sum$Expression <- "down"

topGO = rbind(up_sum, down_sum)

ggplot(topGO, aes(x=-log(Pvalue), y=reorder(Term, Pvalue), size = Count)) +
  geom_point(aes(color=Expression), alpha=0.7) +
  scale_colour_manual(values=c("#D55E00", "#009E73")) +
  scale_size(range=c(1,10)) +
  theme_bw() +
  ggtitle("GO Terms Biological Pathways \n HospitilizedwithoutMV.vs.healthy") +
  theme(plot.title = element_text(size=14, lineheight = .9, family = "Verdana", colour = "black"),
        axis.title = element_text(size=14, lineheight = .9, family = "Verdana", colour = "black"),
        axis.text = element_text(size=12, colour = "black")) +
  xlab("-log(p value)") + ylab("GO Terms BP") # + theme(legend.position="none")


#HospitilizedwitMV.vs.HospitilizedwithoutMV

Allgenestable <- topTable(lmFit, coef=6, number=22593)

Allgenestable <- Allgenestable %>% mutate(significant0.05 = adj.P.Val < 0.05)
Allgenestable <- Allgenestable %>% mutate(significant0.05_up = (adj.P.Val < 0.05 & logFC > 0))
Allgenestable <- Allgenestable %>% mutate(significant0.05_down = (adj.P.Val < 0.05 & logFC < 0))

P1 <- 
  ggplot(Allgenestable, aes(x=logFC , y=-log10(adj.P.Val))) + 
  geom_point(aes(color=significant0.05), size=1, alpha=0.8) + 
  scale_colour_manual(values=c("#999999", "#0072B2")) +
  #geom_vline(xintercept=1, colour= "#56B4E9", linetype="dashed", size=1) + 
  #geom_vline(xintercept = -1,  colour= "#56B4E9", linetype="dashed", size=1) +
  geom_hline(yintercept=-log10(0.05), colour= "#56B4E9", linetype="dashed", size=1) +
  #ggtitle("HospitilizedwitMV.vs.HospitilizedwithoutMV") +
  theme(plot.title = element_text(size=14, lineheight = .9, family = "Verdana", colour = "black"),
        axis.title = element_text(size=14, lineheight = .9, family = "Verdana", colour = "black"),
        axis.text = element_text(size=12, colour = "black")) +
  theme(legend.position="none") + xlab("logFC") + ylab("-log(adjusted p value)")


top <- Allgenestable %>% filter(Allgenestable$adj.P.Val < 0.05) # & (Allgenestable$logFC > 1 | Allgenestable$logFC < -1))
P1+geom_text_repel(data=head(top, 10), aes(label=SYMBOL)) 

#GO analysis
down <- Allgenestable %>% filter(significant0.05_down == TRUE) %>% dplyr::select(ENTREZID)
selectGenesDown <- unlist(down)
up <- Allgenestable %>% filter(significant0.05_up == TRUE) %>% dplyr::select(ENTREZID)
selectGenesUp <- unlist(up)

universeGenes = unique(Allgenestable$ENTREZID)
cutOff=0.05

upParams = new("GOHyperGParams",
               geneIds=selectGenesUp,
               universeGeneIds =universeGenes,
               annotation="org.Hs.eg.db",
               ontology="BP",
               pvalueCutoff=cutOff,
               conditional=FALSE,
               testDirection="over")

downParams = new("GOHyperGParams",
                 geneIds=selectGenesDown,
                 universeGeneIds =universeGenes,
                 annotation="org.Hs.eg.db",
                 ontology="BP",
                 pvalueCutoff=cutOff,
                 conditional=FALSE,
                 testDirection="over")


upBP = hyperGTest(upParams)
summary(upBP)[1:10,]
up_sum <- summary(upBP)[1:5,]
up_sum$Expression <- "up"


downBP = hyperGTest(downParams)
summary(downBP)[1:10,]
down_sum <- summary(downBP)[1:5,]
down_sum$Expression <- "down"

topGO = rbind(up_sum, down_sum)

ggplot(topGO, aes(x=-log(Pvalue), y=reorder(Term, Pvalue), size = Count)) +
  geom_point(aes(color=Expression), alpha=0.7) +
  scale_colour_manual(values=c("#D55E00", "#009E73")) +
  scale_size(range=c(1,10)) +
  theme_bw() +
  ggtitle("GO Terms Biological Pathways \n HospitilizedwitMV.vs.HospitilizedwithoutMV") +
  theme(plot.title = element_text(size=14, lineheight = .9, family = "Verdana", colour = "black"),
        axis.title = element_text(size=14, lineheight = .9, family = "Verdana", colour = "black"),
        axis.text = element_text(size=12, colour = "black")) +
  xlab("-log(p value)") + ylab("GO Terms BP") # + theme(legend.position="none")

############################################################################################################################################################
##Regression analysis for the effect of age

df_eSet <- exprs(eset_Subset_rsvv)
fDataeset <- fData(eset_Subset_rsvv)
pDataeset <- pData(eset_Subset_rsvv)
annot <- fDataeset

mfile  <- merge(annot, df_eSet, by = "row.names", all = TRUE)
mfile <- mfile[,-c(1:3, 5:13)]

#genes in medium orchid module were selected
sel_genes <- c("ARG1", "LCN2", "MMP8", "OLFM4", "AZU1", "ELANE", "CEACAM8", "BPI")

sel <- mfile %>% dplyr::filter(SYMBOL %in% sel_genes)

sel[8, 1] <- "BPI2"
rownames(sel) <- sel$SYMBOL
sel$SYMBOL <- NULL 
str(sel)
sel <- as.matrix(sel)

melt <- reshape::melt(sel)
colnames(melt)[2] <- paste("sampleNames")

df <- tibble::rownames_to_column(pDataeset, "sampleNames")

melt <- merge(melt, df, by="sampleNames", all.x=TRUE)

melt <- melt %>% dplyr::filter(X1 != "BPI2")

desired_order_severity <- c("Healthy", "Non_hospitalized", "Hospitalized_wo_IV", "Invasive_ventilation")
desired_order_age <- c("Below_three_months", "Between_three_six_months", "Above_six_months")

melt$desired_order_severity <- factor(melt$DS_Invasiveventilation, levels = desired_order_severity)
melt$age_ordered <- factor(melt$Age_groups, levels = desired_order_age)

levels(melt$age_ordered) <- c("< 3m", "3-6m", "> 6m")

plot_data <- melt %>% dplyr::select(desired_order_severity, value, X1, age_ordered)


#apply anova tests
anova_model_onegene_x <- plot_data %>% 
  dplyr::select(desired_order_severity, value, X1, age_ordered) %>%
  dplyr::filter(X1=="OLFM4") %>% dplyr::filter(age_ordered == "< 3m")

# Performing the ANOVA
aov_model_x <- aov(value ~ desired_order_severity, data = anova_model_onegene_x)

# Running the post-hoc Tukey HSD test
posthoc_x <- TukeyHSD(aov_model_x, which = "desired_order_severity", conf.level=.95)

print(posthoc_x)

plot_olfm4 <- ggplot(anova_model_onegene_x, aes(x=age_ordered, y=value, fill=desired_order_severity)) +
  geom_boxplot() +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2, binwidth=2, position = position_dodge(width = 0.9)) +
  labs(y= "Log2 intensity", x= "Age in days") +
  facet_wrap(. ~ X1, ncol = 4) +
  scale_fill_manual(values=c("#FDDBC7", "#F4A582", "#D6604D", "#B2182B"), 
                    name="Hospitalisation status",
                    breaks=c("Healthy", "Non_hospitalized", "Hospitalized_wo_IV", "Invasive_ventilation"),
                    labels=c("Healthy controls", "Mild", "Moderate", "Severe")) +
  theme_minimal() +
  theme(legend.title = element_text(face = "bold")) + 
  showSignificance( c(0.72,0.9), 12.5, -0.05, "n.s.") + # Mild vs. healthy
  showSignificance( c(0.72,1.1), 13.7, -0.05, "**") + # Moderate vs. healthy
  showSignificance( c(0.72,1.28), 14.5, -0.05, "***") + # Severe vs. healthy
  showSignificance( c(0.9,1.1), 2.3, +0.05, "**") + # Moderate vs. mild
  showSignificance( c(0.9,1.28), 1.5, +0.05, "***") + # Severe vs. mild
  showSignificance( c(1.1,1.28), 4, +0.05, "n.s.") # Severe vs. moderate 

plot_olfm4


#ARG1

anova_model_onegene_x <- plot_data %>% 
  dplyr::select(desired_order_severity, value, X1, age_ordered) %>%
  dplyr::filter(X1=="ARG1") %>% dplyr::filter(age_ordered == "< 3m")


# Performing the  ANOVA
aov_model_x <- aov(value ~ desired_order_severity, data = anova_model_onegene_x)

# Running the post-hoc Tukey HSD test
posthoc_x <- TukeyHSD(aov_model_x, which = "desired_order_severity", conf.level=.95)

print(posthoc_x)

plot_arg1 <- ggplot(anova_model_onegene_x, aes(x=age_ordered, y=value, fill=desired_order_severity)) +
  geom_boxplot() +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2, binwidth=2, position = position_dodge(width = 0.9)) +
  labs(y= "Log2 intensity", x= "Age in days") +
  facet_wrap(. ~ X1, ncol = 4) +
  scale_fill_manual(values=c("#FDDBC7", "#F4A582", "#D6604D", "#B2182B"), 
                    name="Hospitalisation status",
                    breaks=c("Healthy", "Non_hospitalized", "Hospitalized_wo_IV", "Invasive_ventilation"),
                    labels=c("Healthy controls", "Mild", "Moderate", "Severe")) +
  theme_minimal() + theme(legend.position = "none") +
  theme(legend.title = element_text(face = "bold")) + 
  showSignificance( c(0.72,0.9), 10, -0.05, "n.s.") + # Mild vs. healthy
  showSignificance( c(0.72,1.1), 11, -0.05, "***") + # Moderate vs. healthy
  showSignificance( c(0.72,1.28), 12, -0.05, "***") + # Severe vs. healthy
  showSignificance( c(0.9,1.1), 5.5, +0.05, "n.s.") + # Moderate vs. mild
  showSignificance( c(0.9,1.28), 3.5, +0.05, "***") + # Severe vs. mild
  showSignificance( c(1.1,1.28), 4.5, +0.05, "***") # Severe vs. moderate 

plot_arg1

#LCN2

anova_model_onegene_x <- plot_data %>% 
  dplyr::select(desired_order_severity, value, X1, age_ordered) %>%
  dplyr::filter(X1=="LCN2") %>% dplyr::filter(age_ordered == "< 3m")


# Performing the  ANOVA
aov_model_x <- aov(value ~ desired_order_severity, data = anova_model_onegene_x)

# Running the post-hoc Tukey HSD test
posthoc_x <- TukeyHSD(aov_model_x, which = "desired_order_severity", conf.level=.95)

print(posthoc_x)

plot_lcn2 <- ggplot(anova_model_onegene_x, aes(x=age_ordered, y=value, fill=desired_order_severity)) +
  geom_boxplot() +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2, binwidth=2, position = position_dodge(width = 0.9)) +
  labs(y= "Log2 intensity", x= "Age in days") +
  facet_wrap(. ~ X1, ncol = 4) +
  scale_fill_manual(values=c("#FDDBC7", "#F4A582", "#D6604D", "#B2182B"), 
                    name="Hospitalisation status",
                    breaks=c("Healthy", "Non_hospitalized", "Hospitalized_wo_IV", "Invasive_ventilation"),
                    labels=c("Healthy controls", "Mild", "Moderate", "Severe")) +
  theme_minimal() + theme(legend.position = "none") +
  theme(legend.title = element_text(face = "bold")) + 
  showSignificance( c(0.72,0.9), 11, -0.05, "n.s.") + # Mild vs. healthy
  showSignificance( c(0.72,1.1), 12.2, -0.05, "n.s.") + # Moderate vs. healthy
  showSignificance( c(0.72,1.28), 13, -0.05, "***") + # Severe vs. healthy
  showSignificance( c(0.9,1.1), 6, +0.05, "n.s.") + # Moderate vs. mild
  showSignificance( c(0.9,1.28), 5.3, +0.05, "***") + # Severe vs. mild
  showSignificance( c(1.1,1.28), 4.5, +0.05, "***") # Severe vs. moderate 

plot_lcn2


#MMP8

anova_model_onegene_x <- plot_data %>% 
  dplyr::select(desired_order_severity, value, X1, age_ordered) %>%
  dplyr::filter(X1=="MMP8") %>% dplyr::filter(age_ordered == "< 3m")


# Performing the  ANOVA
aov_model_x <- aov(value ~ desired_order_severity, data = anova_model_onegene_x)

# Running the post-hoc Tukey HSD test
posthoc_x <- TukeyHSD(aov_model_x, which = "desired_order_severity", conf.level=.95)

print(posthoc_x)

plot_mmp8 <- ggplot(anova_model_onegene_x, aes(x=age_ordered, y=value, fill=desired_order_severity)) +
  geom_boxplot() +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2, binwidth=2, position = position_dodge(width = 0.9)) +
  labs(y= "Log2 intensity", x= "Age in days") +
  facet_wrap(. ~ X1, ncol = 4) +
  scale_fill_manual(values=c("#FDDBC7", "#F4A582", "#D6604D", "#B2182B"), 
                    name="Hospitalisation status",
                    breaks=c("Healthy", "Non_hospitalized", "Hospitalized_wo_IV", "Invasive_ventilation"),
                    labels=c("Healthy controls", "Mild", "Moderate", "Severe")) +
  theme_minimal() + theme(legend.position = "none") +
  theme(legend.title = element_text(face = "bold")) + 
  showSignificance( c(0.72,0.9), 11.5, -0.05, "n.s.") + # Mild vs. healthy
  showSignificance( c(0.72,1.1), 13.2, -0.05, "**") + # Moderate vs. healthy
  showSignificance( c(0.72,1.28), 14, -0.05, "***") + # Severe vs. healthy
  showSignificance( c(0.9,1.1), 5, +0.05, "**") + # Moderate vs. mild
  showSignificance( c(0.9,1.28), 4.3, +0.05, "***") + # Severe vs. mild
  showSignificance( c(1.1,1.28), 2.5, +0.05, "***") # Severe vs. moderate 

plot_mmp8



#AZU1

anova_model_onegene_x <- plot_data %>% 
  dplyr::select(desired_order_severity, value, X1, age_ordered) %>%
  dplyr::filter(X1=="AZU1") %>% dplyr::filter(age_ordered == "< 3m")


# Performing the  ANOVA
aov_model_x <- aov(value ~ desired_order_severity, data = anova_model_onegene_x)

# Running the post-hoc Tukey HSD test
posthoc_x <- TukeyHSD(aov_model_x, which = "desired_order_severity", conf.level=.95)

print(posthoc_x)

plot_azu1 <- ggplot(anova_model_onegene_x, aes(x=age_ordered, y=value, fill=desired_order_severity)) +
  geom_boxplot() +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2, binwidth=2, position = position_dodge(width = 0.9)) +
  labs(y= "Log2 intensity", x= "Age in days") +
  facet_wrap(. ~ X1, ncol = 4) +
  scale_fill_manual(values=c("#FDDBC7", "#F4A582", "#D6604D", "#B2182B"), 
                    name="Hospitalisation status",
                    breaks=c("Healthy", "Non_hospitalized", "Hospitalized_wo_IV", "Invasive_ventilation"),
                    labels=c("Healthy controls", "Mild", "Moderate", "Severe")) +
  theme_minimal() + theme(legend.position = "none") +
  theme(legend.title = element_text(face = "bold")) + 
  showSignificance( c(0.72,0.9), 7, -0.05, "n.s.") + # Mild vs. healthy
  showSignificance( c(0.72,1.1), 9, -0.05, "n.s.") + # Moderate vs. healthy
  showSignificance( c(0.72,1.28), 10, -0.05, "***") + # Severe vs. healthy
  showSignificance( c(0.9,1.1), 2.5, +0.05, "n.s.") + # Moderate vs. mild
  showSignificance( c(0.9,1.28), 1.6, +0.05, "***") + # Severe vs. mild
  showSignificance( c(1.1,1.28), 0.8, +0.05, "**") # Severe vs. moderate 

plot_azu1

#elane


anova_model_onegene_x <- plot_data %>% 
  dplyr::select(desired_order_severity, value, X1, age_ordered) %>%
  dplyr::filter(X1=="ELANE") %>% dplyr::filter(age_ordered == "< 3m")


# Performing the  ANOVA
aov_model_x <- aov(value ~ desired_order_severity, data = anova_model_onegene_x)

# Running the post-hoc Tukey HSD test
posthoc_x <- TukeyHSD(aov_model_x, which = "desired_order_severity", conf.level=.95)

print(posthoc_x)

plot_elane <- ggplot(anova_model_onegene_x, aes(x=age_ordered, y=value, fill=desired_order_severity)) +
  geom_boxplot() +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2, binwidth=2, position = position_dodge(width = 0.9)) +
  labs(y= "Log2 intensity", x= "Age in days") +
  facet_wrap(. ~ X1, ncol = 4) +
  scale_fill_manual(values=c("#FDDBC7", "#F4A582", "#D6604D", "#B2182B"), 
                    name="Hospitalisation status",
                    breaks=c("Healthy", "Non_hospitalized", "Hospitalized_wo_IV", "Invasive_ventilation"),
                    labels=c("Healthy controls", "Mild", "Moderate", "Severe")) +
  theme_minimal() + theme(legend.position = "none") +
  theme(legend.title = element_text(face = "bold")) + 
  showSignificance( c(0.72,0.9), 8, -0.05, "n.s.") + # Mild vs. healthy
  showSignificance( c(0.72,1.1), 11, -0.05, "*") + # Moderate vs. healthy
  showSignificance( c(0.72,1.28), 11.7, -0.05, "***") + # Severe vs. healthy
  showSignificance( c(0.9,1.1), 3.2, +0.05, "**") + # Moderate vs. mild
  showSignificance( c(0.9,1.28), 2.6, +0.05, "***") + # Severe vs. mild
  showSignificance( c(1.1,1.28), 1.8, +0.05, "**") # Severe vs. moderate 

plot_elane

#CEACAM8

anova_model_onegene_x <- plot_data %>% 
  dplyr::select(desired_order_severity, value, X1, age_ordered) %>%
  dplyr::filter(X1=="CEACAM8") %>% dplyr::filter(age_ordered == "< 3m")


# Performing the  ANOVA
aov_model_x <- aov(value ~ desired_order_severity, data = anova_model_onegene_x)

# Running the post-hoc Tukey HSD test
posthoc_x <- TukeyHSD(aov_model_x, which = "desired_order_severity", conf.level=.95)

print(posthoc_x)

plot_ceacam8 <- ggplot(anova_model_onegene_x, aes(x=age_ordered, y=value, fill=desired_order_severity)) +
  geom_boxplot() +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2, binwidth=2, position = position_dodge(width = 0.9)) +
  labs(y= "Log2 intensity", x= "Age in days") +
  facet_wrap(. ~ X1, ncol = 4) +
  scale_fill_manual(values=c("#FDDBC7", "#F4A582", "#D6604D", "#B2182B"), 
                    name="Hospitalisation status",
                    breaks=c("Healthy", "Non_hospitalized", "Hospitalized_wo_IV", "Invasive_ventilation"),
                    labels=c("Healthy controls", "Mild", "Moderate", "Severe")) +
  theme_minimal() + theme(legend.position = "none") +
  theme(legend.title = element_text(face = "bold")) + 
  showSignificance( c(0.72,0.9), 9.5, -0.05, "n.s.") + # Mild vs. healthy
  showSignificance( c(0.72,1.1), 11.7, -0.05, "n.s.") + # Moderate vs. healthy
  showSignificance( c(0.72,1.28), 12.5, -0.05, "**") + # Severe vs. healthy
  showSignificance( c(0.9,1.1), 3.5, +0.05, "*") + # Moderate vs. mild
  showSignificance( c(0.9,1.28), 2.8, +0.05, "***") + # Severe vs. mild
  showSignificance( c(1.1,1.28), 2, +0.05, "*") # Severe vs. moderate 

plot_ceacam8

#BPI

anova_model_onegene_x <- plot_data %>% 
  dplyr::select(desired_order_severity, value, X1, age_ordered) %>%
  dplyr::filter(X1=="BPI") %>% dplyr::filter(age_ordered == "< 3m")


# Performing the  ANOVA
aov_model_x <- aov(value ~ desired_order_severity, data = anova_model_onegene_x)

# Running the post-hoc Tukey HSD test
posthoc_x <- TukeyHSD(aov_model_x, which = "desired_order_severity", conf.level=.95)

print(posthoc_x)

plot_bpi <- ggplot(anova_model_onegene_x, aes(x=age_ordered, y=value, fill=desired_order_severity)) +
  geom_boxplot() +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2, binwidth=2, position = position_dodge(width = 0.9)) +
  labs(y= "Log2 intensity", x= "Age in days") +
  facet_wrap(. ~ X1, ncol = 4) +
  scale_fill_manual(values=c("#FDDBC7", "#F4A582", "#D6604D", "#B2182B"), 
                    name="Hospitalisation status",
                    breaks=c("Healthy", "Non_hospitalized", "Hospitalized_wo_IV", "Invasive_ventilation"),
                    labels=c("Healthy controls", "Mild", "Moderate", "Severe")) +
  theme_minimal() + theme(legend.position = "none") +
  theme(legend.title = element_text(face = "bold")) + 
  showSignificance( c(0.72,0.9), 10.5, -0.05, "n.s.") + # Mild vs. healthy
  showSignificance( c(0.72,1.1), 12, -0.05, "n.s.") + # Moderate vs. healthy
  showSignificance( c(0.72,1.28), 12.7, -0.05, "***") + # Severe vs. healthy
  showSignificance( c(0.9,1.1), 5.2, +0.05, "n.s.") + # Moderate vs. mild
  showSignificance( c(0.9,1.28), 4.5, +0.05, "***") + # Severe vs. mild
  showSignificance( c(1.1,1.28), 3.8, +0.05, "**") # Severe vs. moderate 

plot_bpi

## running the pls-da analysis
fdata_check = fData(eset_Subset_rsvv)
fdata_check = fdata_check[!duplicated(fdata_check$SYMBOL),]
fdata_check = fdata_check[!is.na(fdata_check$GENENAME),]


X <- t(eset_Subset_rsvv@assayData$exprs) # use the gene expression data as the X matrix
Y <- as.factor(pData$Resvinet_severity_groups)

X <- X[,rownames(fdata_check)]
colnames(X) <- fdata_check$SYMBOL

result.plsda.srbct <- plsda(X, Y,scale = T, ncomp = 2) # run the method
plotIndiv(result.plsda.srbct) # plot the samples
plotVar(result.plsda.srbct) # plot the variables

splsda.result <- splsda(X, Y, keepX = c(100,100), ncomp=2)
plotIndiv(splsda.result) # plot the samples
plotVar(splsda.result) # plot the variables

plot_df = data.frame(splsda.result$variates$X)
plot_df$group = Y

# Create a scatterplot
scatter_splsda = ggplot(plot_df, aes(x = comp1, y = comp2, color=group)) +
  geom_point() +
  labs(x = "Component 1", y = "Component 2") +
  ggtitle("Sample Projections for Components 1 and 2") + theme_bw()
scatter_splsda

loading_plot_df = data.frame(loadings_x = splsda.result$loadings$X[,1], loadings_y = splsda.result$loadings$X[,2])
loading_plot_df$genesymbol = rownames(loading_plot_df)
filtered_loading_plot <- loading_plot_df %>%
  filter(loadings_x != 0 | loadings_y != 0)

filtered_loading_plot = ggplot(filtered_loading_plot, aes(x=loadings_x, y=loadings_y, label=genesymbol)) + geom_point() + theme_bw() + geom_label_repel()
