title: "WGCNA - Cell population deconvolution"
Project: "RESCEU - case control cohort transcriptomics 01"
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
libraries <- c("WGCNA", "Biobase", "dplyr", "reshape2", "tidyverse", "rstudioapi")

install_packages(libraries)

#nstall libraries
library("WGCNA")
library("Biobase")
library("dplyr")
library("reshape2")
library("tidyverse")
library("rstudioapi")

setwd(dirname(getActiveDocumentContext()$path))

#load exprs data
datExpr <- read.csv(".../data/datExpr.txt", sep="")
#load traits data
datTraits <- read.csv(".../data/datTraits.txt", sep="")

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
enableWGCNAThreads()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.85;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8, col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


allowWGCNAThreads()
#about 30 sec!
softPower = 3;

#We now calculate the adjacencies, using the soft thresholding power:
#we are using a signed matrix, here is some literature on the matter https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/ 
#a signed matrix will give you genes that have positive correlations, if you use unsigned you might end up with genes with pos and neg corretation on the same module
#this becomes confusing to interpret, but can be useful in certain cases.
adjacency = adjacency(datExpr, power = softPower, type = "unsigned");

#To minimize effects of noise and spurious associations, we transform the adjacency into Topological Overlap Matrix, 
#and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType = "unsigned", verbose=5);
dissTOM = 1-TOM
dim(TOM)
#save(TOM, file = "TOMfile_allsamples.txt")

###
#TOM <-read.table("TOMfile_allsamples.txt") 

geneTree = hclust(as.dist(dissTOM), method = "average");

# We like bigger modules, so we set the minimum module size relatively high:
minModuleSize = 25;

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            method = "hybrid", pamStage = TRUE, #PAM allows for more genes to be assigned to clusters
                            deepSplit = 2, maxDistToLabel = 0, #deepslit controls how fine the clusters should be split
                            minClusterSize = minModuleSize);

#dynamicMods1 = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
table(dynamicMods)

#write.table(dynamicMods, file="dynamicMods_10.txt", sep="\t")

#######################
# Convert numeric lables into colors
#dynamicMods <- read.table("dynamicMods_10.txt") 

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
#write.table(dynamicColors, file = "number-genes-per-module.txt")
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


# Calculate eigengenes
MEList = moduleEigengenes(cov, colors = dynamicColors)
MEs = MEList$eigengenes

# Save module colors and labels for use in subsequent parts
#save(MEs, dynamicMods, dynamicColors, geneTree, file = "CC_MEs.RData")

############################################################
samples = rownames(datExpr);

traitRows = match(samples, datTraits$Aros.code);
traitRows <- traitRows[!is.na(traitRows)]


datTraits1 = datTraits_gh[traitRows, -1];
rownames(datTraits1) = datTraits_gh[traitRows, 1];

dim(datTraits1)
str(datTraits1)
names(datTraits1)

#Converting dattraits to numeric values
#convert dat Traits to numeric values

datTraits_num <- datTraits1[, -c(1:2)]  

#datTraits_num <- as.numeric(datTraits_num)
str(datTraits_num)

###############################3
datExpr_cut <- subset(datExpr, row.names(datExpr) %in% row.names(datTraits_num))

collectGarbage();

# Re-cluster samples

sampleTree2 = hclust(dist(datExpr_cut), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry

traitColors = numbers2colors(datTraits_num, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.

plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits_num), 
                    main = "Sample dendrogram and trait heatmap")


#save(datExpr_cut, datTraits, file = "WGCNA_casecontrol_dataInput.RData")


##
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr_cut, dynamicColors)$eigengenes
MEs = orderMEs(MEs0)

moduleTraitCor = cor(MEs, datTraits_num, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#write.table(moduleTraitCor, file="moduleTraitCor.txt")
# moduleTraitCor <-read.table("moduleTraitCor.txt") 
# 
#write.table(datTraits_num, file="datTraits_num.txt")
# datTraits_num <-read.table("datTraits_num.txt") 
# 
#write.table(MEs, file="MEs.txt")
# MEs <-read.table("MEs.txt") 


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), " (",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits_num),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               cex.lab.x = 1,
               xLabelsAngle = 20,
               zlim = c(-1,1))


#####################################################################################
#Correlation with cell counts (cluster-based analysis)

#datExpr_cut <-read.table(".../data/datExpr_cut.txt") 
datTraits_counts <-read.table(".../data/datTraits_counts.txt") 
datTraits_counts$SUBJECT.VISIT <- NULL

nSamples = nrow(datExpr_cut)

#dynamicMods <- read.table("dynamicMods_10.txt")
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Calculate eigengenes
MEList = moduleEigengenes(datExpr_cut, colors = dynamicColors)
moduleEigengenes <- MEList$eigengenes

MEs = MEList$eigengenes

#use data_wide
moduleTraitCor = cor(MEs, datTraits_counts, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#write.table(moduleTraitCor, file="moduleTraitCor_counts.txt")

#moduleTraitCor

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), " (",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits_counts),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = FALSE,
               setStdMargins = FALSE,
               cex.text = 0.8,
               cex.lab.x = 1,
               xLabelsAngle = 20,
               zlim = c(-1,1))
