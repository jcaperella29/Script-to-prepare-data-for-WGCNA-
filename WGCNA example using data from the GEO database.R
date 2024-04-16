library(GEOquery)

library(dplyr)
library(cluster)
library(hu6800.db)
library(tidyverse)
library(WGCNA)
# load series and platform data from GEO

gset <- getGEO("GSE205129", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL23126", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)

Gene_df<-data.frame(ex)

IAC=cor(Gene_df,use="p")
hist(IAC,sub=paste("Mean=",format(mean(IAC[upper.tri(IAC)]),digits=3)))
cluster1=hclust(as.dist(1-IAC)
                
keepGenesExpr = rank(-rowMeans(Gene_df))<=1000
mydata<-Gene_df[keepGenesExpr,]
keepGenesExpr = rank(-rowMeans(Gene_df))<=1000
mydata<-Gene_df[keepGenesExpr,]
dataExpr<-t(mydata)

net = blockwiseModules(dataExpr, power = 7,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold =10, mergeCutHeight = 0.10,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase="TOM", verbose=3, ds=3)
genes=colnames(dataExpr)
moduleColors=labels2colors(net$colors)
mymodules=cbind(genes,moduleColors)
moduleColors=labels2colors(net$colors)
mymodules=cbind(genes,moduleColors)
save
pdat<-phenoData(gset)
pdat<-as.list(pdat@data$characteristics_ch1)
phenoinfo <- lapply(pdat, function(x) replace(x, x=="0", "subject status : ADHD patient"))
phenonumber<-list()
phenonumber[1:3]<-"0"
phenonumber[4:6]<-"1"                
phenotype<-as.factor(unlist(phenonumber)) 
#get Eigengenes
MEs0=moduleEigengenes(dataExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, unlist(phenotype), use = "p")
moduleTraitPvalue=corPvalueStudent(moduleTraitCor, 11)
exportNetworkToCytoscape(
  TOM,
  edgeFile = "C:/Users/ccape/OneDrive/Documents/edgea.tsv",
  nodeFile = "C:/Users/ccape/OneDrive/Documents/nodea.tsv",
  weighted = TRUE,
  threshold = 0.25,
  nodeNames =genes ,
  altNodeNames = NULL,
  nodeAttr = NULL,
  includeColNames = TRUE)

