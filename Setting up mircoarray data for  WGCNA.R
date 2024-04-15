library(cancerclass)
library(golubEsets)
library(dplyr)
library(cluster)
library(hu6800.db)
library(tidyverse)
library(WGCNA)
#reading in and organizing data
data(Golub_Merge)
Gene_df<-data.frame(exprs(Golub_Merge))


#Check For Outliers
IAC=cor(Gene_df,use="p")
hist(IAC,sub=paste("Mean=",format(mean(IAC[upper.tri(IAC)]),digits=3)))
cluster1=hclust(as.dist(1-IAC))
#filter your genes – usually try cutting your data in a half or third, depending on what

keepGenesExpr = rank(-rowMeans(Gene_df))<=1000
mydata<-Gene_df[keepGenesExpr,]
#we need to transpose our matrix for the format WGCNA expects
dataExpr<-t(mydata)
#Create Network these defaults will work for many situations, but you can try #changing
#the minimum module size if you get many, small modules. You can also #change the
#split level to anything between 1 and 4. Merged cut height determines #where in the
#dendrogram you cut a branch into a module
net = blockwiseModules(dataExpr, power = 7,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold =10, mergeCutHeight = 0.10,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase="TOM", verbose=3, ds=3)

mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


genes=colnames(dataExpr)
moduleColors=labels2colors(net$colors)
mymodules=cbind(genes,moduleColors)
moduleColors=labels2colors(net$colors)
mymodules=cbind(genes,moduleColors)
save

df3<-data.frame(Golub_Merge$ALL.AML)

Main_df<-data.frame(exprs(Golub_Merge))
Main_df<-t(Main_df)
Main_df<-data.frame((Main_df))
Main_df$cancer_type<-df3$Golub_Merge.ALL.AML
Genes<-Main_df %>% dplyr::select(-'cancer_type')
GT<-Genes
GT<-data.frame(GT)


GT$cancertype<- as.factor(Main_df$cancer_type)

GT$number<-c(1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,0
             ,0
             ,0
             ,0
             ,0
             ,0
             ,0
             ,0
             ,0
             ,0
             ,0
             ,0
             ,0
             ,0
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,1
             ,0
             ,0
             ,0
             ,0
             ,0
             ,0
             ,0
             ,0
             ,0
             ,0
             ,0
)


phenotype<-as.factor(GT$number) #if you have 12 arrays at timepoint 0, 1 and 7 your list would be:
#(0,0,0,0,1,1,1,7,7,7,7)
#get Eigengenes
MEs0=moduleEigengenes(dataExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, unlist(phenotype), use = "p")
moduleTraitPvalue=corPvalueStudent(moduleTraitCor, 11)


Genes<-t(Genes)
gene.names<- row.names(Genes)

library("AnnotationDbi")


PROBES<- as.character(genes)
OUT_ <- na.omit(AnnotationDbi::select(hu6800.db,keys= PROBES, columns=c("SYMBOL", "ENTREZID", "GENENAME"),keytype="PROBEID"))




exportNetworkToCytoscape(
  TOM,
  edgeFile = "C:/Users/ccape/OneDrive/Documents/edgea.tsv",
  nodeFile = "C:/Users/ccape/OneDrive/Documents/nodea.tsv",
  weighted = TRUE,
  threshold = 0.25,
  nodeNames =OUT_$SYMBOL ,
  altNodeNames = NULL,
  nodeAttr = NULL,
  includeColNames = TRUE)




