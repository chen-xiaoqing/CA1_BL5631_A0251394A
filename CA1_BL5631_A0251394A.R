library(GEOquery)
library(oligo)
library(limma)
library(tidyverse)
library(ggplot2)
BiocManager::install('org.Hs.eg.db',force=TRUE)
library(org.Hs.eg.db)

gse50737<-getGEO("GSE50737",GSEMatrix=TRUE)
gse50737
class(gse50737)
gse50737<-gse50737[[1]]
class(gse50737)

library(dplyr)
pData(gse50737)
pData(gse50737)[['treatment:ch1']]

design<-model.matrix(~0+gse50737[['treatment:ch1']]) 
colnames(design)<-levels(as.factor(gse50737[['treatment:ch1']]))
design

#Fit in linear model
fit <- lmFit(gse50737,design)
names(pData(gse50737))
names(pData(gse50737))[34]<-'treatment'
pd<-pData(gse50737)
pd$treatment<-as.factor(pd$treatment)
pd$treatment
levels(pd$treatment)<-c("benzene_exposed","benzene_poisoning","control")
pd$treatment
colnames(design)<-levels(pd$treatment)
design

#Make contrasts matrix
contrasts_matrix<-makeContrasts(benzene_poisoning-control,
                                benzene_exposed-control,
                                benzene_poisoning-benzene_exposed,
                                levels=design)
contrasts_matrix
library(knitr)
table(contrasts_matrix)

# Do contrasts fit
gse50737_fit <- lmFit(gse50737,design)
gse50737_fit2<-contrasts.fit(gse50737_fit,contrasts = contrasts_matrix)
gse50737_fit2 <- eBayes(gse50737_fit2)
topTable(gse50737_fit2)

topTable(gse50737_fit2,coef=1)
topTable(gse50737_fit2,coef=2)
topTable(gse50737_fit2,coef=3)
summary(decideTests(gse50737_fit2,lfc=1))

##Downstream analysis of Microarray data
#Volcano plot
#benzene_poisoning - control
interesting_genes1 <- topTable(gse50737_fit2,lfc=1,coef=1,
                              p.value = 0.05,number = Inf)
volcanoplot(gse50737_fit2,coef=1)
points(interesting_genes1[['logFC']],-log10(interesting_genes1[['P.Value']]),col='red')

#benzene_exposed - control
interesting_genes2 <- topTable(gse50737_fit2,lfc=1,coef=2,
                              p.value = 0.05,number = Inf)
volcanoplot(gse50737_fit2,coef=2)
points(interesting_genes2[['logFC']],-log10(interesting_genes2[['P.Value']]),col='red')

#benzene_poisoning - benzene_exposed
interesting_genes3 <- topTable(gse50737_fit2,lfc=1,coef=3,
                               p.value = 0.05,number = Inf)
volcanoplot(gse50737_fit2,coef=3)
points(interesting_genes3[['logFC']],-log10(interesting_genes3[['P.Value']]),col='red')

#Heatmap
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(RColorBrewer)

interesting_genes <- names(which((rowSums(abs(decideTests(gse50737_fit2,lfc=1)))) > 0))
eset_of_interest <- gse50737[interesting_genes,]
gse50737[['treatment']]
treatments<-gse50737[['treatment']]
h_treatments=HeatmapAnnotation(treatments=treatments,
                               col=list(treatments=c("benzene exposed"="cornsilk",
                                                     "benzene poisoning"="darkseagreen3",
                                                     "control"="coral")),
                               annotation_name_side='left',
                               gp = gpar(col = "black"),border=TRUE)

Heatmap(exprs(eset_of_interest), name = "exprssion",
        top_annotation = h_treatments,
        row_names_gp=gpar(fontsize=2),column_names_gp=gpar(fontsize=6),
        row_split=6)

heatmap(exprs(eset_of_interest), 
        labCol=gse50737[['treatment']] ,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))
