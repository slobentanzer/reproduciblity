rm(list = ls())
library(WGCNA)
library(factoextra)
library(ggplot2)
library(ggpubr)
library(plyr)
library(flashClust)
options(scipen=999)
library(flashClust)

list.files()

memory_mirnas_human = data.frame(c("hsa-miR-101-3p",
                                   "hsa-miR-143-3p",
                                   "hsa-miR-374c-5p"))
load("./data/ANONMetaData.RData") ##load normalized microRNA expression  data and phenotypic data 


print("in ANON_AD_data")
my_mirnas_in_ANON_AD_data=ANONExpr[match(memory_mirnas_human[,1], rownames(ANONExpr)),]
rownames(my_mirnas_in_ANON_AD_data) = memory_mirnas_human[,1]

data=data.frame(my_mirnas_in_ANON_AD_data)
colnames(data) = colnames(my_mirnas_in_ANON_AD_data)

exprData = data
source("./auxFunc/auxFuncxx.R")
eig.val = data.frame(eigenval)
rownames(eig.val) =  colnames(data)


eig.val$diagnosis = datPheno$diagnosis
ANONData = eig.val
ANONData$factor_diagnosis = as.factor(ANONData$diagnosis)
names(ANONData) = c("normalized_eigenvalue", "diagnosis", "factor_diagnosis")

#####outliers removed due to being processed in unique center and having missing nationalities of the participants
ANONData= ANONData[-c(91,143,47,62,144),] 

#####

plot = ANONData
plot_data= ggboxplot(plot, x = "diagnosis", y = "normalized_eigenvalue", add = "dotplot",
                     fill = "diagnosis", 
                     short.panel.labs = FALSE, outlier.shape = NA,
                     ylab = "eigen-expression")

my_comparison = list( c("MCI", "control"))

plot_final = plot_data+ stat_compare_means(comparisons = my_comparison) 
plot_final+ scale_y_continuous(breaks=seq(-2,2,0.05))
dev.off()

