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

# what is loaded here?
# what are miRNAs? how many are there in human? how many are there in this dataset?
# is the count distribution within normal parameters?


print("in ANON_AD_data")
my_mirnas_in_ANON_AD_data=ANONExpr[match(memory_mirnas_human[,1], rownames(ANONExpr)),]
rownames(my_mirnas_in_ANON_AD_data) = memory_mirnas_human[,1]

data=data.frame(my_mirnas_in_ANON_AD_data)
colnames(data) = colnames(my_mirnas_in_ANON_AD_data)

exprData = data
source("./auxFunc/auxFuncxx.R")
# what is done here? does it make sense?

eig.val = data.frame(eigenval)
rownames(eig.val) =  colnames(data)


eig.val$diagnosis = datPheno$diagnosis
ANONData = eig.val
ANONData$factor_diagnosis = as.factor(ANONData$diagnosis)
names(ANONData) = c("normalized_eigenvalue", "diagnosis", "factor_diagnosis")

#####outliers removed due to being processed in unique center and having missing nationalities of the participants
ANONData= ANONData[-c(91,143,47,62,144),] 

# is there a sufficient explanation for the removal of these samples?
# how are the eigenvalues distributed?
# what is the effect of the exclusion of these particular samples, and why?

#####

plot = ANONData
# pdf("./results/Fig_5B.pdf", useDingbats = FALSE)
plot_data= ggboxplot(plot, x = "diagnosis", y = "normalized_eigenvalue", add = "dotplot",
                     fill = "diagnosis", 
                     short.panel.labs = FALSE, outlier.shape = NA,
                     ylab = "eigen-expression")

my_comparison = list( c("MCI", "control"))

plot_final = plot_data+ stat_compare_means(comparisons = my_comparison) 
plot_final+ scale_y_continuous(breaks=seq(-2,2,0.05))

# where does the p-value that is given in the plot come from?
# which statistical procedure underlies the hypothesis test?
# is the statistical procedure adequate for the question that is asked?
# is it performed correctly with respect to statistical principles?

# print(plot_final)
dev.off()

