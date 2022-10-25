rm(list = ls())
library(WGCNA)
library(factoextra)
library(ggplot2)
library(ggpubr)
library(plyr)
library(flashClust)
options(scipen=999)
library(flashClust)
library(patchwork)

memory_mirnas_human = data.frame(c("hsa-miR-101-3p",
                                   "hsa-miR-143-3p",
                                   "hsa-miR-374c-5p"))
load("./data/ANONMetaData.RData") ##load normalized microRNA expression  data and phenotypic data 

# closer look at expression data
ANONExpr[1:10, 1:10]
library(vsn)
meanSdPlot(as.matrix(ANONExpr))

print("in ANON_AD_data")
my_mirnas_in_ANON_AD_data=ANONExpr[match(memory_mirnas_human[,1], rownames(ANONExpr)),]
rownames(my_mirnas_in_ANON_AD_data) = memory_mirnas_human[,1]

data=data.frame(my_mirnas_in_ANON_AD_data)
colnames(data) = colnames(my_mirnas_in_ANON_AD_data)

dim(data) # there are 2588 mature miRNAs in miRBase v21
# only 145 in here...
exprData = data

# get the code from the unnecessary aux function script
# source("./auxFunc/auxFuncxx.R")
#!/usr/bin/R
###this script calculates a decomposed value from large gene expression datasets. 
##the script is modified from WGCNA source scripts. 
###### Alternatively a decomposed value can be calculated using modduleEigengenes function from WGCNA module. 
##in both cases the output should be same. 


expr = t(exprData)

maxVarExplained = 3
nPC = 1
nVarExplained = min (nPC, maxVarExplained)
colors = "blue"
modLevels = levels(factor(colors))
PrinComps = data.frame(matrix(NA, nrow = dim(expr)[[1]], ncol = length(modLevels)))
averExpr = data.frame(matrix(NA, nrow = dim(expr)[[1]], ncol = length(modLevels)))
varExpl = data.frame(matrix(NA, nrow = dim(expr)[[1]], ncol = length(modLevels)))

modulename = modLevels[1]
restrict1 = as.character(colors) == as.character(modulename)
datModule1 = as.matrix(t(expr[, restrict1]))
n = dim(datModule1)[1]
p = dim(datModule1)[2]

###scaling
datModule = t(scale(t(datModule1)))
svd1 = svd(datModule, nu = min(n, p, nPC), nv = min (n, p, nPC))
varExpl[,1] = (svd1$d[1:min(n,p, nVarExplained)])^2/sum(svd1$d^2)
veMat = cor(svd1$v[, c(1:min(n, p, nVarExplained))], t(datModule), use = "p")


##first principle component 
eigenval = svd1$v[,1]

scaledExpr = scale(t(datModule1))
averExpr[,1] = rowMeans(scaledExpr, na.rm = TRUE)

corAve = cor(averExpr[,1], eigenval, use = "p")
corTest = cor.test(averExpr[,1], eigenval, use = "p")

if(!is.finite(corAve)){
  corAve = 0;
}
if (corAve<0){
  eigenval = -1* eigenval
}

# continuation of original script

eig.fun <- moduleEigengenes(t(data), colors = rep("blue", nrow(data)))

eig.val = data.frame(eigenval)
rownames(eig.val) =  colnames(data)

eig.fun$eigengenes == eig.val

# compute arbitrary collection of three-miRNA signatures as permutation test
p <- list()
for (i in 1:12) {
  set.seed(i) # make it reproducible
  # randomly select three miRNAs from the set
  data <- ANONExpr[sample(1:nrow(ANONExpr), 3), ]
  message(rownames(data))
  
  #do the original thing
  eig.fun <- moduleEigengenes(t(data), colors = rep("blue", nrow(data)))
  eig.val <- eig.fun$eigengenes
  eig.val$diagnosis = datPheno$diagnosis
  ANONData = eig.val
  ANONData$factor_diagnosis = as.factor(ANONData$diagnosis)
  names(ANONData) = c("normalized_eigenvalue", "diagnosis", "factor_diagnosis")
  
  #####outliers removed due to being processed in unique center and having missing nationalities of the participants
  ##cited  in line 779-780
  ANONData= ANONData[-c(91,143,47,62,144),] 
  plot = ANONData
  # pdf("./results/Fig_5B.pdf", useDingbats = FALSE)
  plot_data = ggboxplot(plot, x = "diagnosis", y = "normalized_eigenvalue", add = "dotplot",
                       fill = "diagnosis", 
                       short.panel.labs = FALSE, outlier.shape = NA,
                       ylab = "eigen-expression", title = paste(rownames(data), collapse = ", "))
  
  my_comparison = list( c("MCI", "control"))
  
  plot_final = plot_data + stat_compare_means(comparisons = my_comparison) + 
    # what are they actually doing here to compute the p-value?
    scale_y_continuous(breaks=seq(-2,2,0.05))
  p[[i]] <- plot_final
}

wrap_plots(p, nrow = 4)
ggsave("results/permutation_plots.png", width = 20, height = 25)
# how many of these "signatures" selected completely at random are similarly 
# significant as their signature?

#####
