library(limma)
library(tidyverse)

####################################################################################################
# 8/31/2020 MPM Create this script as a quick exmaple on how to subset and run diff expression using LIMMA Voom method. 
#
#
#
####################################################################################################


#Read in the MAGNet expressionSet
dge <- readRDS(file = 'data/MAGNET_eset.RDS')


## Show Phenodata data in the eSet object
pData 


### sub select the samples you want.. This example is AA DCM vs NF.. 
sub.dge <- dge[,dge$race=='AA' & dge$etiology %in% c('NF','DCM')]
sub.dge

#Factors can be annoying in R and when you sub select, the other levels of the factor are still there, HCM, etc.. So I just refactor
sub.dge$etiology<- factor(sub.dge$etiology)

########## Define model and run voom ############
formula <- as.formula(paste0('~0+etiology+age+gender+',paste0("SVA",rep(1:24),collapse='+')))
(design <-model.matrix(formula,data=pData(sub.dge)))
#Clean up the colnames of the design matrix, don't need the colname in. 
colnames(design)<-gsub('etiology','',colnames(design))
v <- voom(sub.dge,design)


### Fit models ####
(contrast.matrix <- makeContrasts(NF-DCM,levels=design))
fit <- lmFit(v,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#Return results 
topTable(fit2,coef=1,n=25,p.value=0.05)
