library(tidyverse)
library(sva)
library(RUVSeq)
library(RColorBrewer)
library("factoextra")
library(FactoMineR)
require(patchwork)
require(limma)
library(peer)
require(edgeR)
library("ggsci")

##################################################################################################################################
#
# 8/24/2020 MPM Script is broken, files cannot be found. This SCRipt computed multiple normilization methods methods and plotted. 
# Plots can be found in QC/normmethods. 
##################################################################################################################################



svaBatchCor <- function(dat, mmi, mm0,n.sv=NULL){
  dat <- as.matrix(dat)
  Y <- t(dat)
  #library(sva)
  if(is.null(n.sv))   n.sv <- num.sv(dat,mmi,method="leek")
  o <- svaseq(dat,mmi,mm0,n.sv=n.sv)
  W <- o$sv
  alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
  o$corrected <- t(Y - W %*% alpha)
  return(o)
}


exprs <- readRDS('data/subread_counts_allgood.RDS')
exprs$annotation

p<- read_csv('data/phenoData.csv')



dge <- DGEList(counts=exprs$counts)
cpms = cpm(dge)
keep = rowSums(cpms>1)>=.25*dim(dge)[2]
dge <- dge[keep,]
dge <- calcNormFactors(dge)
cpms <- cpm(dge)
dim(dge)


cpms.pca <- PCA(t(cpms), graph = FALSE)


########################## PEERS #####################################
# 8/27/2020 Cannot install PEERS

  
mod <- model.matrix(~age+etiology+race+gender,p)

 model = PEER()
 PEER_setPhenoMean(model,t(cpms))
 PEER_setCovariates(model, as.matrix(mod))
 PEER_setAdd_mean(model, TRUE)
 PEER_setNk(model,10)
 PEER_getNk(model)
 PEER_update(model)
 factors = PEER_getX(model)
 weights = PEER_getW(model)
 precision = PEER_getAlpha(model)
 residuals = PEER_getResiduals(model)
 plot(precision)
 write_tsv(as.data.frame(t(factors)),'data/PEERS_factors')
 
 
 peer_res.pca <- PCA(residuals, graph = FALSE)


 
 Y <-t(cpms)
 W <- as.matrix(factors[,10-18])
 W
 alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
 peer.cor <- t(Y - W %*% alpha)
 peer_cor.pca <-PCA(t(peer.cor), graph = FALSE)
 
 plot_grid(fviz_pca_ind(peer_res.pca, habillage = p$Library_Pool,geom='point')+ggtitle('PEER Residuals'),
           fviz_pca_ind(peer_cor.pca, habillage = p$Library_Pool,geom='point')+ggtitle('PEER corrected'),
           align='v',labels='AUTO',ncol=1
 ) 
 
 
#########################################################
 
################### SVAseq ############################## 
 
 mod0 <- model.matrix(~1,data=p)
 
 num.sv(cpms,mod,method="leek")
 svseq = svaseq(cpms,mod,mod0,n.sv=10)
 unsup_sva = svaseq(cpms,mod,mod0)
 
 sva5 <- svaBatchCor(cpms,mod,mod0,n.sv=5)
 sva10 <- svaBatchCor(cpms,mod,mod0,n.sv=10)
 sva24 <- svaBatchCor(cpms,mod,mod0,n.sv=24)
 
 sva24_factors <- sva24$sv
 rownames(sva24_factors) <- p$Sample
 saveRDS(sva24_factors,file = 'MAGnet_SVA24_factors.RDS')
 
 
 sva5.pca <- PCA(t(sva5$corrected),graph=FALSE)
 sva10.pca <- PCA(t(sva10$corrected),graph=FALSE)
 sva24.pca <- PCA(t(sva24$corrected),graph=FALSE)
  ### Write 
 write_tsv(as.data.frame(t(sva10$sv)),'data/SVA10_factors')
 write_tsv(as.data.frame(t(sva24$sv)),'data/SVA24_factors')
 saveRDS(sva24$corrected,file='data/sav24_corrected.RDS')
 
 
 
 
  
 ###############################################################################################
 ###################################### Limma remove batch effects ############################
 
 batcheffects <- removeBatchEffect(cpms, batch=p$Library_Pool, design=mod,covariates = p$TIN.median.) 
 
 batcheffects.pca <- PCA(t(batcheffects),graph=FALSE)
#####################################################
 
 
#### RUVseq with empirical controls 
 
 
### USing the previouslt dfined dge 
 y <- calcNormFactors(dge, method="upperquartile")
 y <- estimateGLMCommonDisp(y, mod)
 y <- estimateGLMTagwiseDisp(y, mod)
 
 
 fit <- glmFit(y, mod)
 lrt <- glmLRT(fit, coef=2)
 
 controls=rank(lrt$table$LR) <= 500
 batch_ruv_emp <- RUVg(dge$counts, controls, k=10)
 batch_ruv_emp$normalizedCounts
 
 ruv.pca <- PCA(t(batch_ruv_emp$normalizedCounts), graph = FALSE)
 
 
 ############################# Plot all results ########################
 
 ## Do orignal Data. 

 p$etiology <- as.factor(p$etiology)
 p$Library.Pool <- as.factor(p$Library.Pool)
 
 wrap_plots(
   fviz_pca_ind(cpms.pca, habillage = p$etiology,geom='point')+ggtitle('Orig-Status'),
   fviz_pca_ind(cpms.pca, habillage = p$Library.Pool,geom='point')+ggtitle('Orig-Library'),
   ncol=1
 )
  

 ggsave('PCA_original.pdf')
 
 
 
  
 plot_grid(fviz_pca_ind(cpms.pca, habillage = p$Library.Pool,geom='point')+ggtitle('Orig'),
           fviz_pca_ind(batcheffects.pca, habillage = p$Library.Pool,geom='point')+ggtitle('Batch Correct')+theme(legend.position="none"),
           fviz_pca_ind(sva10.pca, habillage = p$Library.Pool,geom='point')+ggtitle('SVA10')+ theme(legend.position="none"),
           fviz_pca_ind(ruv.pca, habillage = p$Library.Pool,geom='point')+ggtitle('RUV')+theme(legend.position="none"),
           align='v',labels='AUTO',ncol=1
 )
 
 ggsave('PCA_mmethods_lib.pdf')
 

 plot_grid(fviz_pca_ind(cpms.pca, habillage = p$etiology,geom='point')+ggtitle('Orig'),
           fviz_pca_ind(batcheffects.pca, habillage = p$etiology,geom='point')+ggtitle('Batch Correct')+theme(legend.position="none"),
           fviz_pca_ind(sva24.pca, habillage = p$etiology,geom='point')+ggtitle('SVA10')+ theme(legend.position="none"),
           #fviz_pca_ind(peer_res.pca, habillage = p$etiology,geom='point')+ggtitle('PEER10')+theme(legend.position="none"),
           fviz_pca_ind(ruv.pca, habillage = p$etiology,geom='point')+ggtitle('RUV')+theme(legend.position="none"),
           align='v',labels='AUTO',ncol=1,axis='l'
 )
 ggsave('PCA_mmethods_etiology.pdf')
 
 
 # SVA 5,10,24
 
 plot_grid(fviz_pca_ind(cpms.pca, habillage = p$etiology,geom='point')+ggtitle('Orig'),
           fviz_pca_ind(sva5.pca, habillage = p$etiology,geom='point')+ggtitle('SVA5')+ theme(legend.position="none"),
           fviz_pca_ind(sva10.pca, habillage = p$etiology,geom='point')+ggtitle('SVA10')+theme(legend.position="none"),
           fviz_pca_ind(sva24.pca, habillage = p$etiology,geom='point')+ggtitle('SVA24')+theme(legend.position="none"),
           align='v',labels='AUTO',ncol=1,axis='l'
 )
 ggsave('PCA_sva_etiology.pdf')
 
 
 
 plot_grid(fviz_pca_ind(cpms.pca, habillage = p$Gender,geom='point')+ggtitle('Orig'),
           fviz_pca_ind(sva5.pca, habillage = p$Gender,geom='point')+ggtitle('SVA5')+ theme(legend.position="none"),
           fviz_pca_ind(sva10.pca, habillage = p$Gender,geom='point')+ggtitle('SVA10')+theme(legend.position="none"),
           fviz_pca_ind(sva24.pca, habillage = p$Gender,geom='point')+ggtitle('SVA24')+theme(legend.position="none"),
           align='v',labels='AUTO',ncol=1,axis='l'
 )
 
 
 
  
  
