library(readr)
library(dplyr)
library(sva)
library(RUVSeq)
library(RColorBrewer)
library("factoextra")
library(FactoMineR)
require(cowplot)
require(tibble)
require(limma)
require(edgeR)
library("factoextra")
library(FactoMineR)
library("ggsci")
library(ExpressExtras)
library(SPIA)
library(data.table)
library(sva)
###Exports out a correct cpm matrix
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
###

################# 
#Load count data and sample data
pData <- read_csv('data/phenoData.csv') %>% arrange(sample_name)

#load('data/subread_counts_allgood.Rdata')
cts.dedup= readRDS('data/subread_counts_allgood.RDS')

exprs <- cts.dedup$counts
# colnames(exprs)=gsub("X.home.mmorley.magnet2.MAGNET_final_corrected.bamfiles.","",colnames(exprs))
# colnames(exprs)=gsub("_dedupped.bam","",colnames(exprs))
# 
### RENAME COLS to correspond to pheno file####
pData=pData[order(pData$sample_name),]
exprs=exprs[,order(colnames(exprs))]
exprs=exprs[order(rownames(exprs)),]
colnames(exprs) <- pData$sample_name

############# Get a data.frame of gene annotations #######################
genenames <- GeneAnnotate2(as.character(rownames(exprs)),organism = "human")
rownames(genenames)=genenames$ENSEMBL
# genes= genenames %>% separate(geneloc,c("chr","start","end"))
# genes$genelength=genes$end-genes$start

#################### Create a count matrix and filter CPM ###########

#dge <- DGEList(counts=exprs)
dge <- DGEList(counts=exprs[rownames(exprs) %in% genenames$ENSEMBL,], genes=genenames)

#dge <- DGEList(counts=exprs[rownames(exprs) %in% genenames$ENSEMBL,], genes=genenames)
dim(dge)
cpms = cpm(dge)

###########-------------------------------------------#################
keep1 = rowSums(cpms>1)>=.25*dim(dge)[2]
keep2 = rowSums(cpms>.5)>=.25*dim(dge)[2]
keep3 = rowMeans(cpms) >=.5
keep4 = rowSums(cpms>0) >=.8*dim(dge)[2]
sum(keep)
sum(keep2)
sum(keep3)
sum(keep4)
###############Calculate FPKM ########################
len= fread("~/NGSshare/hg19_data/GeneLengths.txt")
genes= left_join(genenames, len, by="ENSEMBL")
y <- DGEList(counts=exprs,genes=genes) 
y=y[keep4,]
y <- calcNormFactors(y) 
FPKM <- rpkm(y)

############CREATE DGE OBJECT ##################
# In DGE object, keep only the genes that have > 0 cpm in 80% of samples
dge <- dge[keep4,]
genenames = genenames[keep4,]
dge <- calcNormFactors(dge)
cpms <- cpm(dge)
dim(dge)

################### SVAseq ############################## 
mod <- model.matrix(~age+etiology+race+gender,pData)
mod0 <- model.matrix(~1,data=pData)

## write out a covar file ###
cpms=cpms+1
unsup_sva = svaseq(as.matrix(cpms),mod,mod0)
cpms.corr <- svaBatchCor(cpms,mod,mod0,n.sv=24)$corrected
saveRDS(cpms.corr,file="data/CPMS_SVA_corrected.RDS")

unsup_sva.fpkm = svaseq(as.matrix(FPKM),mod,mod0)
fpkm.corr <- svaBatchCor(FPKM,mod,mod0,n.sv=24)$corrected
saveRDS(fpkm.corr,file="data/FPKM_SVA_corrected.RDS")


#### Add the SVs to the pData table and export them to be used later. 
#
sv.data <- unsup_sva$sv 
colnames(sv.data) = paste0('SVA',1:dim(unsup_sva$sv)[2])


pdata.sva <- cbind(pData,sv.data)

write_csv(pdata.sva,'data/pdata.sva.csv')


#### Need to make a new variable for RACE and Etiology
pData <- pData %>% mutate(disease_race = paste0(race,'_',etiology))



design = model.matrix(~0+pData$age+pData$disease_race+pData$gender + unsup_sva$sv)
design
(colnames(design) <- gsub('pData[:$:]','',colnames(design)))
(colnames(design) <- gsub('unsup_sva\\$','',colnames(design)))

contr.matrix <- makeContrasts(
  AADCMvsNF = disease_raceAA_DCM-disease_raceAA_NF,
  CauDCMvsNF = disease_raceCaucasian_DCM-disease_raceCaucasian_NF,
  CauHCMvsNF = disease_raceCaucasian_HCM-disease_raceCaucasian_NF,
  CauNFvsAANF=disease_raceCaucasian_NF-disease_raceAA_NF,
  CauDCMvsAADCM=disease_raceCaucasian_DCM-disease_raceAA_DCM,
  levels = colnames(design))
contr.matrix



v <- voom(dge,design,plot=TRUE)
pData$minexpr=abs(min(v$E))+1
v$E <- v$E + abs(min(v$E))+1
#### Save voom counts as a csv file for use later 

saveRDS(v$E,file = 'data/Voom_Matrix.RDS')

############### Create the eset #######################
## Create phenoData
rownames(pData) <- (pData$sample_name)
if(all(rownames(pData)!=colnames(v$E))){
  stop('Sample names not equal')
}
phenoData <- new("AnnotatedDataFrame",
                 data=pData)
fData <- new("AnnotatedDataFrame",
             data=genenames)
all(rownames(genenames)==rownames(v$E))



######fit models#######
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean−variance trend")

######## load and prepare all the MSigDB sets for camera ######

data('human_H_v5',package="ExpressExtras")
h.indices <- ids2indices(Hs.H,genenames$ENTREZID)
data('human_c2_v5',package="ExpressExtras")
c2.indices <- ids2indices(Hs.c2,genenames$ENTREZID)
data('human_c3_v5',package="ExpressExtras")
c3.indices <- ids2indices(Hs.c3,genenames$ENTREZID)
data('human_c5_v5',package="ExpressExtras")
GO.indices <- ids2indices(Hs.c4,genenames$ENTREZID)

##################################################################
#Loop over all contrasts and run limma, camera, topgo and spia 
(contrastnames <-gsub('-','_vs_',colnames(contr.matrix)))
#Create list to hold the results for limma,togo and camera for all contrasts
limma <-  vector(mode="list", length=length(contrastnames))
names(limma) <- contrastnames
camera <- vector(mode="list", length=length(contrastnames))
names(camera) <- contrastnames
topgo <-vector(mode="list", length=length(contrastnames))
names(topgo) <- contrastnames
spia<-  vector(mode="list", length=length(contrastnames))
names(spia) <- contrastnames

# Cleanup <-function(tt){
#   tt$ENSEMBL=rownames(tt)
#   tt=left_join(tt,genenames,by="ENSEMBL")
#   res <- tt %>% mutate(fc = ifelse(logFC<0, -1*2^abs(logFC),2^logFC)) %>%
#     dplyr::select(ENSEMBL,SYMBOL,ENTREZID,biotype,geneloc,logFC,fc,P.Value,adj.P.Val,t)
#   rownames(res) <- res$ENSEMBL
#   return(res)
# }

for(i in 1:length(contrastnames)){
  print(contrastnames[i])
  limma[[contrastnames[i]]] <- Cleanup2(topTable(efit,coef=i,n=Inf,p.value=1))
  topgo[[contrastnames[i]]] <- runTopGO(limma[[contrastnames[i]]],organism ="human")

  k=limma[[contrastnames[i]]]

  #for each limma data (corresponding to the contrast), run SPIA
  limma_sel <- k[which(abs(k$fc) > 2 & k$adj.P.Val < 0.05),]
  if(nrow(limma_sel)>0){
    all_genes = as.numeric(k$ENTREZID)
    sig_genes = limma_sel$fc
    names(sig_genes) = limma_sel$ENTREZID
    sig_genes = sig_genes[complete.cases(names(sig_genes))]
    sig_genes = sig_genes[unique(names(sig_genes))]
    spia[[contrastnames[i]]] <- spia(de=sig_genes, all=all_genes, organism='hsa')
  }else{
    spia[[contrastnames[i]]] <- data.frame()
  }
  
  #run camera
  res.h <- camera(v, h.indices, design,contr.matrix[,i],inter.gene.cor=0.01)
  res.c2 <- camera(v, c2.indices, design,contr.matrix[,i],inter.gene.cor=0.01)
  res.GO <- camera(v, GO.indices, design,contr.matrix[,i],inter.gene.cor=0.01)
  #res.c3 <- camera(v, c3.indices, design,i,inter.gene.cor=0.01)
  #res.c4 <- camera(v, c4.indices, design,i,inter.gene.cor=0.01)
  camera[[contrastnames[i]]] <- list(Hallmark=list(camera_result=res.h,indices=h.indices),Curated=list(camera_result=res.c2,indices=c2.indices),GO=list(camera_result=res.GO,indices=GO.indices))
}

AA_NF_DCM=read.csv("../eQTL/QTLtools/results/Race_stratified_perm_cond/AA_NF_DCM_conditional_annotated.csv")
Cau_NF_DCM=read.csv("../eQTL/QTLtools/results/Race_stratified_perm_cond/Cau_NF_DCM_conditional_annotated.csv")
eQTL=list(AA_NF_DCM=AA_NF_DCM,Cau_NF_DCM=Cau_NF_DCM)
############ Save list of results ##############

eset<- ExpressionSet(assayData=as.matrix(v$E),phenoData=phenoData,featureData=fData,annotation="hg19")

results <- list(eset=eset,limma=limma,camera=camera, topgo=topgo,spia=spia,eQTL=eQTL,fpkm=fpkm.corr)
save(results,file="data/magnet_final.RData")

############################

tt.CauNFvsAANF <- topTable(efit,coef = 'CauNFvsAANF',number = Inf)

tt.CauNFvsAANF[id,]

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)
  
nfvsdcm <- topTreat(tfit, coef=1, n=Inf)

cpms.corr['ENSG00000259207',] 

####gene plot ###
id=tt.CauNFvsAANF$ENSEMBL[1]
 as.data.frame(cpms.corr[id,]) %>% 
  {colnames(.)[1] = "signal"; .}%>%
   rownames_to_column(var='Sample') %>% 
   inner_join(.,p) %>%
   filter(CHF_Etiology %in% c('DCM','NF')) %>% 
   mutate(CHF_Etiology=factor(CHF_Etiology,levels=c('NF','DCM'))) %>%
   ggplot(aes(x=CHF_Etiology,color=Race,y=signal)) +  geom_point(position=position_jitterdodge(dodge.width=0.9)) + 
   geom_boxplot(alpha=0,outlier.colour = NA, 
                 position = position_dodge(width=0.9)) + xlab('Status') + ylab('CPM')

ggsave('~/dsdata/investigators1/kiran/MME_boxplot.png') 
 
 as.data.frame(cpms[id,]) %>% 
 {colnames(.)[1] = "signal"; .}%>%
   rownames_to_column(var='Sample') %>% 
   inner_join(.,p) %>%
   ggplot(aes(x=CHF_Etiology,color=Race,y=signal)) +geom_boxplot()








