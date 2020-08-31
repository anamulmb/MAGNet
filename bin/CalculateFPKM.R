library(readr)
library(dplyr)
library(tidyr)
library(data.table)
library(ExpressExtras)
library(GenomicFeatures)
library(plyr)
library(edgeR)

#Load count data
load('data/subread_counts_allgood.Rdata')
exprs <- cts.dedup$counts
colnames(exprs)=gsub("X.home.mmorley.magnet2.MAGNET_final_corrected.bamfiles.","",colnames(exprs))
colnames(exprs)=gsub("_dedupped.bam","",colnames(exprs))

############# Get a data.frame of gene annotations #######################
genenames <- GeneAnnotate(as.character(rownames(exprs)),organism = "human")
genes= genenames %>% separate(geneloc,c("chr","start","end"))

############# Get genelengths #######################
gtffile="~/Desktop/Apoorva/NGSshare/hg19_data/Homo_sapiens.GRCh37.75.gtf"
txdb <- makeTxDbFromGFF(gtffile,format="gtf")
# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
df <- ldply (exonic.gene.sizes,data.frame)
colnames(df)=c("ENSEMBL","genelength")
genes=inner_join(genes,df,by="ENSEMBL")

#################### Create a count matrix and filter CPM ###########
dge <- DGEList(counts=exprs)
dim(dge)
cpms = cpm(dge)

###########-------------------------------------------#################
keep1 = rowSums(cpms>1)>=.25*dim(dge)[2]
keep2 = rowSums(cpms>.5)>=.25*dim(dge)[2]
keep3 = rowMeans(cpms) >=.5
keep4 = rowSums(cpms>0) >=.8*dim(dge)[2]
sum(keep1)
sum(keep2)
sum(keep3)
sum(keep4)

###############Calculate FPKM ########################
y <- DGEList(counts=exprs,genes=genes) 
y <- calcNormFactors(y) 
FPKM <- rpkm(y)
write.csv(FPKM,file="data/FPKM_allgenes_corrExonLen.csv")

y <- DGEList(counts=exprs,genes=genes) 
y=y[keep4,]
y <- calcNormFactors(y) 
write.csv(FPKM,file="data/FPKM_filteredgenes_corrExonLen.csv")


write.table(df,file="~/NGSshare/hg19_data/GeneLengths.txt",row.names = F,quote = F,sep="\t")

