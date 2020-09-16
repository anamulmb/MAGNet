![logo](https://github.com/mpmorley/MAGNet/blob/master/MAGnet_logo_heart.png)

MAGNet is a collaborative group of investigators who use genomic approaches to understand human myocardial disease. Current projects include human myocardial expression, quantitative trait mapping and systems genetics.

NHLBI Grant number 1R01HL105993-01A1


## RNAseq
The Myocardial Applied Genomics Network [MAGNet](www.med.upenn.edu/magnet), collects and banks human cardiac tissue for genomic research.  All subjects or next of kin provided written informed consent for tissue donation and analyses and all study protocols were approved by relevant institutional review boards. Left ventricular free-wall tissue was harvested at the time of cardiac surgery from subjects with heart failure undergoing transplantation and from unused donor hearts with apparently normal function. The heart was perfused with cold cardioplegia prior to cardiectomy to arrest contraction and prevent ischemic damage, and tissue specimens were frozen in liquid nitrogen. Total RNA was extracted using the miRNeasy Kit (Qiagen) including DNAse treatment. RNA concentration and quality was determined using the NanoVue Plus™ spectrophotometer (GE Healthcare) and the Agilent 2100 RNA Nano Chip (Agilent).	

RNA sequencing libraries were prepared using the Illumina TruSeq stranded mRNA kit followed by the Nugen Ovation amplification kit. To avoid confounding by batch effects, libraries were randomly selected into pools of 32, and pools were sequenced on a Hiseq2500 to a depth of ~30 million 100-bp paired-end reads per biological sample. Fastq files were aligned against human reference (hg19/hGRC37) using the [STAR aligner](https://github.com/alexdobin/STAR). Duplicate reads were removed using MarkDuplicates from Picard tools, and per gene read counts for Ensembl (v75) gene annotations were computed. 
Expression levels in counts per million (CPM) were normalized and transformed using the VOOM procedure in the LIMMA R package. Surrogate variables to account sources of latent variation such as batch were calculated using the svaseq function from the [SVA package](https://bioconductor.org/packages/release/bioc/html/sva.html). Differential gene expression between races was performed with the LIMMA R package using the following linear model:

Y = β0 + β1×race + β2×sex + β3×Age + β4-14×SVA1:SVA11

where Y is log2 transformed gene expression, race is either African Americans or European Americans plus adjustments for sex,age and 11 surrogate variables




### Data Download

* [Raw Counts](https://www.dropbox.com/s/i5dthgl5c5ij5gd/Counts.csv?dl=0)
* [SVA Corrected CPMS](https://www.dropbox.com/s/mpgbhujqezts998/CPMS_SVA_corrected.RDS?dl=0)
* [Final eSet](https://www.dropbox.com/s/797rft3a7iihhmc/MAGNET_eset.RDS?dl=0)
* [Sample information (phenodata file)](https://www.dropbox.com/s/eihem5fbnkg7bpm/phenoData.csv?dl=0)


## eQTL

Expression quantitative trait locus analysis was performed using the [QTLtools](https://qtltools.github.io/qtltools/) package with adjustment for sex, race, and the first 3 genetic principal components and the 24 SVA-computed covariates.

```bash
for j in $(seq 1 30); do
QTLtools cis  --vcf final_maf10.vcf.bz --bed Phenotypes.bed.gz --out chunk_$j --cov $BASE/covars.txt --perm 1000 --chunk $j 30
done
```
