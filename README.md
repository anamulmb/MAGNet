# MAGNet

MAGNet is a collaborative group of investigators who use genomic approaches to understand human myocardial disease. Current projects include human myocardial expression, quantitative trait mapping and systems genetics.

NHLBI Grant number 1R01HL105993-01A1


# RNAseq
The Myocardial Applied Genomics Network (MAGNet; www.med.upenn.edu/magnet), collects and banks human cardiac tissue for genomic research.  All subjects or next of kin provided written informed consent for tissue donation and analyses and all study protocols were approved by relevant institutional review boards. Left ventricular free-wall tissue was harvested at the time of cardiac surgery from subjects with heart failure undergoing transplantation and from unused donor hearts with apparently normal function. The heart was perfused with cold cardioplegia prior to cardiectomy to arrest contraction and prevent ischemic damage, and tissue specimens were frozen in liquid nitrogen. Total RNA was extracted using the miRNeasy Kit (Qiagen) including DNAse treatment. RNA concentration and quality was determined using the NanoVue Plus™ spectrophotometer (GE Healthcare) and the Agilent 2100 RNA Nano Chip (Agilent).	

RNA sequencing libraries were prepared using the Illumina TruSeq stranded mRNA kit followed by the Nugen Ovation amplification kit. To avoid confounding by batch effects, libraries were randomly selected into pools of 32, and pools were sequenced on a Hiseq2500 to a depth of ~30 million 100-bp paired-end reads per biological sample. Fastq files were aligned against human reference (hg19/hGRC37) using the STAR aligner. (PMID: 23104886) Duplicate reads were removed using MarkDuplicates from Picard tools, and per gene read counts for Ensembl (v75) gene annotations were computed. Expression levels in counts per million (CPM) were normalized and transformed using VOOM in the LIMMA R package. Surrogate variables to account sources for latent variation such as batch were calculated using the svaseq function from the R SVA package. 




# eQTL
