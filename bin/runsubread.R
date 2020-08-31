library(Rsubread)

gtffile='/home/mmorley/NGSshare/hg19_data/Homo_sapiens.GRCh37.75.gtf'

(bamfiles <- list.files('~/magnet2/MAGNET_final_corrected/bamfiles',pattern=".dedupped.bam$",full.names=TRUE))


cts.dedup <- featureCounts(files=bamfiles,annot.ext=gtffile,strandSpecific=2,ignoreDup=T,isGTFAnnotationFile=TRUE,GTF.featureType="exon",
                           GTF.attrType="gene_id",nthreads=40,isPairedEnd=TRUE)
save(cts.dedup,file='./data/subread_counts_allgood.Rdata')

