library(dplyr)
library(stringr)

projdir <- "/tscc/projects/PairedTF"

# ATAC peak OGC vs OPC
ATAC_mtx <- read.table(file.path(projdir, "OGCvsOPC/ATAC_OPC_OGC.tab"), comment.char="@", header=TRUE, sep="\t")
ATAC_mtx %>% as.data.frame %>% rename(chr=X.chr, ATAC_OPC=X326_OPC_NN.srt, ATAC_OGC=X327_Oligo_NN.srt) %>% mutate(peak=paste0(chr, "_", start, "_", end)) %>% distinct %>% tibble::column_to_rownames("peak") %>% select(ATAC_OGC, ATAC_OPC) -> ATAC_mat
ATAC_cpm <- edgeR::cpm(ATAC_mat)
ATAC_cpm %>% as.data.frame %>% mutate(log2FC=log2((ATAC_OGC+1e-16)/(ATAC_OPC+1e-16))) %>% arrange(desc(abs(log2FC))) -> ATAC_cpm_log2FC
ATAC_cpm_log2FC %>% filter(log2FC>2) -> ATAC_cpm_log2FC_pos
ATAC_cpm_log2FC %>% filter(log2FC<(-2)) -> ATAC_cpm_log2FC_neg
write.csv(ATAC_cpm_log2FC, file=file.path(projdir, "OGCvsOPC/ATAC_peak_cpm_log2FC.csv"), row.names=TRUE, quote=FALSE)
write.csv(ATAC_cpm_log2FC_pos, file=file.path(projdir, "OGCvsOPC/ATAC_peak_cpm_log2FC_pos.csv"), row.names=TRUE, quote=FALSE)
write.csv(ATAC_cpm_log2FC_neg, file=file.path(projdir, "OGCvsOPC/ATAC_peak_cpm_log2FC_neg.csv"), row.names=TRUE, quote=FALSE)
ATAC_cpm_log2FC_pos %>% tibble::rownames_to_column("peak") %>% tidyr::separate(peak, into=c("chr", "start", "end"), sep="_", remove=FALSE) %>% dplyr::select(chr, start, end, peak, log2FC) -> ATAC_cpm_log2FC_pos_bed
ATAC_cpm_log2FC_neg %>% tibble::rownames_to_column("peak") %>% tidyr::separate(peak, into=c("chr", "start", "end"), sep="_", remove=FALSE) %>% dplyr::select(chr, start, end, peak, log2FC) -> ATAC_cpm_log2FC_neg_bed
write.table(ATAC_cpm_log2FC_pos_bed, file=file.path(projdir, "OGCvsOPC/ATAC_peak_cpm_log2FC_pos.bed"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(ATAC_cpm_log2FC_neg_bed, file=file.path(projdir, "OGCvsOPC/ATAC_peak_cpm_log2FC_neg.bed"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

# add expression 
ATAC_cpm_log2FC_pos_proximal <- read.table(file=file.path(projdir, "OGCvsOPC/ATAC_peak_cpm_log2FC_pos.proximal.peaks"), header=TRUE)
ATAC_RNA_cpm <- readRDS(file="/tscc/projects/CEMBA2/allen.l2.cpm.ds200.rds")
ATAC_cpm_log2FC_pos_proximal$OPC_RNAexp <- ATAC_RNA_cpm[match(ATAC_cpm_log2FC_pos_proximal$gene, rownames(ATAC_RNA_cpm)), "OPC_NN"]
ATAC_cpm_log2FC_pos_proximal$OGC_RNAexp <- ATAC_RNA_cpm[match(ATAC_cpm_log2FC_pos_proximal$gene, rownames(ATAC_RNA_cpm)), "Oligo_NN"]
write.csv(ATAC_cpm_log2FC_pos_proximal, file=file.path(projdir, "OGCvsOPC/ATAC_peak_cpm_log2FC_pos_proximal_withRNAexp.csv"), row.names=FALSE, quote=FALSE)
