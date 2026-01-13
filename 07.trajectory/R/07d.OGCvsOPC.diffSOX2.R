library(dplyr)
library(stringr)
# SOX2 peak OGC vs OPC 
projdir="/tscc/projects/PairedTF/"
SOX2_mtx <- read.table(file.path(projdir, "OGCvsOPC/mtx/SOX2_peak_q01.tab"), comment.char="@", header=TRUE, sep="\t")
SOX2_mtx %>% as.data.frame %>% rename(chr=X.chr) %>% mutate(peak=paste0(chr, "_", start, "_", end)) %>% tibble::column_to_rownames("peak") %>% select(SOX2_OGC, SOX2_OPC) -> SOX2_mat
SOX2_cpm <- edgeR::cpm(SOX2_mat)
SOX2_cpm %>% as.data.frame %>% mutate(log2FC=log2((SOX2_OGC+1e-16)/(SOX2_OPC+1e-16))) %>% arrange(desc(abs(log2FC))) -> SOX2_cpm_log2FC
SOX2_cpm_log2FC %>% filter(log2FC>2) -> SOX2_cpm_log2FC_pos
SOX2_cpm_log2FC %>% filter(log2FC<(-2)) -> SOX2_cpm_log2FC_neg
write.csv(SOX2_cpm_log2FC, file=file.path(projdir, "OGCvsOPC/SOX2_peakq01_cpm_log2FC.csv"), row.names=TRUE, quote=FALSE)
write.csv(SOX2_cpm_log2FC_pos, file=file.path(projdir, "OGCvsOPC/SOX2_peakq01_cpm_log2FC_pos.csv"), row.names=TRUE, quote=FALSE)
write.csv(SOX2_cpm_log2FC_neg, file=file.path(projdir, "OGCvsOPC/SOX2_peakq01_cpm_log2FC_neg.csv"), row.names=TRUE, quote=FALSE)
SOX2_cpm_log2FC_pos %>% tibble::rownames_to_column("peak") %>% tidyr::separate(peak, into=c("chr", "start", "end"), sep="_", remove=FALSE) %>% dplyr::select(chr, start, end, peak, log2FC) -> SOX2_cpm_log2FC_pos_bed
SOX2_cpm_log2FC_neg %>% tibble::rownames_to_column("peak") %>% tidyr::separate(peak, into=c("chr", "start", "end"), sep="_", remove=FALSE) %>% dplyr::select(chr, start, end, peak, log2FC) -> SOX2_cpm_log2FC_neg_bed
write.table(SOX2_cpm_log2FC_pos_bed, file=file.path(projdir, "OGCvsOPC/SOX2_peakq01_cpm_log2FC_pos.bed"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(SOX2_cpm_log2FC_neg_bed, file=file.path(projdir, "OGCvsOPC/SOX2_peakq01_cpm_log2FC_neg.bed"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

# add expression
SOX2_cpm_log2FC_pos_proximal <- read.table(file=file.path(projdir, "OGCvsOPC/SOX2_peakq01_cpm_log2FC_pos.proximal.peaks"), header=TRUE)
SOX2_RNA_combined_cpm <- readRDS(file=file.path(projdir, "SOX2.final.celltype.cpm.rds"))
SOX2_cpm_log2FC_pos_proximal$OPC_RNAexp <- SOX2_RNA_combined_cpm[match(SOX2_cpm_log2FC_pos_proximal$gene, rownames(SOX2_RNA_combined_cpm)), "OPC"]
SOX2_cpm_log2FC_pos_proximal$OGC_RNAexp <- SOX2_RNA_combined_cpm[match(SOX2_cpm_log2FC_pos_proximal$gene, rownames(SOX2_RNA_combined_cpm)), "OGC"]
write.csv(SOX2_cpm_log2FC_pos_proximal, file=file.path(projdir, "OGCvsOPC/SOX2_peakq01_cpm_log2FC_pos_proximal_withRNAexp.csv"), row.names=FALSE, quote=FALSE)

