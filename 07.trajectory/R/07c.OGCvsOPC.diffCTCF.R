library(dplyr)
library(stringr)

# CTCF peak OGC vs OPC
projdir="/tscc/projects/PairedTF/"
CTCF_mtx <- read.table(file.path(projdir, "OGCvsOPC/mtx/CTCF_peak_q01.tab"), comment.char="@", header=TRUE, sep="\t")
CTCF_mtx %>% as.data.frame %>% rename(chr=X.chr) %>% mutate(peak=paste0(chr, "_", start, "_", end)) %>% tibble::column_to_rownames("peak") %>% select(CTCF_OGC, CTCF_OPC) -> CTCF_mat
CTCF_cpm <- edgeR::cpm(CTCF_mat)
CTCF_cpm %>% as.data.frame %>% mutate(log2FC=log2((CTCF_OGC+1e-16)/(CTCF_OPC+1e-16))) %>% arrange(desc(abs(log2FC))) -> CTCF_cpm_log2FC
CTCF_cpm_log2FC %>% filter(log2FC>2) -> CTCF_cpm_log2FC_pos
CTCF_cpm_log2FC %>% filter(log2FC<(-2)) -> CTCF_cpm_log2FC_neg
write.csv(CTCF_cpm_log2FC, file=file.path(projdir, "OGCvsOPC/CTCF_peakq01_cpm_log2FC.csv"), row.names=TRUE, quote=FALSE)
write.csv(CTCF_cpm_log2FC_pos, file=file.path(projdir, "OGCvsOPC/CTCF_peakq01_cpm_log2FC_pos.csv"), row.names=TRUE, quote=FALSE)
write.csv(CTCF_cpm_log2FC_neg, file=file.path(projdir, "OGCvsOPC/CTCF_peakq01_cpm_log2FC_neg.csv"), row.names=TRUE, quote=FALSE)
CTCF_cpm_log2FC_pos %>% tibble::rownames_to_column("peak") %>% tidyr::separate(peak, into=c("chr", "start", "end"), sep="_", remove=FALSE) %>% dplyr::select(chr, start, end, peak, log2FC) -> CTCF_cpm_log2FC_pos_bed
CTCF_cpm_log2FC_neg %>% tibble::rownames_to_column("peak") %>% tidyr::separate(peak, into=c("chr", "start", "end"), sep="_", remove=FALSE) %>% dplyr::select(chr, start, end, peak, log2FC) -> CTCF_cpm_log2FC_neg_bed
write.table(CTCF_cpm_log2FC_pos_bed, file=file.path(projdir, "OGCvsOPC/CTCF_peakq01_cpm_log2FC_pos.bed"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(CTCF_cpm_log2FC_neg_bed, file=file.path(projdir, "OGCvsOPC/CTCF_peakq01_cpm_log2FC_neg.bed"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

# add expression 
CTCF_cpm_log2FC_pos_proximal <- read.table(file=file.path(projdir, "OGCvsOPC/CTCF_peakq01_cpm_log2FC_pos.proximal.peaks"), header=TRUE)
CTCF_RNA_cpm <- readRDS(file=file.path(projdir, "CTCF.final.celltypes.cpm.rds"))
CTCF_cpm_log2FC_pos_proximal$OPC_RNAexp <- CTCF_RNA_cpm[match(CTCF_cpm_log2FC_pos_proximal$gene, rownames(CTCF_RNA_cpm)), "OPC"]
CTCF_cpm_log2FC_pos_proximal$OGC_RNAexp <- CTCF_RNA_cpm[match(CTCF_cpm_log2FC_pos_proximal$gene, rownames(CTCF_RNA_cpm)), "OGC"]
write.csv(CTCF_cpm_log2FC_pos_proximal, file=file.path(projdir, "OGCvsOPC/CTCF_peakq01_cpm_log2FC_pos_proximal_withRNAexp.csv"), row.names=FALSE, quote=FALSE)
