library(dplyr)
library(ggplot2)
library(stringr)


# refine gene hubs with at least two genes
projdir <- "/tscc/projects/PairedTF"
celltype2subclass <- read.table(file.path(projdir, "celltype_subclass.txt"), header=TRUE, sep="\t")
celltype2subclass2color <- read.csv(file.path(projdir, "celltype2subclass2color.csv"))
celltype2subclass$HexColor <- celltype2subclass2color$HexColor[match(celltype2subclass$celltype1, celltype2subclass2color$celltype)]
celltype2subclass$HexColor19 <- celltype2subclass2color$HexColor19[match(celltype2subclass$celltype1, celltype2subclass2color$celltype)]

genehub_list <- list()
genehub_multigenes_list <- list()
for (celltype in celltype2subclass$celltype1) {
  celltype_name <- celltype2subclass$celltype[which(celltype2subclass$celltype1==celltype)]
  genehubs <- read.table(file.path(projdir, paste0("RNAPII/rose/", celltype_name, "/", celltype_name, "_SuperStitched_REGION_TO_GENE.txt")), sep="\t", header=TRUE)
  genehubs %>% mutate(Region_size=STOP-START, Celltype=celltype) %>% rowwise() %>% mutate(Gene_num=sum(str_split(OVERLAP_GENES, ",", simplify=TRUE)!="")) -> genehubs
  genehubs %>% filter(Gene_num>=2) -> genehubs_multigenes
  genehub_list[[celltype]] <- genehubs
  genehub_multigenes_list[[celltype]] <- genehubs_multigenes
}

genehubs_multigenes_all <- do.call(rbind, genehub_multigenes_list) %>% as.data.frame %>% mutate(Celltype=factor(Celltype, levels=unique(Celltype)))
write.csv(genehubs_multigenes_all, file=file.path(projdir, "RNAPII/rose/genehubs_multigenes_all.csv"), row.names=FALSE, quote=TRUE)

pt_hubs_num <- ggplot(genehubs_multigenes_all, aes(x=Celltype, fill=Celltype)) +  
  geom_bar() + scale_y_continuous(breaks=c(0, 500, 1000, 1500), limits=c(0, 2000), labels=c(0, 500, 1000, 1500)) +
  scale_fill_manual(values=celltype2subclass$HexColor19[which(celltype2subclass$celltype1 %in% genehubs_multigenes_all$Celltype)]) +
  theme_classic() + ggtitle("") + xlab("") + ylab("# of gene hubs") + 
  theme(legend.position="none", axis.title=element_text(color="black", size=10), axis.text=element_text(color="black", size=10))
ggsave(pt_hubs_num, file=file.path(projdir, "RNAPII/rose/hubs_multigenes_num.pdf"), height=5, width=15)


pt_hubs_size <- ggplot(genehubs_multigenes_all, aes(x=Celltype, y=Region_size, fill=Celltype)) + 
  scale_y_continuous(trans=scales::pseudo_log_trans(sigma=1, base=10), breaks=(10**c(1:6)), labels=scales::trans_format('log10',scales::math_format(10^.x)), expand=expansion(0)) + 
  geom_violin(width=0.8) + geom_boxplot(width=0.2, color="grey", alpha=0.2, outlier.shape = NA) + 
  scale_fill_manual(values=celltype2subclass$HexColor19[which(celltype2subclass$celltype1 %in% genehubs_multigenes_all$Celltype)]) + 
  theme_classic() + ggtitle("") + xlab("") + ylab("Size of hub region") + 
  theme(legend.position="none", plot.title=element_text(size=15, color="black", hjust=0.5), axis.title.y=element_text(size=10, color="black"), axis.text=element_text(size=10, color="black")) 
ggsave(pt_hubs_size, file=file.path(projdir, "RNAPII/rose/hubs_multigenes_size.pdf"), height=5, width=15)

pt_hubs_genenum <- ggplot(genehubs_multigenes_all, aes(x=Celltype, y=Gene_num, fill=Celltype)) + 
  scale_y_continuous(trans=scales::pseudo_log_trans(sigma=1, base=2), breaks=(2**c(1:6)), labels=scales::trans_format('log2',scales::math_format(2^.x)), expand=expansion(0)) + 
  geom_violin(width=0.8) + geom_boxplot(width=0.2, color="grey", alpha=0.2, outlier.shape = NA) + 
  scale_fill_manual(values=celltype2subclass$HexColor19[which(celltype2subclass$celltype1 %in% genehubs_multigenes_all$Celltype)]) + 
  theme_classic() + ggtitle("") + xlab("") + ylab("# of overlapped genes") + 
  theme(legend.position="none", plot.title=element_text(size=15, color="black", hjust=0.5), axis.title.y=element_text(size=10, color="black"), axis.text=element_text(size=10, color="black")) 
ggsave(pt_hubs_genenum, file=file.path(projdir, "RNAPII/rose/hubs_multigenes_genenum.pdf"), height=5, width=15)

# expresion of genes in hubs
suppressPackageStartupMessages(library("clusterProfiler"))
library("org.Mm.eg.db")

genes_list <- list()
for (celltype in names(genehub_multigenes_list)) {
  genehubs <- genehub_multigenes_list[[celltype]]
  genes <- unique(unlist(str_split(genehubs$OVERLAP_GENES, ",")))
  genes <- genes[which(genes!="")]
  genes_list[[celltype]] <- genes
}

genesymbol_list <- list()
for (celltype in names(genes_list)) {
  genes <- genes_list[[celltype]]
  genesymbols <- unique(bitr(genes, fromType="REFSEQ", toType=c("SYMBOL"), OrgDb="org.Mm.eg.db")$SYMBOL)
  genesymbol_list[[celltype]] <- genesymbols
}

RNAPII_seurat <- readRDS(file.path(projdir, "21.seurat/cluster/CTCF_RNAPolII.rds"))
RNAPII_cells <- colnames(RNAPII_seurat)[grepl("0[7-9]$|1[0-2]$", colnames(RNAPII_seurat))]
RNAPII_RNA <- subset(RNAPII_seurat, cells=RNAPII_cells)
RNAPII_RNA_count <- AggregateExpression(RNAPII_RNA, assays="RNA", slot="count", group.by="celltype")$RNA
RNAPII_RNA_cpm <- edgeR::cpm(RNAPII_RNA_count)
RNAPII_RNA_cpm_filterCT <- RNAPII_RNA_cpm[,c("CLAGL", "ITL6GL", "ITL5GL", "ITL45GL", "ITL23GL", "PPPGL", "ETL5GL", "L6bGL", "CTL6GL", "NPL5GL", "VIPGA", "PVGA", "SSTGA", "ASC", "OPC", "OGC", "VLMC", "Endo", "MGL")]
RNAPII_RNA_logcpm <- log1p(RNAPII_RNA_cpm_filterCT)


geneexpress_list <- list()
for (celltype in names(genesymbol_list)) {
  geneexpress <- data.frame(Celltype=celltype, Gene=genesymbol_list[[celltype]])
  geneexpress$RNA <- RNAPII_RNA_cpm_filterCT[match(geneexpress$Gene, rownames(RNAPII_RNA_cpm_filterCT)), celltype]
  geneexpress$RNAlog <- RNAPII_RNA_logcpm[match(geneexpress$Gene, rownames(RNAPII_RNA_logcpm)), celltype]
  geneexpress_list[[celltype]] <- geneexpress
}

geneexpress_all <- do.call(rbind, geneexpress_list) %>% as.data.frame %>% mutate(Celltype=factor(Celltype, levels=unique(Celltype)))

pt_hubs_geneexpress <- ggplot(geneexpress_all, aes(x=Celltype, y=RNAlog, fill=Celltype)) + 
  geom_violin(width=0.8) + geom_boxplot(width=0.2, color="grey", alpha=0.2, outlier.shape = NA) + 
  scale_fill_manual(values=celltype2subclass$HexColor19[which(celltype2subclass$celltype1 %in% genehubs_multigenes_all$Celltype)]) + 
  theme_classic() + ggtitle("") + xlab("") + ylab("Expression (logCPM) of overlapped genes") + 
  theme(legend.position="none", plot.title=element_text(size=15, color="black", hjust=0.5), axis.title.y=element_text(size=10, color="black"), axis.text=element_text(size=10, color="black")) 
ggsave(pt_hubs_geneexpress, file=file.path(projdir, "RNAPII/rose/hubs_multigenes_geneexpress.pdf"), height=5, width=15)

gene_combns <- list()
for (i in 1:length(genehubs_multigenes$OVERLAP_GENES)) {
  ovlp_genes <- genehubs_multigenes$OVERLAP_GENES[i]
  gene_combn <- t(combn(unlist(str_split(ovlp_genes, ",")),2))
  gene_combn_sort <- t(apply(gene_combn,1,sort))
  gene_combns[[i]] <- gene_combn_sort
}
gene_combns_all <- do.call(rbind, gene_combns) %>% as.data.frame %>% dplyr::rename("Gene1"="V1", "Gene2"="V2") %>% distinct()
gene_symbol <- bitr(c(gene_combns_all$Gene1, gene_combns_all$Gene2), fromType="REFSEQ", toType=c("SYMBOL"), OrgDb="org.Mm.eg.db")
gene_combns_all$Symbol1 <- gene_symbol$SYMBOL[match(gene_combns_all$Gene1, gene_symbol$REFSEQ)]
gene_combns_all$Symbol2 <- gene_symbol$SYMBOL[match(gene_combns_all$Gene2, gene_symbol$REFSEQ)]
gene_combns_all %>% filter(!is.na(Symbol1) & !is.na(Symbol2)) %>% filter((Symbol1 %in% rownames(RNAPII_RNA_logcpm)) & (Symbol2 %in% rownames(RNAPII_RNA_logcpm))) -> gene_combns_filter
gene_combns_filter %>% dplyr::select(Symbol1, Symbol2) %>% filter(Symbol1!=Symbol2) %>% distinct() -> gene_combns_filter
gene_combns_filter <- t(apply(gene_combns_filter, 1, sort)) %>% as.data.frame %>% dplyr::rename("Symbol1"="V1", "Symbol2"="V2")

gene_combns_cor <- data.frame(Gene1=gene_combns_filter$Symbol1, Gene2=gene_combns_filter$Symbol2, 
  Cor=apply(gene_combns_filter[,c("Symbol1", "Symbol2")], 1, function(x) cor(RNAPII_RNA_logcpm[x[1],], RNAPII_RNA_logcpm[x[2],])), 
  Pval=apply(gene_combns_filter[,c("Symbol1", "Symbol2")], 1, function(x) cor.test(RNAPII_RNA_logcpm[x[1],], RNAPII_RNA_logcpm[x[2],])$p.value))

set.seed(2025)
gene_combns_filter %>% mutate(Pair=paste(Symbol1, Symbol2, sep="_")) -> gene_combns_filter

genes_sample <- sample(rownames(RNAPII_RNA_logcpm), length(unique(c(gene_combns_filter$Symbol1, gene_combns_filter$Symbol2))))
allgene_combn <- t(combn(genes_sample,2))
combn_sample <- allgene_combn[sample(1:nrow(allgene_combn), nrow(gene_combns_filter)),]
combn_sample <- t(apply(combn_sample,1,sort)) %>% as.data.frame %>% dplyr::rename("Symbol1"="V1", "Symbol2"="V2") %>% mutate(Pair=paste(Symbol1, Symbol2, sep="_")) %>% filter(!(Pair %in% gene_combns_filter$Pair))
gene_combns_cor_sample <- data.frame(Gene1=combn_sample$Symbol1, Gene2=combn_sample$Symbol2, 
  Cor=apply(combn_sample[,c("Symbol1", "Symbol2")], 1, function(x) cor(RNAPII_RNA_logcpm[x[1],], RNAPII_RNA_logcpm[x[2],])), 
  Pval=apply(combn_sample[,c("Symbol1", "Symbol2")], 1, function(x) cor.test(RNAPII_RNA_logcpm[x[1],], RNAPII_RNA_logcpm[x[2],])$p.value))

gene_combns_cor$Type <- "Gene pairs within hubs"
gene_combns_cor_sample$Type <- "Random gene pairs"
gene_combns_cor_all <- rbind(gene_combns_cor, gene_combns_cor_sample)
cor_pval <- t.test(gene_combns_cor$Cor, gene_combns_cor_sample$Cor)$p.value
cor_pval <- wilcox.test(gene_combns_cor$Cor, gene_combns_cor_sample$Cor)$p.value

pt_hubs_cor <- ggplot(gene_combns_cor_all, aes(x=Type, y=Cor, fill=Type)) + 
  geom_violin(width=0.8) + geom_boxplot(width=0.2, color="grey", alpha=0.2, outlier.shape = NA) + 
  scale_fill_manual(values=celltype2subclass$HexColor19[sample(1:19,2)]) + 
  theme_classic() + ggtitle("") + xlab("") + ylab("Correlation of gene expression across cell types") + 
  theme(legend.position="none", plot.title=element_text(size=15, color="black", hjust=0.5), axis.title.y=element_text(size=10, color="black"), axis.text=element_text(size=10, color="black")) 
pt_hubs_cor_label <- pt_hubs_cor + annotate(geom="segment", x=1, xend=2, y=1.1, yend=1.1, colour="black") + geom_text(aes(label="p<2.2e-16", x=1.5, y=1.15, size=5)) 

ggsave(pt_hubs_cor_label, file=file.path(projdir, "RNAPII/rose/hubs_multigenes_gene_cor.pdf"), height=6, width=6)


genehubs_multigenes_all %>% filter(Gene_num>=2) %>% group_by(Celltype,OVERLAP_GENES) %>% summarise(num=n()) %>% group_by(OVERLAP_GENES) %>% summarise(sum=n()) %>% arrange(desc(sum)) -> genehubs_multigenes_num

multigenes_list <- list()
for (i in 1:nrow(genehubs_multigenes_num)) {
  ovlp_genes <- unlist(str_split(genehubs_multigenes_num$OVERLAP_GENES[i], ","))
  multigenes <- data.frame(ID=rep(i, length(ovlp_genes)), Gene=ovlp_genes, Freq=rep(genehubs_multigenes_num$sum[i], length(ovlp_genes)), Num=rep(length(ovlp_genes), length(ovlp_genes)))
  multigenes_list[[i]] <- multigenes
}

multigenes_all <- do.call(rbind, multigenes_list) %>% as.data.frame 
multigenes_symbol <- bitr(multigenes_all$Gene, fromType="REFSEQ", toType=c("SYMBOL"), OrgDb="org.Mm.eg.db")
multigenes_all$Symbol <- multigenes_symbol$SYMBOL[match(multigenes_all$Gene, multigenes_symbol$REFSEQ)]
write.csv(multigenes_all, file=file.path(projdir, "RNAPII/rose/multigenes_multigenes_all.csv"), row.names=FALSE, quote=FALSE)

genehubs_multigenes_all %>% dplyr::select(CHROM, START, STOP, Celltype) %>% dplyr::rename("chr"="CHROM", "start"="START", "end"="STOP") %>% arrange(chr, start, end) -> genehubs_multigenes_all_bed

df2gr <- function(df){
  gr <- GRanges(df$chr, IRanges(df$start, df$end))
  mcols(gr)$rank <- df$rank
  mcols(gr)$celltype <- df$celltype
  return(gr)
}

suppressPackageStartupMessages(library("GenomicRanges"))
peak.list <- lapply(names(genehub_multigenes_list), function(i){
  df <- genehub_multigenes_list[[i]] %>% dplyr::select(CHROM, START, STOP, stitchedPeakRank, Celltype) %>% dplyr::rename("chr"="CHROM", "start"="START", "end"="STOP", "rank"="stitchedPeakRank", "celltype"="Celltype")
  p.gr <- df2gr(df)
})

merged.gr <- do.call(c, peak.list)
#merged.gr <- nonOverlappingGR(merged.gr, by = "rank", decreasing = FALSE)
merged.gr <- GenomicRanges::reduce(merged.gr)
outUnion <- as.data.frame(merged.gr)
write.table(outUnion, file=file.path(projdir, "RNAPII/rose/genehubs_multigenes_merge.bed"), sep="\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# RNAPII super enhacer for all celltype
ucsc_refseq <- read.table("/home/kaw033/softwares/ROSE/annotation/mm10_refseq.ucsc", header=TRUE, comment.char="@")
ucsc_refseq %>% dplyr::select(chrom, txStart, txEnd, name2) %>% filter(!grepl("random|chrUn", chrom)) %>% dplyr::rename("start"="txStart", "end"="txEnd", "gene"="name2") %>% mutate(size=end-start) -> ucsc_refseq_bed
ucsc_refseq_bed %>% arrange(gene, desc(size)) %>% distinct(gene, .keep_all=TRUE) -> ucsc_refseq_bed_filter

RNAPII_list <- list()
for (celltype in celltype2subclass$celltype1) {
  celltype_name <- celltype2subclass$celltype[which(celltype2subclass$celltype1==celltype)]
  RNAPII_file <- file.path(projdir, paste0("RNAPII/rose/", celltype_name, "/", celltype_name, "_SuperStitched_REGION_TO_GENE.txt"))
  if (file.exists(RNAPII_file) & file.exists(ATAC_file)) {
    RNAPII_rose <- read.table(RNAPII_file, sep="\t", header=TRUE)
    RNAPII_list[[celltype]] <- RNAPII_rose
  }
}

RNAPII_SE_list <- list()
for (celltype in names(ATAC_list)) {
  RNAPII_rose <- RNAPII_list[[celltype]]
  RNAPII_rose %>% dplyr::select(CHROM, START, STOP) %>% dplyr::rename("chrom"="CHROM", "start"="START", "end"="STOP") -> RNAPII_rose_bed
  RNAPII_refseq_ovlp <- bedtoolsr::bt.intersect(RNAPII_rose_bed, ucsc_refseq_bed_filter, wo=TRUE) %>% distinct %>% dplyr::rename(chrom1=V1, start1=V2, end1=V3, chrom2=V4, start2=V5, end2=V6, gene=V7, size=V8, ovlp=V9) %>% 
    mutate(RNAPII=paste0(chrom1, ":", start1, "-", end1)) %>% arrange(chrom1, start1)
  RNAPII_refseq_ovlp %>% group_by(RNAPII) %>% summarise(total_ovlp=sum(as.numeric(ovlp))) %>% mutate(chrom=str_split(RNAPII, ":|-", simplify=TRUE)[,1], 
    start=str_split(RNAPII, ":|-", simplify=TRUE)[,2], end=str_split(RNAPII, ":|-", simplify=TRUE)[,3], size=as.numeric(end)-as.numeric(start), ratio=total_ovlp/size) %>% 
    arrange(ratio) %>% as.data.frame -> RNAPII_refseq_ovlp_sum
  RNAPII_refseq_ovlp_sum %>% filter(ratio>0.05|size<1000) %>% pull(RNAPII) -> RNAPII_refseq_RMratio_peaks
  RNAPII_rose_bed %>% mutate(RNAPII=paste0(chrom, ":", start, "-", end)) %>% filter(!(RNAPII %in% RNAPII_refseq_RMratio_peaks)) %>% dplyr::select(chrom, start, end) -> RNAPII_rose_RMratio_bed
  RNAPII_SE_list[[celltype]] <- RNAPII_rose_RMratio_bed
}

RNAPII_SE_form <- list()
for (celltype in names(RNAPII_SE_list)) {
  RNAPII_SE <- RNAPII_SE_list[[celltype]] %>% mutate(SE=paste0(chrom, ":", start, "-", end), celltype=celltype)
  RNAPII_SE_form[[celltype]] <- RNAPII_SE
}
RNAPII_SEs <- do.call(rbind, RNAPII_SE_form) %>% as.data.frame
write.table(RNAPII_SEs, file=file.path(projdir, "RNAPII/rose/RNAPII_SEs.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

