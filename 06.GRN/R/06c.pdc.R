library(optparse)

op <- list(
  make_option(c("--tf"), type = "character"),
)

args <- parse_args(OptionParser(option_list = op))

# * configs
tf_name <- args$tf

# prepare data for pdc
suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))
options(stringsAsFactors = FALSE)

projdir <- "/tscc/projects/PairedTF"

TF_RNA <- readRDS(file=file.path(projdir, paste0(tf_name, ".final.reidents.rds")))
message("RNA:", "\t", nrow(TF_RNA), "\t", ncol(TF_RNA))

RNA_genes <- which(rowSums(GetAssayData(TF_RNA, slot="count")>0)>=1)
RNA <- TF_RNA[RNA_genes,]    
message("RNA filter:", "\t", nrow(RNA), "\t", ncol(RNA))

RNA_meta <- RNA@meta.data
RNA_meta$barcode <- str_split(rownames(RNA_meta), "_", simplify=TRUE)[,2]
RNA_meta$cell <- paste(RNA_meta$sampleID_DNA, RNA_meta$barcode, sep=".")
rownames(RNA_meta) <- RNA_meta$cell
RNA_meta %>% select(orig.ident, nCount_RNA, nFeature_RNA, percent.mt, seurat_clusters, sampleID, sampleID_DNA, cell, celltype, barcode) %>% mutate(logUMI=log(nCount_RNA)) -> meta_table

#RNA_mat <- RNA@assays$RNA@counts
RNA_mat <- RNA@assays$RNA@layers$counts
colnames(RNA_mat) <- meta_table$cell
rownames(RNA_mat) <- rownames(RNA)

TF_input <- file.path(projdir, paste0(tf_name, ".pmat.seurat4.rds"))
TF <- readRDS(TF_input)
message("TF:", "\t", nrow(TF), "\t", ncol(TF))

TF_meta <- read.csv(file.path(projdir, paste0(tf_name,".meta.table.forDNA.sampleID.csv")))
TF_meta %>% mutate(cell=paste(sampleID, barcode, sep=":")) -> TF_meta
TF$celltype <- TF_meta$celltype[match(colnames(TF), TF_meta$cell)]
Idents(TF) <- TF$celltype
#TF$sampleID <- str_split(colnames(TF), "[:_]", simplify=TRUE)[,1]
TF$barcode <- str_split(colnames(TF), ":", n=2, simplify=TRUE)[,2]
TF$cell <- paste(TF$sampleID, TF$barcode, sep=".")
#TF_meta <- TF@meta.data

TF_mat <- TF@assays$RNA@counts
colnames(TF_mat) <- as.character(TF$cell)

message("RNA and DNA cell same: ", identical(RNA_meta$cell, TF$cell))
message("DNA cell not in RNA cell: ", setdiff(TF$cell, RNA_meta$cell))
message("RNA cell not in DNA cell: ", paste(setdiff(RNA_meta$cell, TF$cell), collapse="\t"))
RNA_meta <- RNA_meta[match(TF$cell, rownames(RNA_meta)), ]
message("RNA meta filter and DNA cell same: ", identical(as.character(RNA_meta$cell), as.character(colnames(TF_mat))))
meta_table <- meta_table[match(TF$cell, rownames(meta_table)),]
message("meta table and DNA cell same: ", identical(as.character(meta_table$cell), as.character(colnames(TF_mat))))
RNA_mat <- RNA_mat[, match(TF$cell, colnames(RNA_mat))]
message("RNA filter and DNA cell same: ", identical(as.character(colnames(RNA_mat)), as.character(colnames(TF_mat))))

input_gene_peak <- file.path(projdir, tf_name, paste0(tf_name ,".peaks.union.tssUpDn1000k.gene.pairs"))
all_gene_peak <- read.table(input_gene_peak)
colnames(all_gene_peak) <- c("gene","peak")
message("gene peak pairs: ", nrow(all_gene_peak))
all_gene_peak <- all_gene_peak[which(all_gene_peak$gene %in% rownames(RNA_mat)),]
message("gene peak pairs filter: ", nrow(all_gene_peak))

dir.create(file.path(projdir, tf_name, "data"), recursive=TRUE, showWarnings=FALSE)
saveRDS(RNA_mat, file=file.path(projdir, paste0(tf_name, ".rna.mtx.rds")))
saveRDS(meta_table, file=file.path(projdir, paste0(tf_name, ".rna.meta.rds")))
saveRDS(TF_mat, file=file.path(projdir, paste0(tf_name, ".TF.mtx.rds")))
saveRDS(all_gene_peak, file=file.path(projdir, paste0(tf_name, ".gene_peak_pairs.rds")))


# prepare cpm matrix for pdc
RNA_mat <- readRDS(file.path(projdir, paste0(tf_name, ".rna.mtx.rds")))
meta_table <- readRDS(file.path(projdir, paste0(tf_name, ".rna.meta.rds")))
TF_mat <- readRDS(file.path(projdir, paste0(tf_name, ".TF.mtx.rds")))
#all_gene_peak <- readRDS(file.path(projdir, paste0(tf_name, ".gene_peak_pairs.rds")))

message("meta table and DNA cell same: ", identical(as.character(meta_table$cell), as.character(colnames(TF_mat))))
message("meta table and RNA cell same: ", identical(as.character(meta_table$cell), as.character(colnames(RNA_mat))))

RNA_match <- CreateSeuratObject(counts=RNA_mat, meta.data=meta_table)
TF_match <- CreateSeuratObject(counts=TF_mat, meta.data=meta_table)

RNA_match_count <- AggregateExpression(RNA_match, assays="RNA", slot="count", group.by="celltype")$RNA
RNA_match_cpm <- edgeR::cpm(RNA_match_count)

TF_match_count <- AggregateExpression(TF_match, assays="RNA", slot="count", group.by="celltype")$RNA
TF_match_cpm <- edgeR::cpm(TF_match_count)

rnaMat <- RNA_match_cpm[,celltypes]
TFMat <- TF_match_cpm[,celltypes]

message("RNA and DNA cell type same: ", identical(colnames(rnaMat), colnames(TFMat)))

saveRDS(rnaMat, file=file.path(projdir, paste0(tf_name, ".RNA.cpm.rds")))
saveRDS(TFMat, file=file.path(projdir, paste0(tf_name, ".DNA.cpm.rds")))


# for pearson pdc
### functions copy from Songpeng
scaleByRow <- function(mat,
                       is.center = TRUE,
                       is.scale = TRUE,
                       rmZeroSdRows = TRUE,
                       is.capped = TRUE,
                       up.hard = NULL,
                       down.hard = NULL,
                       up.quantile = 0.999,
                       down.quantile = 0.001
                       ) {
  sdRows <- apply(mat, 1, sd)
  zeroRows <- (sdRows == 0)
  if(any(zeroRows)) {
    message(sum(zeroRows), " has zero standard deviation.")
    if (rmZeroSdRows) {
      message("Will remove them.")
      mat <- mat[!zeroRows, ]
    } else {
      message("Data may have nan or NA after scaling.")
    }
  }
  r <- t(scale(t(mat), center = is.center, scale = is.scale))
  rownames(r) <- rownames(mat)
  colnames(r) <- colnames(mat)
  if (is.capped) {
    message("Will cap the data")
    up.thres <- if(!is.null(up.hard)) {
      up.hard
    } else {
      quantile(r, up.quantile)
    }
    down.thres <- if(!is.null(down.hard)) {
      down.hard
    } else {
      quantile(r, down.quantile)
    }
    message(
      paste("Restrict data between", down.thres, "and", up.thres, sep = " "))
    r[r >= up.thres] <- up.thres
    r[r <= down.thres] <- down.thres
  }
  return(r)
}

scaleUseRankByRow <- function(mat, tiesMethod = "average") {
  rank_mat <- apply(mat, 1, rank, ties.method = tiesMethod)
  r <- scaleByRow(mat = t(rank_mat), is.center = TRUE,
    is.scale = TRUE, rmZeroSdRows = TRUE,
    is.capped = FALSE)
  return(r)
}

fastPearCorOnRow.scaled <- function(mat1, mat2) {
  r <- rowSums(mat1 * mat2) / ncol(mat1)
  # ideally we do not need this
  # but our mat may be use hard threshold to
  # limit the too high / too low values
  # So add this just in case
  r[r >=1 ] <- 1
  r[r <= (-1)] <- (-1)
  return(r)
}

fastSpearCorOnRow.distinctRank <- function(mat1, mat2, tiesMethod = "average") {
  rankMat1 <- t(apply(mat1, 1, rank, ties.method = tiesMethod))
  rankMat2 <- t(apply(mat2, 1, rank, ties.method = tiesMethod))
  n <- ncol(mat1)
  diffRank <- 1 -  6 * rowSums((rankMat1-rankMat2)^2) / (n * (n^2-1))
  return(diffRank)
}

fastCorOnRow.chunk <- function(p1byc, p2byc, p12pairs,
                                chunkSize = 10e+05,
                                corMethod = c("spearman", "pearson")) {
  method <- match.arg(corMethod)
  # * filter p12pairs based on rownames of p1byc, p2byc
  upairs <- unique(p12pairs)
  message("Unique pairs: ", nrow(upairs))
  p1s <- unique(upairs[,1])
  p2s <- unique(upairs[,2])
  p1s.used <- p1s[p1s %in% rownames(p1byc)]
  p2s.used <- p2s[p2s %in% rownames(p2byc)]
  message(paste(length(p1s.used), "/", length(p1s),"ids are in p1byc."))
  message(paste(length(p2s.used), "/", length(p2s),"ids are in p2byc."))
  
  uupairs <- upairs[(upairs[,1] %in% p1s.used) & (upairs[,2] %in% p2s.used), ]
  message("After filtering names: ", nrow(uupairs), " pairs left.")

  # * use common columns in p1byc, p2byc
  com.cols <- intersect(colnames(p1byc), colnames(p2byc))
  message(length(com.cols), " common columns in p1byc and p2byc.")
  p1byc <- p1byc[, com.cols]
  p2byc <- p2byc[, com.cols]
  
  index.all <- seq_len(nrow(uupairs))
  chunks <- split(index.all, ceiling(seq_along(index.all) / chunkSize))
  if (method == "pearson") {
    message("Remind: p1byc and p2byc should be scaled for pearson.")
  }
  if (method == "spearman") {
    message("Remind: p1byc and p2byc should be scaled by rank for spearman")
    ## p1byc <- scaleUseRankByRow(p1byc)
    ## p2byc <- scaleUseRankByRow(p2byc)
  }
  rs <- lapply(chunks, function(index.chunk) {
    p12 <- uupairs[index.chunk, ]
    mat1 <- p1byc[p12[,1], ]
    mat2 <- p2byc[p12[,2], ]
    fastPearCorOnRow.scaled(mat1 = mat1, mat2 = mat2)
  })
  message("Merge results from all the chunks.")
  r <- unlist(rs)
  names(r) <- paste(uupairs[,1], uupairs[,2], sep = "@")
  return(r)
}

cor.method <- "pearson"
#cor.method <- "spearman"

rnaMat <- readRDS(file=file.path(projdir, paste0(tf_name, ".RNA.cpm.rds")))
TFMat <- readRDS(file=file.path(projdir, paste0(tf_name, ".DNA.cpm.rds")))
all_gene_peak <- readRDS(file.path(projdir, paste0(tf_name, ".gene_peak_pairs.rds")))

message("RNA and DNA cell type same: ", identical(colnames(rnaMat), colnames(TFMat)))

rnaMat <- log1p(rnaMat)
TFMat <- log1p(TFMat)

s.rna.mat <- if(cor.method == "pearson") {
	scaleByRow(mat = rnaMat)
} else {
	scaleUseRankByRow(mat = rnaMat)
}

s.TF.mat <- if(cor.method == "pearson") {
	scaleByRow(mat = TFMat)
} else {
	scaleUseRankByRow(mat = TFMat)
}

r <- fastCorOnRow.chunk(
	p1byc = s.rna.mat, p2byc = s.TF.mat, p12pairs = all_gene_peak,
	chunkSize = 10000, corMethod = cor.method)

rdf <- data.frame(
	pdc = names(r),
	cor = r
)
dir.create(file.path(projdir, tf_name, "results"), recursive=TRUE, showWarnings=FALSE)
write.table(x = rdf, file = file.path(projdir, tf_name, "/results/", paste0(tf_name, ".PGC.real.csv")),
sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)

# for shulffle
message("Will shuf the columns of TFMat")
set.seed(2025)
origCols <- colnames(TFMat)
TFMat <- TFMat[, sample(ncol(TFMat))]
colnames(TFMat) <- origCols

s.TF.mat <- if(cor.method == "pearson") {
	scaleByRow(mat = TFMat)
} else {
	scaleUseRankByRow(mat = TFMat)
}

r <- fastCorOnRow.chunk(
p1byc = s.rna.mat, p2byc = s.TF.mat, p12pairs = all_gene_peak,
chunkSize = 10000, corMethod = cor.method)

rdf <- data.frame(
	pdc = names(r),
	cor = r
)
write.table(x = rdf, file = file.path(projdir, tf_name, "/results/", paste0(tf_name, ".PGC.shuf.csv")),
	sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)


# summary cor
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(fitdistrplus)

density.theme <-   theme_bw() +
  theme(axis.text.x=element_text(colour = "black", size = 8),
    axis.text.y=element_text(colour = "black", size = 8),
    axis.title.x=element_text(colour = "black", size = 8),
    axis.title.y=element_text(colour = "black", size = 8))

loadCor <- function(f) {
  r <- data.table::fread(file = f, header = TRUE, sep = ",",
    data.table = FALSE)
  gene.peak  <- vapply(r$pdc, function(i) {
    t <- strsplit(i, split = "@", fixed = TRUE)
    t[[1]]
  }, c("a", "peak1"))
  r$gene <- gene.peak[1, ]
  r$peak <- gene.peak[2, ]
  return(r)
}

calSingleTailPvalue <- function(x, m, sd)  {
  p_lower <- pnorm(x, mean = m, sd = sd, lower.tail = TRUE, log.p = FALSE)
  p_1_lower <- 1-p_lower
  ind_pos <- (x >= 0)
  p <- ind_pos * p_1_lower + (1-ind_pos) * p_lower
  return(p)
}

corRealPear <- loadCor(file.path(projdir, tf_name, "/results/", paste0(tf_name, ".PGC.real.csv")))
corRdmPear <- loadCor(file.path(projdir, tf_name, "/results/", paste0(tf_name, ".PGC.shuf.csv")))

withr::with_pdf(new = file.path(projdir, tf_name, "/results/", paste0(tf_name, ".PGC.cor.pdf")),
code =  {
  data.pear <- rbind(corRealPear, corRdmPear)
  data.pear$class <- c(rep("Real.Pearson", nrow(corRealPear)),
    rep("RdmShuf.Pearson", nrow(corRdmPear)))
  p.pear <- ggplot(data.pear, aes(x = cor, colour = class)) +
    geom_density() +
    ggtitle("Pearson correlation") +
    density.theme
  print(p.pear)
}
)

fitRdmCor <- function(corRealFile, corRdmFile, corMethod = "spearman") {
	corReal <- loadCor(corRealFile)
	corReal <- corReal[!is.na(corReal$cor), ]
	corRdm <- loadCor(corRdmFile)
	corRdm <- corRdm[!is.na(corRdm$cor), ]
	fit.rdm.shuf <- fitdistrplus::fitdist(corRdm$cor, "norm", keepdata = FALSE)
	m.rdm.shuf <- fit.rdm.shuf$estimate[[1]]
	sd.rdm.shuf <- fit.rdm.shuf$estimate[[2]]
	para <- data.frame(class = "RdmShufPearson", mean = m.rdm.shuf, sd = sd.rdm.shuf)
	data.table::fwrite(para,
	  file.path(projdir, tf_name, "/results/", paste0(tf_name, ".fitnorm.rdm.shuf.", corMethod,".txt")),
	  quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
	# * get pvalue and fdr
	corReal$Pval <- calSingleTailPvalue(
	  x = corReal$cor, m = m.rdm.shuf, sd = sd.rdm.shuf)
	corReal$FDR <- p.adjust(corReal$Pval, method = "BH")
	data.table::fwrite(
	  corReal,
	  file = file.path(projdir, tf_name, "/results/", paste0(tf_name, ".PGC.cor.",corMethod,".fdr.tsv")),
	  quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
	r.alignv1 <- data.frame(
	  conns = gsub("@", "|", corReal$pdc),
	  pcc = corReal$cor,
	  class = rep("real", nrow(corReal)),
	  Pval = corReal$Pval,
	  FDR = corReal$FDR
)

data.table::fwrite(
  r.alignv1,
  file = file.path(projdir, tf_name, "/results/", paste0(tf_name, ".PGC.cor.",corMethod, ".fdr.alignv1.tsv")),
  quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
}

fitRdmCor(
	corRealFile = file.path(projdir, tf_name, "/results/", paste0(tf_name, ".PGC.real.csv")),
	corRdmFile = file.path(projdir, tf_name, "/results/", paste0(tf_name, ".PGC.shuf.csv")),
	corMethod = "pearson")
}

r.align1 <- data.table::fread(file.path(projdir, tf_name, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".fdr.alignv1.tsv")), header = TRUE, sep = "\t", quote = "", data.table = FALSE)
r.align1 %>% filter(Pval<0.005&abs(pcc)>0.6) %>% mutate(gene=str_split(conns, "\\|", simplify=TRUE)[,1], peak=str_split(conns, "\\|", simplify=TRUE)[,2]) -> r.align1_filt
r.align1_filt %>% filter(pcc>0) -> pos.pdc
r.align1_filt %>% filter(pcc<0) -> neg.pdc
message("Total sig pairs: ", nrow(r.align1_filt), "\t", "pos pairs: ", nrow(pos.pdc), "\t", "neg pairs: ", nrow(neg.pdc))
data.table::fwrite(r.align1_filt, file = file.path(projdir, tf_name, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.tsv")), quote = FALSE, sep = "\t",
col.names = TRUE, row.names = FALSE)
data.table::fwrite(pos.pdc, file = file.path(projdir, tf_name, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".pos.tsv")), quote = FALSE, sep = "\t",
col.names = TRUE, row.names = FALSE)
data.table::fwrite(neg.pdc, file = file.path(projdir, tf_name, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".neg.tsv")), quote = FALSE, sep = "\t",
col.names = TRUE, row.names = FALSE)

# vocano
cols <- c("pos"="#ffad73", "neg"="#26b3ff", "ns"="grey") 
sizes <- c("pos"=3, "neg"=3, "ns"=3) 
alphas <- c("pos"=1, "neg"=1, "ns"=1)

r.align1 <- r.align1 %>% mutate(Type=case_when(pcc>0.6 & Pval<0.005 ~ "pos", pcc<(-0.6) & Pval<0.005 ~ "neg", TRUE ~ "ns")) %>% mutate(Type=forcats::fct_relevel(Type, "pos", "neg", "ns")) 

pt_TF_PGC <- ggplot(r.align1, aes(x=pcc, y=-log10(Pval))) + 
geom_point(aes(fill=Type, color=Type, size=Type, alpha=Type)) + scale_fill_manual(values=cols) + scale_color_manual(values=cols) + scale_size_manual(values=sizes) + scale_alpha_manual(values=alphas) + 
geom_hline(yintercept=-log10(0.005), linetype="dashed") + geom_vline(xintercept=c(0.6, -0.6), linetype="dashed") + 
scale_x_continuous(breaks=c(-1, -0.6, 0, 0.6, 1), limits=c(-1, 1)) + ggtitle(paste0(tf_name, " peak and gene correlation")) + 
xlab("PCC") + ylab("-log10(Pval)") + theme_classic() + 
theme(legend.title=element_text(size=12), legend.text=element_text(size=10), axis.text.x=element_text(colour="black", size=16), 
      axis.text.y=element_text(colour="black", size=16), axis.title=element_text(colour="black", size=18), plot.title=element_text(hjust=0.5, size=20, color="black"))
pt_TF_PGC.label <- pt_TF_PGC +  
annotate("text", x=0.5, y=1.2, label=paste0("# of pos pairs: ", nrow(pos.pdc)), color="black", size=5) + 
annotate("text", x=-0.5, y=1.2, label=paste0("# of neg pairs: ", nrow(neg.pdc)), color="black", size=5)

library("ggrastr")
pt_TF_PGC.raster <- rasterize(pt_TF_PGC.label, layers='Point', dpi=300)

ggsave(pt_TF_PGC.raster, file=file.path(projdir, tf_name, "/results/", paste0("pt_", tf_name, "_PGC.label.pdf")), width=10, height=8)
ggsave(pt_TF_PGC.label, file=file.path(projdir, tf_name, "/results/", paste0("pt_", tf_name, "_PGC.label.png")), width=10, height=8)

# generate bedpe file
gene_tss <- read.table("/tscc/projects/annot/gencode.vM23.gene.tssUpDn1k.bed")
r.align1 %>% tidyr::separate(conns, into=c("gene", "peak"), sep="\\|", remove=FALSE) %>% tidyr::separate(peak, into=c("chr2", "start2", "end2"), sep="[:-]", remove=FALSE) -> r.align1_sep
r.align1_sep$chr1 <- gene_tss$V1[match(r.align1_sep$gene, gene_tss$V7)]
r.align1_sep$start1 <- gene_tss$V2[match(r.align1_sep$gene, gene_tss$V7)]
r.align1_sep$end1 <- gene_tss$V3[match(r.align1_sep$gene, gene_tss$V7)]
r.align1_sep$strand1 <- "."
r.align1_sep$strand2 <- "."

r.align1_sep %>% filter(Pval<0.005&abs(pcc)>0.6) %>% dplyr::select(chr1,start1,end1,chr2,start2,end2,conns,Pval,strand1,strand2) -> r.align1_bedpe
r.align1_sep %>% filter(Pval<0.005, pcc>0.6) %>% dplyr::select(chr1,start1,end1,chr2,start2,end2,conns,Pval,strand1,strand2) -> pos.pdc_bedpe
r.align1_sep %>% filter(Pval<0.005, pcc<(-0.6)) %>% dplyr::select(chr1,start1,end1,chr2,start2,end2,conns,Pval,strand1,strand2) -> neg.pdc_bedpe
write.table(r.align1_bedpe, file=file.path(projdir, tf_name, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.bedpe")), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(pos.pdc_bedpe, file=file.path(projdir, tf_name, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".pos.bedpe")), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(neg.pdc_bedpe, file=file.path(projdir, tf_name, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".neg.bedpe")), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)



