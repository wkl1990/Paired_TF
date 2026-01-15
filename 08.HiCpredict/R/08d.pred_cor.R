# compare HiC prediction
library(dplyr)
library(stringr)
library(pheatmap)
options(scipen=100)

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

projdir <- "/tscc/projects/PairedTF"

marker_regions <- read.table(file.path(projdir, "hic_pred/marker_regions.txt"), header=TRUE)
celltype_subclass <- read.table(file.path(projdir, "CTCF/celltype_subclass.txt"), header=TRUE, sep="\t")

# get real and predicted HiC matrix
hic_mtx_real10K_list <- list()
hic_mtx_cv10k_list <- list()

for (i in 1:nrow(marker_regions)) {
    chr <- marker_regions$chr[i]
    start <- marker_regions$region1[i]
    real_end <- as.integer(start) + 2100000
    pred_end <- as.integer(start) + 2097152
    gene <- marker_regions$gene[i]
  for (j in 1:nrow(celltype_subclass)) {
    subclass <- celltype_subclass$subclass1[j]
    celltype <- celltype_subclass$celltype1[j]
    mtx_dir <- file.path(projdir, "hic_pred/outputs/marker_regions", "mtx", celltype)
    mtx_real_10K_file <- paste0(mtx_dir, "/Real_", chr, "_", start, "_", real_end, "_res10K.txt")
    mtx_real_10K <- read.table(mtx_real_10K_file)
    hic_mtx_real10K_list[[gene]][[celltype]] <- as.numeric(as.matrix(mtx_real_10K[1:209,1:209])) 
    mtx_bam2bw <- file.path(projdir, "hic_pred/outputs/marker_regions", "mtx", paste0(celltype, "_bam2bw"))
    mtx_bam2bw_cv10k_file <- paste0(mtx_bam2bw, "/bam2bw_cv10K_", chr, "_", start, "_", pred_end, ".txt")
    mtx_bam2bw_cv10k <- read.table(mtx_bam2bw_cv10k_file)
    hic_mtx_cv10k_list[[gene]][[celltype]] <- as.numeric(as.matrix(mtx_bam2bw_cv10k[1:209,1:209]))
  }
}

hic_mtx_real10K_mtx <- list()
hic_mtx_cv10k_mtx <- list()
for (i in 1:nrow(marker_regions)) {
  gene <- marker_regions$gene[i]
  hic_mtx_real10K_mtx[[gene]] <- do.call(cbind, hic_mtx_real10K_list[[gene]])
  hic_mtx_cv10k_mtx[[gene]] <- do.call(cbind, hic_mtx_cv10k_list[[gene]])
}

# correlation between real and predicted HiC matrix
hic_mtx_realVScv10k_cor <- list()
for (i in 1:nrow(marker_regions)) {
  gene <- marker_regions$gene[i]
  hic_mtx_realVScv10k_cor[[gene]] <- cor(hic_mtx_real10K_mtx[[gene]], hic_mtx_cv10k_mtx[[gene]], method="spearman")
}

hic_mtx_realVScv10k_cor_zscore <- list()
for (i in 1:nrow(marker_regions)) {
  gene <- marker_regions$gene[i]
  hic_mtx_realVScv10k_cor_zscore[[gene]] <- scaleByRow(hic_mtx_realVScv10k_cor[[gene]], is.capped=FALSE)
}

# plot heatmap of mean correlation matrix
hic_mtx_realVScv10k_cor_zscore_mean <- Reduce('+', hic_mtx_realVScv10k_cor_zscore)
hic_mtx_realVScv10k_cor_zscore_mean_qn <- preprocessCore::normalize.quantiles(hic_mtx_realVScv10k_cor_zscore_mean)
rownames(hic_mtx_realVScv10k_cor_zscore_mean_qn) <- rownames(hic_mtx_realVScv10k_cor_zscore_mean)
colnames(hic_mtx_realVScv10k_cor_zscore_mean_qn) <- colnames(hic_mtx_realVScv10k_cor_zscore_mean)
pheatmap(hic_mtx_realVScv10k_cor_zscore_mean_qn, cluster_rows=FALSE, cluster_cols=FALSE)
pdf(file=file.path(projdir, "hic_pred/outputs/marker_regions/mtx/correlation/hic_mtx_realVScv10k_cor_zscore_mean_qn.pdf"), width=9, height=8)

