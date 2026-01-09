#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input RDS")
parser$add_argument("-a", "--assay", default="RNA", help="assay")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

input = as.character(args$input)
assay = as.character(args$assay)
outF = as.character(args$output)

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("stringr"))

print("Step1: read rds data!")
RNAdata <- readRDS(input)
DefaultAssay(RNAdata) <- assay
geneMarkerV2 <- read.csv("/tscc/projects/fromCEMBA/geneMarkerV2.csv")

#marker gene
print("Step2: run pairwise markers!")
RNA.clusters <- levels(Idents(RNAdata))

if (length(RNA.clusters)>1) {
 cluster.markers <- list()
  marker.table <- c("cluster1", "cluster2", "marker_num")
  for (i in RNA.clusters) {
    cluster1.markers <- list()
    for (j in RNA.clusters[which(RNA.clusters != i)]) {
      markers <- FindMarkers(RNAdata, ident.1=i, ident.2=j, only.pos=TRUE, min.pct=0.5, logfc.threshold=1, min.diff.pct=0.3, min.cells.feature=10, return.thresh=0.01, verbose=FALSE)
      if (nrow(markers) != 0) {
        markers$cluster1 <- i 
        markers$cluster2 <- j
        marker_num <- nrow(markers)
        marker.table <- rbind(marker.table, c(i, j, marker_num))
        cluster1.markers[[paste(i, j, sep="_")]] <- markers %>% tibble::rownames_to_column("gene") %>% mutate(pct.diff=pct.1-pct.2, pct.diff.FC=(pct.1-pct.2)/pct.1) %>% select(gene, everything()) %>% subset(p_val_adj<0.01)       
      }
      rm(marker_num, markers)  
    }
    cluster.markers[[i]] <- cluster1.markers
    rm(cluster1.markers)
  }
  if (!is.null(dim(marker.table))) {
    colnames(marker.table) <- marker.table[1,]
    marker.table1 <- matrix(as.numeric(marker.table[-1,]), ncol=ncol(marker.table))
    colnames(marker.table1) <- colnames(marker.table)
    marker.gene <- do.call(rbind, lapply(cluster.markers, function(x) do.call(rbind, x))) 
    marker.gene$celltype <- geneMarkerV2$celltype[match(marker.gene$gene, geneMarkerV2$gene)]
    marker.gene$subclass <- geneMarkerV2$subclass[match(marker.gene$gene, geneMarkerV2$gene)]
    marker.gene$class <- geneMarkerV2$class[match(marker.gene$gene, geneMarkerV2$gene)]
    marker.gene$PMID <- geneMarkerV2$PMID[match(marker.gene$gene, geneMarkerV2$gene)]
    write.csv(marker.table1, file=paste0(outF, "_marker.table.csv"))
    saveRDS(cluster.markers, file=paste0(outF, "_cluster.markers.rds"))
    write.csv(marker.gene, file=paste0(outF, "_pairmarker.gene.csv"))
  }
} else {
  print(paste0("Warning:", outF, " has only 1 cluster!"))
}

print("Step3: cluster merge statistics!")
reduce_redundant <- function(x){
  a <- as.numeric(str_split(x, "_", simplify=TRUE)[1,1])
  b <- as.numeric(str_split(x, "_", simplify=TRUE)[1,2])
  pairs <- ifelse(a<b, paste0(a, "_", b), paste0(b, "_", a))  
  return(pairs)
}

if (exists(quote(marker.table1))) {
  all_pairs <- c()
  for (i in RNA.clusters) {
    for (j in RNA.clusters[which(RNA.clusters != i)]) {
      all_pairs <- c(all_pairs, paste0(i, "_", j))
    }
  }
  seq_pairs <- paste(marker.table1[,"cluster1"], marker.table1[,"cluster2"], sep="_")
  merge_pairs <- setdiff(all_pairs, seq_pairs)
  if (length(merge_pairs) != 0){
    write.csv(merge_pairs, file=paste0(outF, "_cluster.merge_pairs.csv"))
    merge_clusters <- str_split(merge_pairs, "_", simplify=TRUE)
    merge_clusters_rev <- paste0(merge_clusters[,2], "_", merge_clusters[,1])
    merge_pairs_share <- intersect(merge_pairs, merge_clusters_rev)
    write.csv(merge_pairs_share, file=paste0(outF, "_cluster.merge_pairs_share.csv"))
    merge_pairs_share_reduce <- unique(unlist(lapply(merge_pairs_share, reduce_redundant)))
    write.csv(merge_pairs_share_reduce, file=paste0(outF, "_cluster.merge_pairs_share_reduce.csv"))
    merge_pairs_seq <- setdiff(merge_pairs, merge_clusters_rev)
    write.csv(merge_pairs_seq, file=paste0(outF, "_cluster.merge_pairs_seq.csv"))
  }
}
print("Job is done!")


