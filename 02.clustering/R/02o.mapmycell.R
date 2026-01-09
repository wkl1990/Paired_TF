#!/usr/bin/env Rscript

library(Seurat)
library(data.table)
library(dplyr)
library(stringr)
library(ComplexHeatmap)

# get h5ad for MapMyCells (ABC atlas)
PairedTF_RNA <- readRDS("21.seurat/cluster/PairedTF.cluster.rds")
PairedTF_counts <- PairedTF_RNA@assays$RNA$counts
PairedTF_counts_trans <- t(as.matrix(PairedTF_counts))
PairedTF_sparse_counts <- as(PairedTF_counts_trans, "dgCMatrix")
library(anndata) # Load library
countAD <- AnnData(X=PairedTF_sparse_counts)
object.size(countAD)
write_h5ad(countAD, "21.seurat/cluster/PairedTF.counts.h5ad")

# MapMyCells online 
#https://knowledge.brain-map.org/mapmycells/process 

# check MapMyCells
RNA_meta <- read.csv("21.seurat/cluster/PairedTF_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping.csv", sep=",", comment.char="#")
if (identical(RNA_meta$cell_id, colnames(PairedTF_RNA))) {PairedTF_meta <- PairedTF_RNA@meta.data}
PairedTF_meta %>% select(orig.ident, nCount_RNA, nFeature_RNA, percent.mt, seurat_clusters, sampleID) -> PairedTF_meta
if (identical(rownames(PairedTF_meta), RNA_meta$cell_id)) {
  PairedTF_meta$subclass <- RNA_meta$subclass_name
  PairedTF_RNA$subclass <- RNA_meta$subclass_name
}

table(PairedTF_meta$seurat_clusters, PairedTF_meta$subclass)

PairedTF_seurat2subclass <- gmodels::CrossTable(PairedTF_meta$seurat_clusters, PairedTF_meta$subclass)

PairedTF_seurat2subclass_row <- PairedTF_seurat2subclass$prop.row
PairedTF_seurat2subclass_row <- PairedTF_seurat2subclass_row[,which(apply(PairedTF_seurat2subclass_row, 2, max)>=0.01)]
col_fun <- circlize::colorRamp2(c(0, 1), c("white", "red"))
cols_order <- union(apply(PairedTF_seurat2subclass_row, 1, which.max), c(1:ncol(PairedTF_seurat2subclass_row)))

pt_seurat2subclass_row <- Heatmap(PairedTF_seurat2subclass_row, name="prop.row", col=col_fun, 
  column_order=cols_order, 
  row_order=c(1:nrow(PairedTF_seurat2subclass_row)),
  #cluster_rows=hclust_rows,
  show_row_names=TRUE, 
  row_names_gp=gpar(fontsize=8),
  show_column_names=TRUE,
  column_names_gp=gpar(fontsize=8),
  column_title = "subclass",
  row_title = "cluster")

PairedTF_seurat2subclass_col <- PairedTF_seurat2subclass$prop.col
PairedTF_seurat2subclass_col <- PairedTF_seurat2subclass_col[,which(apply(PairedTF_seurat2subclass_col, 2, max)>=0.01)]

col_fun <- circlize::colorRamp2(c(0, 1), c("white", "red"))
cols_order <- union(apply(PairedTF_seurat2subclass_col, 1, which.max), c(1:ncol(PairedTF_seurat2subclass_col)))

pt_seurat2subclass_col <- Heatmap(PairedTF_seurat2subclass_col, name="prop.col", col=col_fun, 
  column_order=cols_order, 
  row_order=c(1:nrow(PairedTF_seurat2subclass_col)),
  #cluster_rows=hclust_rows,
  show_row_names=TRUE, 
  row_names_gp=gpar(fontsize=8),
  show_column_names=TRUE,
  column_names_gp=gpar(fontsize=8),
  column_title = "subclass",
  row_title = "cluster")

pdf(file="21.seurat/cluster/PairedTF_seurat2subclass.heatmap.pdf", width=8, height=8)
pt_seurat2subclass_row
pt_seurat2subclass_col
dev.off()
