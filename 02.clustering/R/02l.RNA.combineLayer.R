#!/usr/bin/env Rscript
# example: Rscript /tscc/scripts/R/02l.RNA.combineLayer.R -i PairedTF.merge.combined.sampleID.rds -s sampleID -o cluster/PairedTF.merge &> log/PairedTF.merge.log

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input RDS")
parser$add_argument("-s", "--sample", default="sampleID", help="sample name in meta table")
parser$add_argument("-n", "--nFeature", default=2000, help="number of variable features")
parser$add_argument("-c", "--npc", default=100, help="number of pcs")
parser$add_argument("-m", "--method", default=4, help="clustering method, default Leiden")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

input = as.character(args$input)
sampleID = as.character(args$sample)
nFeature = as.numeric(args$nFeature)
npc = as.numeric(args$npc)
method = as.numeric(args$method)
outF = as.character(args$output)

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("findPC"))


autoPC_rna <- function(data, npc=100, assay="RNA", reduction="pca", method="first derivative") {
  DefaultAssay(data) <- assay
  data.use <- Stdev(data, reduction=reduction)
  pc_len <- length(data.use)
  npc <- ifelse(pc_len>=max(npc), npc, unique(ceiling(seq(pc_len, pc_len/3, length.out=length(npc)))))
  pc <- as.numeric(findPC(sdev=data.use, number=npc, method=method))
  return(pc)
}

cal_silhouette <- function(x, cluster){
  if(length(unique(cluster)) < 2) {
    message("Only one cluster level were found!")
    return(NA)
  }
  #diss <- parDist(x)
  #sil <- cluster::silhouette(cluster, diss)
  #ss <- mean(sil[,3])
  sil_score <- clusterCrit::intCriteria(x, cluster, "Silhouette")
  ss <- sil_score$silhouette
  return(ss)  
}

get_silhouette <- function(data, assay="RNA", matrix="pca"){
  DefaultAssay(data) <- assay
  if (matrix=="pca") {
    mat <- Embeddings(data, reduction="pca")
  } else {
    data_use <- GetAssayData(data, slot="scale.data")
    mat <- t(data_use)
  }
  resolutions <- seq(0.1,2,0.1)
  silhouette_tbl <- data.frame(Resolution=resolutions, Sil_score=rep(NA,length(resolutions)))
  for (resolution in resolutions) {
    message("Calculate silhouette score for resolution ", resolution)
    data <- cluster_rna(data, resolution=resolution, algorithm=4, assay=assay)
    clusters <- as.integer(data$seurat_clusters)
    sil_score <- cal_silhouette(mat, clusters)
    silhouette_tbl$Sil_score[which(silhouette_tbl$Resolution==resolution)] <- sil_score
  }
  silhouette_select <- silhouette_tbl$Resolution[which.max(silhouette_tbl$Sil_score)]
  return(silhouette_select)
}

cluster_rna <- function(RNAdata, resolution=resolution, algorithm=method, assay="integrated") {
  tryCatch(
    {
      DefaultAssay(RNAdata) <- assay
      RNAdata <- FindClusters(object=RNAdata, resolution=resolution, algorithm=algorithm)
      return(RNAdata)
    }, error=function(error_message) {
      message(error_message)
      DefaultAssay(RNAdata) <- assay
      RNAdata <- FindClusters(object=RNAdata, resolution=resolution, algorithm=1)
      return(RNAdata)
    }
  )
}


print("Step1: read rds data!")
RNAdata <- readRDS(input)
if (!is.null(RNAdata@meta.data[[sampleID]])) {
  RNAdata[["RNA"]] <- split(RNAdata[["RNA"]], f=RNAdata@meta.data[[sampleID]])
} else {
  stop(message(sampleID, " is not exist!"))
}

print("Step2: run seurat5 cca!")
RNAdata <- NormalizeData(RNAdata, verbose=FALSE)
RNAdata <- FindVariableFeatures(RNAdata, selection.method="vst", nfeatures=nFeature, verbose=FALSE)
RNAdata <- ScaleData(RNAdata, verbose=FALSE)
RNAdata <- RunPCA(RNAdata, npcs=npc, verbose=FALSE)

RNAdata <- IntegrateLayers(object=RNAdata, method=CCAIntegration, orig.reduction="pca", new.reduction="integrated.cca", verbose=FALSE)
RNAdata[["RNA"]] <- JoinLayers(RNAdata[["RNA"]])

n_pc <- autoPC_rna(RNAdata, assay="RNA")
RNAdata <- FindNeighbors(RNAdata, reduction="integrated.cca", dims=1:n_pc, verbose=FALSE)

sil_resolution <- get_silhouette(RNAdata, assay="RNA")
RNAdata <- cluster_rna(RNAdata, resolution=sil_resolution, algorithm=method, assay="RNA")
RNAdata <- RunUMAP(RNAdata, dims=1:n_pc, a=1.8956, b=0.8006, reduction="integrated.cca", verbose=FALSE)

print("Step3: save data!")
pt_seurat5 <- DimPlot(RNAdata, reduction="umap", group.by="seurat_clusters") + NoLegend()
ggsave(pt_seurat5, file=paste0(outF, ".seurat5.cluster.umap.pdf"), width=16, height=8)
saveRDS(RNAdata, file=paste0(outF, ".seurat5.cluster.rds"))

print("Job is done!")




