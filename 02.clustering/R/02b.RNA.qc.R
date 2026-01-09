#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input raw rds")
parser$add_argument("-n", "--nFeature", default=500, help="nFeature cutoff")
parser$add_argument("-u", "--uFeature", default=NULL, help="max nFeature cutoff")
parser$add_argument("-m", "--mt", default=5, help="mitochondria cutoff")
parser$add_argument("-c", "--mincells", default=1, help="genes detected in at least the cutoff cells")
parser$add_argument("-g", "--organism", default="mouse", help="organism species")
parser$add_argument("-t", "--plot", action='store_true', help="whether to output plot")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

input = as.character(args$input)
organism = as.character(args$organism)
nFeature = as.numeric(args$nFeature)
if (!is.null(args$uFeature)){
  uFeature = as.numeric(args$uFeature)
}
mt = as.numeric(args$mt)
mincells = as.numeric(args$mincells)
plot = ifelse(args$plot, TRUE, FALSE) 
outF = as.character(args$output)

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))


qc_rna <- function(RNAdata, organism="mouse", nFeature=500, uFeature=NULL, mt=5, min.cells=1) {
  if (organism=="mouse") {
    pattern <- "^mt-"
  } else if (organism=="human") {
    pattern <- "^MT-"
  } else {
    stop("Unrecognized organism!")
  }
  RNAdata[["percent.mt"]] <- PercentageFeatureSet(RNAdata, pattern=pattern)
  RNAdata <- subset(RNAdata, subset=nFeature_RNA>nFeature & percent.mt<mt)
  if (!is.null(uFeature)) {RNAdata <- subset(RNAdata, subset=nFeature_RNA<uFeature)}
  RNAdata_counts <- GetAssayData(RNAdata, slot="count")
  RNAdata_genes <- which(rowSums(RNAdata_counts>0)>=min.cells)
  #print(length(RNAdata_genes))
  RNAdata <- RNAdata[RNAdata_genes,]    
  pt_qc_vln <- VlnPlot(RNAdata, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
  pt_qc_UMIvsMT <- FeatureScatter(RNAdata, feature1="nFeature_RNA", feature2="percent.mt")
  qclist <- list(data=RNAdata, vlnplot=pt_qc_vln, scatplot=pt_qc_UMIvsMT)
  return(qclist)
}

print("Step1: read rds data!")
RNAdata <- readRDS(input)

# QC
print("Step2: qc start!")
if (length(RNAdata)==1) {
  if (exists(quote(uFeature))) {
    maxFeature <- uFeature
  } else {
    maxFeature <- NULL
  }
#  maxFeature <- ifelse((exists(quote(uFeature))), uFeature, NULL)
  RNAdata_qclist <- qc_rna(RNAdata, organism=organism, nFeature=nFeature, uFeature=maxFeature, mt=mt, min.cells=mincells)
  RNAsce <- RNAdata_qclist$data
  pt_qc_vln <- RNAdata_qclist$vlnplot
  pt_qc_UMIvsMT <- RNAdata_qclist$scatplot
  RNA_cells_qc <- paste(names(RNAsce), colnames(RNAsce), sep=":")
} else {
  stop("Error: do not support more than one sample currently!")
}

print("Step3: save data!")
if (plot==TRUE){
  ggsave(pt_qc_vln, file=paste0(outF, ".qc.vln.pdf"), width=16, height=8)
  ggsave(pt_qc_UMIvsMT, file=paste0(outF, ".qc.UMIvsMT.pdf"), width=8, height=8)
}
saveRDS(RNAsce, file=paste0(outF, ".qc.rds"))
write.table(RNA_cells_qc, file=paste0(outF, ".cells_qc.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
print("Job is done!")




