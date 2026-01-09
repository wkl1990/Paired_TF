#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input RDS")
parser$add_argument("-a", "--assay", default="RNA", help="assay")
parser$add_argument("-g", "--organism", default="mouse", help="organism genome")
parser$add_argument("-s", "--sample", default="sampleID", help="sample name in meta table")
parser$add_argument("-l", "--list", default=NULL, help="another gene list to plot")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

input = as.character(args$input)
assay = as.character(args$assay)
organism = as.character(args$organism)
sampleID = as.character(args$sample)
if (!is.null(args$list)){
  list = as.character(args$list)
}
outF = as.character(args$output)

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("stringr"))

plot_feature_combine <- function(plot_umap=plot_umap, pt.size=0.1, object=NULL, features=NULL, reduction="umap"){
  height=ceiling(length(features) / 3)
  plot_marker <- FeaturePlot(object=object, features=features, pt.size=pt.size, max.cutoff="q95", ncol=3, reduction=reduction)
  plot_combine <- ggarrange(plot_umap, plot_marker, heights=c(1, height), nrow=2, ncol=1)
  return(plot_combine)
}

print("Step1: read rds data!")
RNAdata <- readRDS(input)

# brain features
ASAP_feature <- c("GFAP", "OLIG2", "CSPG4", "NES", "RBFOX3", "PTBP1", "PTBP2", "LGR5", "FGF17", "FLT1", "TSHZ2", "FOXP2") 
Epithelial_feature <- c("NES", "SOX2", "NOTCH1", "HES1", "HES3", "CDH1", "OCLN", "SOX10") #Neuroepithelial cells
RG_feature <- c("VIM", "NES", "PAX6", "HES1", "HES5", "GFAP", "SLC1A3", "FABP7", "TNC", "CDH2", "SOX2") #Radial glia
IP_feature <- c("EOMES", "ASCL1") #Intermediate progenitors
IMneurons_feature <- c("DCX", "TUBB3", "NEUROD1", "TBR1", "STMN1") #Immature neurons
OPC_feature <- c("PDGFRA", "CSPG4") #Oligodendrocyte precursor cells
MatureOligo_feature <- c("OLIG1", "OLIG2", "OLIG3", "MBP", "CLDN11", "MOG", "SOX10", "PLP1") #Mature oligodendrocytes
Schwann_feature <- c("MPZ", "NCAM1", "GAP43", "S100A1", "S100B", "NGFR") #Schwann cells
Astrocytes_feature <- c("GFAP", "SLC1A3", "SLC1A2", "GLUL", "S100B", "ALDH1L1") #Astrocytes
Microglia_feature <- c("TMEM119", "ITGAM", "PTPRC", "AIF1", "CX3CR1", "ADGRE1", "CD68", "CD40") #Microglia
MatureNeurons_feature <- c("RBFOX3", "MAP2", "NEFM", "NEFH", "SYP", "DLG4") #Mature neurons
GluNeurons_feature <- c("SLC17A7", "SLC17A6", "GRIN1", "GRIN2B", "GLS", "GLUL") #Glutamatergic neurons
GABANeurons_feature <- c("SLC6A1", "GABBR1", "GABBR2", "GAD2", "GAD1") #GABAergic neurons
DopaminNeurons_feature <- c("TH", "SLC6A3", "FOXA2", "KCNJ6", "NR4A2", "LMX1B") #Dopaminergic neurons
SerNeurons_feature  <- c("TPH1", "SLC6A4", "FEV") #Serotonergic neurons
ChNeurons_feature <- c("CHAT", "SLC18A3", "ACHE") #Cholinergic neurons

celltype_namelist <- c("ASAP features", "Neuroepithelial cells", "Radial glia", 
  "Intermediate progenitors", "Immature neurons", "Oligodendrocyte precursor cells", "Mature oligodendrocytes", 
  "Schwann cells", "Astrocytes", "Microglia", "Mature neurons", "Glutamatergic neurons", "GABAergic neurons",
  "Dopaminergic neurons", "Serotonergic neurons", "Cholinergic neurons")

feature_list <- list(ASAP_feature=ASAP_feature, Epithelial_feature=Epithelial_feature, RG_feature=RG_feature, 
  IP_feature=IP_feature, IMneurons_feature=IMneurons_feature, OPC_feature=OPC_feature, MatureOligo_feature=MatureOligo_feature, 
  Schwann_feature=Schwann_feature, Astrocytes_feature=Astrocytes_feature, Microglia_feature=Microglia_feature,
  MatureNeurons_feature=MatureNeurons_feature, GluNeurons_feature=GluNeurons_feature, GABANeurons_feature=GABANeurons_feature,
  DopaminNeurons_feature=DopaminNeurons_feature, SerNeurons_feature=SerNeurons_feature, ChNeurons_feature=ChNeurons_feature)

names(celltype_namelist) <- names(feature_list)

if (exists(quote(list))) {
  Users_feature <- as.character(unlist(str_split(list, ",")))
  celltype_namelist <- c(celltype_namelist, "Users features")
  feature_list[["Users_feature"]] <- Users_feature
} 

if (organism=="mouse"){
  for (i in names(feature_list)){
    feature_list[[i]] <- stringr::str_to_title(feature_list[[i]])
  }
}

print("Step2: feature plot!")
#st plot
DefaultAssay(RNAdata) <- "RNA"
p_umap0 <- DimPlot(object=RNAdata, reduction="umap", label=TRUE) + NoLegend()
p_umap1 <- DimPlot(object=RNAdata, reduction="umap", group.by="orig.ident")
p_umap2 <- DimPlot(object=RNAdata, reduction="umap", group.by=sampleID)
p_umap <- p_umap0 | p_umap1 | p_umap2

if (ncol(RNAdata)<100000) {
  pt_size <- 0.1
} else {
  pt_size <- NULL
}


plot_knownmarker <- list()
for (i in names(feature_list)) {
  if (length(intersect(feature_list[[i]],rownames(RNAdata)))==0){
    feature_list[[i]] <- NULL
    celltype_namelist <- setdiff(celltype_namelist, celltype_namelist[i])
    next
  }
  plot_knownmarker[[i]] <- plot_feature_combine(p_umap, pt_size, RNAdata, feature_list[[i]], "umap")
}

print("Step3: output plot!")
for (i in seq(length(feature_list))) {
  height=4 * (ceiling(length(feature_list[[i]]) / 3) + 1)
  ggsave(file=paste0(outF, ".", names(feature_list)[i], ".pdf"), plot=gridExtra::arrangeGrob(grobs=plot_knownmarker[i], bottom=celltype_namelist[i]), width=12, height=height)
}
print("Job is done!")


