#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input RDS")
parser$add_argument("-n", "--name", default=NULL, help="rename cluster")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

input = as.character(args$input)
if (!is.null(args$name)) {
	name = as.character(args$name)
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

print("Step1: read rds data!")
RNAdata <- readRDS(input)

print("Step2: parse clustering meta!")
if (exists(quote(name))) {
	new.cluster.ids <- unlist(str_split(name, ","))
	names(new.cluster.ids) <- levels(RNAdata)
	RNAdata <- RenameIdents(RNAdata, new.cluster.ids)
	#RNA.clusters <- levels(Idents(RNAdata))
	cluster.table <- table(Idents(RNAdata))
	p_umap0 <- DimPlot(object=RNAdata, reduction="umap", label=TRUE) + NoLegend()
	print("Step3: save data!")
	write.csv(cluster.table, file=paste0(outF, "_reidents.table.csv"))
	RNAdata$celltype <- Idents(RNAdata)
	RNA_meta <- RNAdata@meta.data %>% tibble::rownames_to_column("barcode") %>% select(barcode,orig.ident,nCount_RNA,nFeature_RNA,percent.mt,seurat_clusters,sampleID,celltype) 
	write.csv(RNA_meta, file=paste0(outF, "_meta.table.csv"), row.names=FALSE, quote=FALSE)
	ggsave(p_umap0, file=paste0(outF, ".reidents.umap.pdf"), width=8, height=8)
	saveRDS(RNAdata, file=paste0(outF, ".reidents.rds"))
} else {
	stop("No reidents info is provided!")
}
print("Job is done!")

