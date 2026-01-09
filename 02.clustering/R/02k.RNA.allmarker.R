#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input RDS")
parser$add_argument("-a", "--assay", default="RNA", help="assay")
parser$add_argument("-p", "--min_pct", default=0.5, help="min.pct")
parser$add_argument("-f", "--logfc_threshold", default=1, help="logfc.threshold")
parser$add_argument("-d", "--min_diff_pct", default=0.3, help="min.diff.pct")
parser$add_argument("-t", "--return_thresh", default=0.01, help="return.thresh")
parser$add_argument("-n", "--name", default=NULL, help="rename cluster")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

input = as.character(args$input)
assay = as.character(args$assay)
min_pct = as.numeric(args$min_pct)
logfc_threshold = as.numeric(args$logfc_threshold) 
min_diff_pct = as.numeric(args$min_diff_pct)
return_thresh = as.numeric(args$return_thresh)
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
DefaultAssay(RNAdata) <- assay
geneMarkerV2 <- read.csv("/tscc/projects/fromCEMBA/geneMarkerV2.csv")


print("Step2: parse clustering meta!")
if (exists(quote(name))) {
	new.cluster.ids <- unlist(str_split(name, ","))
	names(new.cluster.ids) <- levels(RNAdata)
	RNAdata <- RenameIdents(RNAdata, new.cluster.ids)
	RNA.clusters <- levels(Idents(RNAdata))
	cluster.table <- table(Idents(RNAdata))
	write.csv(cluster.table, file=paste0(outF, "_recluster.table.csv"))
	p_umap0 <- DimPlot(object=RNAdata, reduction="umap", label=TRUE) + NoLegend()
	ggsave(p_umap0, file=paste0("./plot/", outF, ".recluster.umap.pdf"), width=8, height=8)
	saveRDS(RNAdata, file=paste0(outF, ".recluster.rds"))
} else {
	print("No recluster info is provided. Use original cluster info!")
	RNA.clusters <- levels(Idents(RNAdata))
}


print("Step3: find all markers!")
if (length(RNA.clusters)>1) {
	RNA.positive.markers <- FindAllMarkers(RNAdata, only.pos=TRUE, min.pct=min_pct, logfc.threshold=logfc_threshold, min.diff.pct=min_diff_pct, min.cells.feature=10, return.thresh=return_thresh, verbose=FALSE)
	RNA.positive.all.markers <- RNA.positive.markers %>% mutate(pct.diff=pct.1-pct.2, pct.diff.FC=(pct.1-pct.2)/pct.1) %>% select(gene, everything()) %>% subset(p_val_adj<0.01)
	RNA.positive.all.markers$celltype <- geneMarkerV2$celltype[match(RNA.positive.all.markers$gene, geneMarkerV2$gene)]
    RNA.positive.all.markers$subclass <- geneMarkerV2$subclass[match(RNA.positive.all.markers$gene, geneMarkerV2$gene)]
    RNA.positive.all.markers$class <- geneMarkerV2$class[match(RNA.positive.all.markers$gene, geneMarkerV2$gene)]
    RNA.positive.all.markers$PMID <- geneMarkerV2$PMID[match(RNA.positive.all.markers$gene, geneMarkerV2$gene)]
	positive.top10 <- RNA.positive.all.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
	RNA.positive.all.markers.table <- table(RNA.positive.all.markers$cluster)
	print("Step4: save data!")
	saveRDS(RNA.positive.markers, file=paste0(outF, ".positive.markers.rds"))
	write.csv(RNA.positive.all.markers, file=paste0(outF, ".positive.all.markers.csv"), row.names=FALSE)
	write.csv(positive.top10, file=paste0(outF, ".positive.top10.markers.csv"), row.names=FALSE)
	write.csv(RNA.positive.all.markers.table, file=paste0(outF, ".positive.all.markers.table.csv"), row.names=FALSE)	
} else {
	print("Only 1 cluster in data. No marker can be get!")
}
print("Job is done!")


