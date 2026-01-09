#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input path")
parser$add_argument("-n", "--name", default=NULL, help="sample name")
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
suppressPackageStartupMessages(library("stringr"))

# Load one dataset
options(Seurat.object.assay.version = "v3")

load_1data <- function(folder, sample=NULL){
  if (is.null(sample)) {
    sample <- folder
  }
  names(folder) <- sample
  print(paste0("Load ", sample, " start."))
  data <- Read10X(folder)
  #Create Seurat object
  RNAdata <- CreateSeuratObject(counts=data, project=sample, min.cells=0, min.features=0)
  RNA_cells <- paste(sample, colnames(RNAdata), sep=":")
  RNAlist <- list(data=RNAdata, cell=RNA_cells)
  return(RNAlist)
}

print("Step1: load data!")
if (exists(quote(name))) {
  sample <- name
  RNAlist <- load_1data(input, sample)
} else {
  print("No sample names are provided, will use file names!")
  RNAlist <- load_1data(input)
}
RNAdata <- RNAlist$data
RNAcell <- RNAlist$cell

print("Step2: save data!")
saveRDS(RNAdata, file=paste0(outF, ".raw.rds"))
write.table(RNAcell, file=paste0(outF, ".cells_raw.txt"), row.names=FALSE, quote=FALSE, col.names=FALSE)
print("Job is done!")

