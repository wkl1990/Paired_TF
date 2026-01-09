#!/usr/bin/R

#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input prefix")
parser$add_argument("-t", "--type", default="name", help="type name")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

input = as.character(args$input)
name = as.character(args$type)

# run shell consensus matrix
plotConsensus <- function(consensusFile, type,
                          outPDF, width = 7, height = 7) {
  consensus <- read.table(consensusFile)
  colnames(consensus) <- c("file", "Resolution", "Dispersion",
                           "ProportionOfAmbiguousClustering")
  pdf(file = outPDF, width = width, height = height)
  par(mar = c(5, 4, 4, 4) + 0.3)
  plot(consensus$Resolution, consensus$Dispersion,
       type = "l", pch = 16, col = "red", lwd = 4,
       xlab = "Resolution", ylab = "Dispersion", cex.lab = 1.5,
       main =
         paste("Consensus analysis for Leiden-base clustering on",
               type, "cells")
       )
  par(new = TRUE)
  plot(consensus$Resolution, consensus$ProportionOfAmbiguousClustering,
       pch = 17, col = "blue",
       axes = FALSE, type = "l", lwd = 4, xlab = "", ylab = ""
       )
  axis(side = 4,
       at = pretty(range(consensus$ProportionOfAmbiguousClustering)))
  legend("topright", legend = c("Dispersion (left;red)", "PAC (right;blue)"),
         col = c("red", "blue"), cex = 1, lty = c(3, 3))
  mtext("Proportion Of Ambiguous Clustering (PAC)",
        side = 4, line = 3, cex = 1.5)
  dev.off()
}


plotConsensus(paste0(input, ".consensus.stat.txt"), type=name, outPDF=paste0(input, ".cluster.consensus.pdf"), width=12)





