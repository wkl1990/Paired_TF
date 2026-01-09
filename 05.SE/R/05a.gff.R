library(dplyr)
library(stringr)
library(tidyr)

# prepare peak gff for ROSE
projdir <- "/tscc/projects/"
RNAPII_peaks <- read.csv(file.path(projdir, "RNAPII.peaks.union.csv"))
RNAPII_peaks %>% tidyr::separate(col=Peaks, into=c("chr", "start", "end"), sep="[:-]", remove=FALSE) %>% filter(chr != "chrM") -> RNAPII_peaks
for (celltype in colnames(RNAPII_peaks)[-c(1:4)]) {
  print(celltype)
  peaks_gff <- NULL
  RNAPII_peaks %>% select(Peaks, chr, start, end, celltype) %>% filter(!!as.symbol(celltype)=="True") %>% mutate(strand=".", Peaks2=Peaks, Peaks3=Peaks, score=0, frame=0) %>% select(chr, Peaks, Peaks2, start, end, score, strand, frame, Peaks3) -> peaks_gff
  write.table(peaks_gff, file=file.path(projdir, paste0("peak_gff/", celltype, ".gff")), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
}
