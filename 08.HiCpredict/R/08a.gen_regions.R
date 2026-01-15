library(dplyr)
library(stringr)
library(bedtoolsr)

# generate marker gene regions for hic predicition
projdir <- "/tscc/projects/PairedTF"
marker_bed <- read.table(file.path(projdir, "hic_pred/marker_genes.bed"))
marker_bed %>% dplyr::rename("chr"="V1", "start"="V2", "end"="V3", "gene"="V4") %>% mutate(len=end-start, upstream=start-5000, region1=floor(upstream/10000)*10000, region2=region1+2000000) -> marker_bed
# check gaps
marker_bed %>% select(chr, region1, region2) %>% dplyr::rename("chrom"="chr", "start"="region1", "end"="region2") -> marker_pred_bed
gap_bed <- read.table(file.path(projdir, "/data/mm10/centrotelo.bed"))
gap_bed %>% dplyr::rename("chrom"="V1", "start"="V2", "end"="V3") -> gap_bed
marker_pred_gap <- bedtoolsr::bt.intersect(marker_pred_bed, gap_bed, wo=TRUE)  
write.table(marker_bed, file=file.path(projdir, "hic_pred/marker_regions.txt"), quote=FALSE, row.names=FALSE, sep="\t")
