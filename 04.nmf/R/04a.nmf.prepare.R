library(hdf5r)
library(tidyverse)
library(rlang)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
library(future)
library(data.table)

# read data
projdir <- "/tscc/projects/PairedTF"
cpm_pbysc <- read.csv(file.path(projdir, "TF.cpm_peakByCelltype.csv"), row.names=1)
cpm_pbysc_mat <- as.matrix(cpm_pbysc)
upValue <- quantile(cpm_pbysc_mat, 0.9999)
cpm_capped <- cpm_pbysc_mat
cpm_capped[cpm_capped > upValue] <- upValue
# to celltype by peak
cpm_capped <- t(cpm_capped)
peaks <- colnames(cpm_capped)
clusters <- rownames(cpm_capped)

# save the pmat with all the peak into a hdf5 format
conn <- hdf5r::H5File$new(file.path(projdir, "TF.cpm_peakByCelltype.h5"), mode = "w")
data.grp <- conn$create_group("X")
# NOTE: hdf5r will transpose the mat
# https://github.com/hhoeflin/hdf5r/issues/81
data.grp[["mat"]] <- cpm_capped
# colnames corresponds to cpm.capped
data.grp[["colnames"]] <- peaks
# rownames corresponds to cpm.capped
data.grp[["rownames"]] <- clusters
conn$close_all()

write.table(peaks,
  file = file.path(projdir, "TF.peak.nms.txt"),
  quote = FALSE,
  row.names = FALSE, col.names = FALSE)
write.table(clusters,
  file = file.path(projdir, "TF.celltype.nms.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE)

# set order of the clusters
subclass <- read.csv(file.path(projdir, "TF_subclass.csv"), sep=",", comment.char="#")
subclass %>% select(subclass_label, subclass_name) %>% distinct() %>% arrange(subclass_name) -> subclass_dict
RNA_meta <- read.csv(file.path(projdir, "TF_meta.table.csv"))
subclass$celltype <- RNA_meta$celltype[match(subclass$cell_id, RNA_meta$barcode)]
celltype2subclass <- gmodels::CrossTable(subclass$celltype, subclass$subclass_name)
celltype_map_subclass <- data.frame(table=rep(NA, length(unique(subclass$celltype))), row=rep(NA, length(unique(subclass$celltype))), col=rep(NA, length(unique(subclass$celltype))), tbl=rep(NA, length(unique(subclass$celltype))))
rownames(celltype_map_subclass) <- unique(subclass$celltype)
colnames(celltype_map_subclass) <- names(celltype2subclass)
for (celltype in rownames(celltype_map_subclass)) {
	for (stas in names(celltype2subclass)) {
		match_subclass <- colnames(celltype2subclass[[stas]])[which.max(celltype2subclass[[stas]][celltype,])]
		celltype_map_subclass[celltype, stas] <- match_subclass
	}
}
identical(celltype_map_subclass[,1], celltype_map_subclass[,4])

celltype_map_subclass %>% select(prop.tbl) %>% rownames_to_column("celltype") %>% arrange(prop.tbl, desc(celltype)) %>% pull(celltype) -> celltype_order
setequal(celltype_order, clusters)

write.table(celltype_order,
  file = file.path(projdir, "TF.celltype.srt.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE)


