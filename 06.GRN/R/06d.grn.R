# generate cell type proximal peaks
projdir <- "/tscc/projects/PairedTF"
TFs <- read.table(file.path(projdir, "selected_TFs.txt"))$V1
tf_list <- read.table(file.path(projdir, "reference/tf_lists/allTFs_mm.txt"))

TF_RNA_combined <- readRDS(file=file.path(projdir, "TF.final.reidents.rds"))
celltypes <- levels(TF_RNA_combined)
for (TF in TFs) {
  peak_file <- file.path(projdir, paste0(TF, ".peaks.union.csv"))
  proximal_file <- file.path(projdir, paste0(TF, ".peaks.union.tssUpDn1k.proximal.distal.peaks"))
  if (file.exists(peak_file)&file.exists(proximal_file)){
    message(TF)
    dir.create(file.path(projdir, "celltype"), recursive=TRUE, showWarnings=FALSE)
    TF_peaks <- read.csv(peak_file)
    TF_peaks %>% tidyr::separate(Peaks, into=c("chr", "start", "end"), sep="[:-]", remove=FALSE) %>% mutate(Peaks=paste0(chr, "_", start, "_", end)) %>% filter(chr != "chrM") -> TF_peaks
    proximal_peaks <- read.table(proximal_file, header=FALSE)
    for (celltype in celltypes) {
      TF_peaks %>% filter(!!sym(celltype) == "True") %>% select(chr, start, end) -> celltype_peaks
      write.table(celltype_peaks, file=file.path(projdir, "celltype", paste0(TF, ".", celltype, ".peaks.bed")), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
      celltype_peaks %>% mutate(Peaks=paste0(chr, "_", start, "_", end)) -> celltype_peaks
      proximal_peaks %>% filter(V2 == "proximal") %>% filter(V1 %in% celltype_peaks$Peaks) -> celltype_proximal_peaks
      write.table(celltype_proximal_peaks, file=file.path(projdir, "celltype", paste0(TF, ".", celltype, ".proximal.peaks")), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
    }
  }
}

# filter TF in PGC and proximal
for (tf_name in TFs) {
  message(tf_name)
  for (celltype in celltypes) {
    message(celltype)
    celltype_proximal_file <- file.path(projdir, "celltype", paste0(tf_name, ".", celltype, ".proximal.peaks"))
    if (file.exists(celltype_proximal_file)&file.info(celltype_proximal_file)$size>0) {
      celltype_proximal_peaks <- read.table(celltype_proximal_file)
      celltype_proximal_peaks %>% filter(V3 %in% tf_list$V1) -> celltype_proximal_TFs
      message("Total proximal peaks: ", nrow(celltype_proximal_peaks))
      message("Proximal peaks with TFs: ", nrow(celltype_proximal_TFs))
      write.table(celltype_proximal_TFs, file=file.path(projdir, "celltype", paste0(tf_name, ".", celltype, ".proximal.TFs")), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
    } else {
      message("No proximal peaks found for ", celltype)
      write.table(NULL, file=file.path(projdir, "celltype", paste0(tf_name, ".", celltype, ".proximal.TFs")), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
    }
  }
}

for (tf_name in TFs) {
  message(tf_name)
  TF_proximal_gene_list <- list()
  for (celltype_name in celltypes) {
    message(celltype_name)
    celltype_proximal_file <- file.path(projdir, "celltype", paste0(tf_name, ".", celltype_name, ".proximal.peaks"))
    if (file.exists(celltype_proximal_file)&file.info(celltype_proximal_file)$size>0) {
      celltype_proximal_peaks <- read.table(celltype_proximal_file)
      message("Total proximal peaks: ", nrow(celltype_proximal_peaks))
      TF_proximal_genes <- unique(celltype_proximal_peaks$V3)
      message("Total proximal genes: ", length(TF_proximal_genes))
      TF_proximal_gene_list[[celltype_name]] <- data.frame(TF=tf_name, celltype=celltype_name, gene=TF_proximal_genes)
    } else if (file.exists(celltype_proximal_file)&file.info(celltype_proximal_file)$size==0) {
      message("No proximal peaks found for ", celltype_name)
      TF_proximal_gene_list[[celltype_name]] <- data.frame(TF=NULL, celltype=NULL, gene=NULL)
    } else {
      message("No proximal peaks file found for ", celltype_name)
      TF_proximal_gene_list[[celltype_name]] <- data.frame(TF=NULL, celltype=NULL, gene=NULL)
    }
  }

  do.call(rbind, TF_proximal_gene_list) -> TF_proximal_gene_all
  write.csv(TF_proximal_gene_all, file=file.path(projdir, paste0(tf_name, ".celltypes.proximal.genes.csv")), row.names=FALSE, quote=FALSE)
  TF_proximal_gene_num <- data.frame(celltype=factor(names(TF_proximal_gene_list), levels=rev(names(TF_proximal_gene_list))), number=unlist(lapply(TF_proximal_gene_list,nrow)))
  write.csv(TF_proximal_gene_num, file=file.path(projdir, paste0(tf_name, ".celltypes.proximal.genes.num.csv")), row.names=FALSE, quote=FALSE)
 }

for (tf_name in TFs) {
  message(tf_name)
  TF_proximal_TFs_list <- list()
  for (celltype_name in celltypes) {
    message(celltype_name)
    celltype_proximal_file <- file.path(projdir, "celltype", paste0(tf_name, ".", celltype_name, ".proximal.TFs"))
    if (file.exists(celltype_proximal_file)&file.info(celltype_proximal_file)$size>0) {
      celltype_proximal_peaks <- read.table(celltype_proximal_file)
      message("Total proximal peaks: ", nrow(celltype_proximal_peaks))
      TF_proximal_TFs <- unique(celltype_proximal_peaks$V3)
      message("Total proximal TFs: ", length(TF_proximal_TFs))
      TF_proximal_TFs_list[[celltype_name]] <- data.frame(TF=tf_name, celltype=celltype_name, gene=TF_proximal_TFs)
    } else if (file.exists(celltype_proximal_file)&file.info(celltype_proximal_file)$size==0) {
      message("No proximal peaks found for ", celltype_name)
      TF_proximal_TFs_list[[celltype_name]] <- data.frame(TF=NULL, celltype=NULL, gene=NULL)
    } else {
      message("No proximal peaks file found for ", celltype_name)
      TF_proximal_TFs_list[[celltype_name]] <- data.frame(TF=NULL, celltype=NULL, gene=NULL)
    }
  }

  do.call(rbind, TF_proximal_TFs_list) -> TF_proximal_TFs_all
  write.csv(TF_proximal_TFs_all, file=file.path(projdir, paste0(tf_name, ".celltypes.proximal.TFs.csv")), row.names=FALSE, quote=FALSE)
  TF_proximal_TFs_num <- data.frame(celltype=factor(names(TF_proximal_TFs_list), levels=rev(names(TF_proximal_TFs_list))), number=unlist(lapply(TF_proximal_TFs_list,nrow)))
  write.csv(TF_proximal_TFs_num, file=file.path(projdir, paste0(tf_name, ".celltypes.proximal.TFs.num.csv")), row.names=FALSE, quote=FALSE)
}

for (tf_name in TFs) {
  message(tf_name)
  allpgc_file <- file.path(projdir, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.tsv"))
  pospgc_file <- file.path(projdir, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".pos.tsv"))
  negpgc_file <- file.path(projdir, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".neg.tsv"))
  all_pgc <- read.table(allpgc_file, header=TRUE, sep="\t")
  pos_pgc <- read.table(pospgc_file, header=TRUE, sep="\t")
  neg_pgc <- read.table(negpgc_file, header=TRUE, sep="\t")
  message("Total PGC pairs: ", nrow(all_pgc), "\t", "pos pairs: ", nrow(pos_pgc), "\t", "neg pairs: ", nrow(neg_pgc))
  all_pgc %>% filter(gene %in% tf_list$V1) -> all_pgc_TFs
  pos_pgc %>% filter(gene %in% tf_list$V1) -> pos_pgc_TFs
  neg_pgc %>% filter(gene %in% tf_list$V1) -> neg_pgc_TFs
  message("Total PGC pairs with TFs: ", nrow(all_pgc_TFs), "\t", "pos pairs with TFs: ", nrow(pos_pgc_TFs), "\t", "neg pairs with TFs: ", nrow(neg_pgc_TFs))
  write.table(all_pgc_TFs, file=file.path(projdir, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.TFs")), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  write.table(pos_pgc_TFs, file=file.path(projdir, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".pos.TFs")), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  write.table(neg_pgc_TFs, file=file.path(projdir, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".neg.TFs")), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
}


for (tf_name in TFs) {
  message(tf_name)
  allpgc_file <- file.path(projdir, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.tsv"))
  all_pgc <- read.table(allpgc_file, header=TRUE, sep="\t")

  # filter cell type peaks
  TF_union_peak <- read.csv(file.path(projdir, "results/", paste0(tf_name, ".peaks.union.csv")))
  all_pgc %>% left_join(TF_union_peak, by=join_by(peak==Peaks)) -> all_pgc_celltype
  all_pgc_celltype %>% tidyr::separate(peak, into=c("chrom", "start", "end"), sep="[:-]", remove=FALSE) %>% relocate(chrom, start, end) -> all_pgc_celltype_bed

  all_pgc_list <- list()
  for (celltype in celltypes) {
  #  celltype <- celltype2subclass$celltype1[i]
  #  subclass <- celltype2subclass$subclass_label[i]
  #  subclass_name <- gsub(" ", "_", gsub("/", "-", subclass))
    message(celltype)
  #  ATAC_PGC_file <- paste0("/projects/ps-renlab2/kaw033/JX_PairedTag/pub_data/ATAC/wmb_catlas_ahchu/sa2pdc_bedpe/", "sa2subclass.", subclass_name, ".pdc.bedpe")
  #  if (!file.exists(ATAC_PGC_file)) {
  #    message(paste0("ATAC file not found: ", celltype))
  #  }
    all_pgc_celltype %>% filter(!!as.symbol(celltype)=="True") %>% dplyr::select(conns, pcc, class, Pval, FDR, gene, peak) %>% mutate(celltype=celltype) -> all_pgc_filter
    all_pgc_list[[celltype]] <- all_pgc_filter
    message("Total PGC pairs in ", celltype, ": ", nrow(all_pgc_filter))
    write.table(all_pgc_filter, file=file.path(projdir, "/results/celltype", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.", celltype, ".tsv")), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  }

  do.call(rbind, all_pgc_list) -> all_pgc_all
  write.csv(all_pgc_all, file=file.path(projdir, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.celltypes.csv")), quote=FALSE, row.names=FALSE)
  all_pgc_num <- data.frame(celltype=factor(names(all_pgc_list), levels=rev(names(all_pgc_list))), number=unlist(lapply(all_pgc_list,nrow)))
  write.csv(all_pgc_num, file=file.path(projdir, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.celltypes.num.csv")), row.names=FALSE, quote=FALSE)
}


for (tf_name in TFs) {
  message(tf_name)
  TF_pgc_TFs_list <- list()
  for (celltype_name in celltypes) {
    message(celltype_name)
    celltype_all_pgc_file <- file.path(projdir, "/results/celltype", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.", celltype_name, ".tsv"))
    if (file.exists(celltype_all_pgc_file)&file.info(celltype_all_pgc_file)$size>0) {
      celltype_all_pgc <- read.table(celltype_all_pgc_file, header=TRUE, sep="\t")
      message("Total peak and gene pairs: ", nrow(celltype_all_pgc))
      celltype_all_pgc %>% filter(gene %in% tf_list$V1) -> celltype_TF_pgc
      message("Total peak and TF pairs: ", nrow(celltype_TF_pgc))
      write.table(celltype_TF_pgc, file=file.path(projdir, "/results/celltype", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.", celltype_name, ".TFs.tsv")), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
      celltype_pgc_TFs <- unique(celltype_TF_pgc$gene)
      message("Total pgc TFs: ", length(celltype_pgc_TFs))
      if (length(celltype_pgc_TFs)==0) {
        TF_pgc_TFs_list[[celltype_name]] <- data.frame(TF=NULL, celltype=NULL, gene=NULL)
      } else {
        TF_pgc_TFs_list[[celltype_name]] <- data.frame(TF=tf_name, celltype=celltype_name, gene=celltype_pgc_TFs)
      }
    } else if (file.exists(celltype_all_pgc_file)&file.info(celltype_all_pgc_file)$size==0) {
      message("No proximal peaks found for ", celltype_name)
      write.table(NULL, file=file.path(projdir, "/results/celltype", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.", celltype_name, ".TFs.tsv")), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
      TF_pgc_TFs_list[[celltype_name]] <- data.frame(TF=NULL, celltype=NULL, gene=NULL)
    } else {
      message("No proximal peaks file found for ", celltype_name)
      write.table(NULL, file=file.path(projdir, "/results/celltype", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.", celltype_name, ".TFs.tsv")), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
      TF_pgc_TFs_list[[celltype_name]] <- data.frame(TF=NULL, celltype=NULL, gene=NULL)
    }
  }

  do.call(rbind, TF_pgc_TFs_list) -> TF_pgc_TFs_all
  write.csv(TF_pgc_TFs_all, file=file.path(projdir, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.celltypes.TFs.csv")), row.names=FALSE, quote=FALSE)
  TF_pgc_TFs_num <- data.frame(celltype=factor(names(TF_pgc_TFs_list), levels=rev(names(TF_pgc_TFs_list))), number=unlist(lapply(TF_pgc_TFs_list,nrow)))
  write.csv(TF_pgc_TFs_num, file=file.path(projdir, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.celltypes.TFs.num.csv")), row.names=FALSE, quote=FALSE)
}

# check expression of genes and TFs 
for (tf_name in TFs) {
  print(tf_name)
  RNA_cpm <- readRDS(file=file.path(projdir, "/data/", paste0(tf_name, ".RNA.cpm.rds")))
  TF_proximal_gene_all <- read.csv(file=file.path(projdir, paste0(tf_name, ".celltypes.proximal.genes.csv")))
  TF_proximal_TFs_all <- read.csv(file=file.path(projdir, paste0(tf_name, ".celltypes.proximal.TFs.csv")))
  all_pgc_all <- read.csv(file=file.path(projdir, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.celltypes.csv")))
  TF_pgc_TFs_all <- read.csv(file=file.path(projdir, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.celltypes.TFs.csv")))
  TF_pgc_gene_all <- all_pgc_all %>% mutate(TF=tf_name) %>% dplyr::select(TF, celltype, gene) %>% distinct()
  write.csv(TF_pgc_gene_all, file=file.path(projdir, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.celltypes.genes.csv")), row.names=FALSE, quote=FALSE)

  RNA_cpm_melt <- as.data.frame(RNA_cpm) %>% rownames_to_column(var="gene") %>% pivot_longer(-gene, names_to="celltype", values_to="cpm")
  TF_proximal_gene_all %>% left_join(RNA_cpm_melt, by=join_by(gene, celltype)) -> TF_proximal_gene_all_exp
  TF_proximal_TFs_all %>% left_join(RNA_cpm_melt, by=join_by(gene, celltype)) -> TF_proximal_TFs_all_exp
  TF_pgc_gene_all %>% left_join(RNA_cpm_melt, by=join_by(gene, celltype)) -> TF_pgc_gene_all_exp
  TF_pgc_TFs_all %>% left_join(RNA_cpm_melt, by=join_by(gene, celltype)) -> TF_pgc_TFs_all_exp
  write.csv(TF_proximal_gene_all_exp, file=file.path(projdir, paste0(tf_name, ".celltypes.proximal.genes.exp.csv")), row.names=FALSE, quote=FALSE)
  write.csv(TF_proximal_TFs_all_exp, file=file.path(projdir, paste0(tf_name, ".celltypes.proximal.TFs.exp.csv")), row.names=FALSE, quote=FALSE)
  write.csv(TF_pgc_gene_all_exp, file=file.path(projdir, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.celltypes.genes.exp.csv")), row.names=FALSE, quote=FALSE)
  write.csv(TF_pgc_TFs_all_exp, file=file.path(projdir, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.celltypes.TFs.exp.csv")), row.names=FALSE, quote=FALSE)
  TF_proximal_gene_all_exp %>% filter(cpm>0) -> TF_proximal_gene_exp_filter
  TF_proximal_TFs_all_exp %>% filter(cpm>0) -> TF_proximal_TFs_exp_filter
  TF_pgc_gene_all_exp %>% filter(cpm>0) -> TF_pgc_gene_exp_filter
  TF_pgc_TFs_all_exp %>% filter(cpm>0) -> TF_pgc_TFs_exp_filter
  write.csv(TF_proximal_gene_exp_filter, file=file.path(projdir, paste0(tf_name, ".celltypes.proximal.genes.exp.filter.csv")), row.names=FALSE, quote=FALSE)
  write.csv(TF_proximal_TFs_exp_filter, file=file.path(projdir, paste0(tf_name, ".celltypes.proximal.TFs.exp.filter.csv")), row.names=FALSE, quote=FALSE)
  write.csv(TF_pgc_gene_exp_filter, file=file.path(projdir, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.celltypes.genes.exp.filter.csv")), row.names=FALSE, quote=FALSE)
  write.csv(TF_pgc_TFs_exp_filter, file=file.path(projdir, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.celltypes.TFs.exp.filter.csv")), row.names=FALSE, quote=FALSE)

  TF_proximal_gene_all_exp %>% filter(cpm>0) %>% group_by(celltype) %>% summarise(number=n()) %>% as.data.frame() -> TF_proximal_gene_expnum
  TF_proximal_TFs_all_exp %>% filter(cpm>0) %>% group_by(celltype) %>% summarise(number=n()) %>% as.data.frame() -> TF_proximal_TFs_expnum
  TF_pgc_gene_all_exp %>% filter(cpm>0) %>% group_by(celltype) %>% summarise(number=n()) %>% as.data.frame() -> TF_pgc_gene_expnum
  TF_pgc_TFs_all_exp %>% filter(cpm>0) %>% group_by(celltype) %>% summarise(number=n()) %>% as.data.frame() -> TF_pgc_TFs_expnum

  TF_proximal_gene_expnum_all <- data.frame(celltype=factor(celltypes, levels=rev(celltypes)))
  TF_proximal_gene_expnum_all <- TF_proximal_gene_expnum_all %>% left_join(TF_proximal_gene_expnum, by=join_by(celltype)) %>% replace_na(list(number=0)) %>% mutate(celltype=factor(celltype, levels=rev(celltypes)))
  TF_proximal_TFs_expnum_all <- data.frame(celltype=factor(celltypes, levels=rev(celltypes)))
  TF_proximal_TFs_expnum_all <- TF_proximal_TFs_expnum_all %>% left_join(TF_proximal_TFs_expnum, by=join_by(celltype)) %>% replace_na(list(number=0)) %>% mutate(celltype=factor(celltype, levels=rev(celltypes)))
  TF_pgc_gene_expnum_all <- data.frame(celltype=factor(celltypes, levels=rev(celltypes)))
  TF_pgc_gene_expnum_all <- TF_pgc_gene_expnum_all %>% left_join(TF_pgc_gene_expnum, by=join_by(celltype)) %>% replace_na(list(number=0)) %>% mutate(celltype=factor(celltype, levels=rev(celltypes)))
  TF_pgc_TFs_expnum_all <- data.frame(celltype=factor(celltypes, levels=rev(celltypes)))
  TF_pgc_TFs_expnum_all <- TF_pgc_TFs_expnum_all %>% left_join(TF_pgc_TFs_expnum, by=join_by(celltype)) %>% replace_na(list(number=0)) %>% mutate(celltype=factor(celltype, levels=rev(celltypes)))
  write.csv(TF_proximal_gene_expnum_all, file=file.path(projdir, paste0(tf_name, ".celltypes.proximal.genes.exp.num.csv")), row.names=FALSE, quote=FALSE)
  write.csv(TF_proximal_TFs_expnum_all, file=file.path(projdir, paste0(tf_name, ".celltypes.proximal.TFs.exp.num.csv")), row.names=FALSE, quote=FALSE)
  write.csv(TF_pgc_gene_expnum_all, file=file.path(projdir, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.celltypes.genes.exp.num.csv")), row.names=FALSE, quote=FALSE)
  write.csv(TF_pgc_TFs_expnum_all, file=file.path(projdir, "/results/", paste0(tf_name, ".PGC.cor.",cor.method, ".sig.celltypes.TFs.exp.num.csv")), row.names=FALSE, quote=FALSE)
}

# network
TF_RNA_combined_count <- AggregateExpression(TF_RNA_combined, assays="RNA", slot="count", group.by="celltype")$RNA
TF_RNA_combined_cpm <- edgeR::cpm(TF_RNA_combined_count)

# TF network
TF_proximal_network_list <- list()
for (TF in TFs) {
  message(TF)
  TF_proximal_network_file <- file.path(projdir, "results/", paste0(TF, ".celltypes.proximal.TFs.exp.filter.csv"))
  
  if (file.exists(TF_proximal_network_file)) {
    TF_proximal_network <- read.csv(TF_proximal_network_file)
    TF_proximal_network_list[[TF]] <- TF_proximal_network
  } else {
    message(paste0("No proximal network file for ", TF))
  }
}

TF_proximal_network_all <- do.call(rbind, TF_proximal_network_list)
write.csv(TF_proximal_network_all, file=file.path(projdir, "results/", "TFs.celltypes.proximal.TFs.exp.filter.csv"), quote=FALSE, row.names=FALSE)

TF_pgc_network_list <- list()
for (TF in TFs) {
  message(TF)
  TF_pgc_network_file <- file.path(projdir, "results", paste0(TF, ".PGC.cor.pearson.sig.celltypes.TFs.exp.filter.csv"))

  if (file.exists(TF_pgc_network_file)) {
    TF_pgc_network <- read.csv(TF_pgc_network_file)
    TF_pgc_network_list[[TF]] <- TF_pgc_network
  } else {
    message(paste0("No pgc network file for ", TF))
  }
}

TF_pgc_network_all <- do.call(rbind, TF_pgc_network_list)
write.csv(TF_pgc_network_all, file=file.path(projdir, "results", "TFs.PGC.celltypes.TFs.exp.filter.csv"), quote=FALSE, row.names=FALSE)

TF_network_all <- rbind(TF_proximal_network_all %>% mutate(type="proximal"), TF_pgc_network_all %>% mutate(type="pgc"))
write.csv(TF_network_all, file=file.path(projdir, "results", "TFs.celltypes.TFs.exp.filter.csv"), quote=FALSE, row.names=FALSE)

# all gene network
TF_proximal_network_list <- list()
for (TF in TFs) {
  message(TF)
  TF_proximal_network_file <- file.path(projdir, "results/", paste0(TF, ".celltypes.proximal.genes.exp.filter.csv"))

  if (file.exists(TF_proximal_network_file)) {
    TF_proximal_network <- read.csv(TF_proximal_network_file)
    TF_proximal_network_list[[TF]] <- TF_proximal_network
  } else {
    message(paste0("No proximal network file for ", TF))
  }
}

TF_proximal_network_all <- do.call(rbind, TF_proximal_network_list)
write.csv(TF_proximal_network_all, file=file.path(projdir, "results/", "TFs.celltypes.proximal.genes.exp.filter.csv"), quote=FALSE, row.names=FALSE)

TF_pgc_network_list <- list()
for (TF in TFs) {
  message(TF)
  TF_pgc_network_file <- file.path(projdir, "results", paste0(TF, ".PGC.cor.pearson.sig.celltypes.genes.exp.filter.csv"))

  if (file.exists(TF_pgc_network_file)) {
    TF_pgc_network <- read.csv(TF_pgc_network_file)
    TF_pgc_network_list[[TF]] <- TF_pgc_network
  } else {
    message(paste0("No pgc network file for ", TF))
  }
}

TF_pgc_network_all <- do.call(rbind, TF_pgc_network_list)
write.csv(TF_pgc_network_all, file=file.path(projdir, "results", "TFs.PGC.celltypes.genes.exp.filter.csv"), quote=FALSE, row.names=FALSE)

TF_network_all <- rbind(TF_proximal_network_all %>% mutate(type="proximal"), TF_pgc_network_all %>% mutate(type="pgc"))
write.csv(TF_network_all, file=file.path(projdir, "results", "TFs.celltypes.genes.exp.filter.csv"), quote=FALSE, row.names=FALSE)


# filter gene and TF with logcpm<1 for all cell types
TF_network_allTFs <- read.csv(file=file.path(projdir, "results", "TFs.celltypes.TFs.exp.filter.csv"))
TF_RNA_combined_logcpm <- log1p(TF_RNA_combined_cpm)
markergenes <- read.csv(file=file.path(projdir, "results/markergene_importantgene.csv"))
TF_network_allgenes <- read.csv(file=file.path(projdir, "results", "TFs.celltypes.genes.exp.filter.csv"))

for (celltype_name in celltypes) {
  message(celltype_name)
  if (celltype_name %in% c("CLAGL", "ITL6GL", "ITL5GL", "ITL45GL", "ITL23GL", "ETL5GL", "L6bGL", "CTL6GL", "NPL5GL")) {
    cell_class <- "GLUT"
  } else if (celltype_name %in% c("VIPGA", "PVGA", "SSTGA")) {
    cell_class <- "GABA"
  } else {
    cell_class <- celltype_name
  }
  message(cell_class)

  markergenes %>% filter(type==cell_class) %>% distinct %>% pull(gene) -> cellclass_markergenes

  TF_network_allTFs %>% filter(celltype==celltype_name) %>% dplyr::select(TF, gene) %>% distinct -> TF_network_celltype_TFs
  TF_network_allgenes %>% filter(celltype==celltype_name&gene%in%cellclass_markergenes) %>% dplyr::select(TF, gene) %>% distinct -> TF_network_celltype_markers
  TF_network_celltype <- rbind(TF_network_celltype_TFs, TF_network_celltype_markers) %>% distinct

  TF_RNA_combined_logcpm %>% as.data.frame %>% dplyr::select(celltype_name) %>% tibble::rownames_to_column("gene") %>% mutate_(logCPM=celltype_name) -> TF_RNA_cpm_celltype

  TF_RNA_cpm_celltype %>% filter(logCPM>1) -> TF_RNA_cpm_celltype_filter

  TF_network_celltype %>% mutate(TF=str_to_title(TF)) %>% filter(gene %in% TF_RNA_cpm_celltype_filter$gene) %>% distinct -> TF_network_celltype_filter

  TF_network_celltype_filter %>% group_by(gene) %>% summarise(number=n()) %>% arrange(desc(number)) -> TF_network_celltype_genenum

  TF_network_celltype_filter_select <- TF_network_celltype_filter
  message(paste0("Number of edges: ", nrow(TF_network_celltype_filter_select)))
  TF_network_celltype_genes_select <- union(TF_network_celltype_filter$TF, TF_network_celltype_filter$gene)
  TF_network_celltype_genes_select %>% as.data.frame %>% dplyr::rename(gene=".") %>% mutate(type=case_when(gene %in% TF_names ~ "TF12", gene %in% cellclass_markergenes ~ "markergene", TRUE ~ "otherTF")) -> TF_network_celltype_genes_select_df
  message(paste0("Number of genes: ", nrow(TF_network_celltype_genes_select_df)))


  write.csv(TF_network_celltype_filter_select, file=file.path(projdir, "results/celltype/logcpm1", paste0(celltype_name, ".TFs.network.select.csv")), quote=FALSE, row.names=FALSE)
  write.csv(TF_RNA_cpm_celltype, file=file.path(projdir, "results/celltype/logcpm1", paste0(celltype_name, ".genes.cpm.csv")), quote=FALSE, row.names=FALSE)
  write.csv(TF_network_celltype_genes_select_df, file=file.path(projdir, "results/celltype/logcpm1", paste0(celltype_name, ".network.genes.type.csv")), quote=FALSE, row.names=FALSE)
  write.csv(TF_network_celltype_genenum, file=file.path(projdir, "results/celltype/logcpm1", paste0(celltype_name, ".network.genes.number.csv")), quote=FALSE, row.names=FALSE)
}


