suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("stringr"))
library(Seurat)
library(Matrix)
library(monocle3)
library(cicero)
library(dyno)
library(tidyverse)
library(slingshot)

projdir <- "/tscc/projects/PairedTF/"

# RNA
RNA_input <- file.path(projdir, "PairedTF.final.reidents.rds")
RNA_PairedTF <- readRDS(RNA_input)

RNA_PairedTF_OPCOGC <- subset(RNA_PairedTF, idents=c("OPC", "OGC"))
RNA_PairedTF_OPCOGC$sampleID <- str_split(colnames(RNA_PairedTF_OPCOGC), "[_]", simplify=TRUE)[,3]
RNA_PairedTF_OPCOGC$barcode <- str_split(colnames(RNA_PairedTF_OPCOGC), "_", simplify=TRUE)[,4]
RNA_PairedTF_OPCOGC$cell <- paste(RNA_PairedTF_OPCOGC$sampleID, RNA_PairedTF_OPCOGC$barcode, sep=".")

# dynverse
RNA_PairedTF_OPCOGC_counts <- RNA_PairedTF_OPCOGC@assays$RNA@counts
RNA_PairedTF_OPCOGC_expression <- RNA_PairedTF_OPCOGC@assays$RNA@data
RNA_PairedTF_OPCOGC_counts_filter <- t(as.matrix(RNA_PairedTF_OPCOGC_counts)[Matrix::rowSums(RNA_PairedTF_OPCOGC_counts) != 0, Matrix::colSums(RNA_PairedTF_OPCOGC_counts) != 0]) 
RNA_PairedTF_OPCOGC_expression_filter <- t(as.matrix(RNA_PairedTF_OPCOGC_expression)[Matrix::rowSums(RNA_PairedTF_OPCOGC_expression) != 0, Matrix::colSums(RNA_PairedTF_OPCOGC_expression) != 0]) 
RNA_PairedTF_OPCOGC_meta_filter <- RNA_PairedTF_OPCOGC_meta[rownames(RNA_PairedTF_OPCOGC_expression_filter), ]

identical(rownames(RNA_PairedTF_OPCOGC_counts_filter), rownames(RNA_PairedTF_OPCOGC_expression_filter)) && identical(colnames(RNA_PairedTF_OPCOGC_counts_filter), colnames(RNA_PairedTF_OPCOGC_expression_filter))

RNA_PairedTF_OPCOGC_dataset <- wrap_expression(
  counts = RNA_PairedTF_OPCOGC_expression_filter,
  expression = RNA_PairedTF_OPCOGC_expression_filter
)

guidelines <- guidelines_shiny(RNA_PairedTF_OPCOGC_dataset)
methods_selected <- guidelines$methods_selected

# MST method
RNA_PairedTF_OPCOGC_model <- infer_trajectory(RNA_PairedTF_OPCOGC_dataset, "MST")

RNA_PairedTF_OPCOGC_model <- RNA_PairedTF_OPCOGC_model %>% add_dimred(dyndimred::dimred_mds, expression_source = RNA_PairedTF_OPCOGC_dataset$expression)
plot_dimred(
  RNA_PairedTF_OPCOGC_model, 
  expression_source = RNA_PairedTF_OPCOGC_dataset$expression, 
  grouping = RNA_PairedTF_OPCOGC_meta_filter$celltype
)

ciliated_genes <- c("Mobp", "Plp1", "Pdgfra", "Vcan")

plot_dimred(
  RNA_PairedTF_OPCOGC_model, 
  expression_source = RNA_PairedTF_OPCOGC_dataset$expression,
  feature_oi = ciliated_genes[1]
)

plot_dimred(
  RNA_PairedTF_OPCOGC_model, 
  expression_source = RNA_PairedTF_OPCOGC_dataset$expression, 
  color_cells = "feature",
  feature_oi = "Plp1",
  color_density = "grouping",
  grouping = RNA_PairedTF_OPCOGC_meta_filter$celltype,
  label_milestones = FALSE
)

RNA_PairedTF_OPCOGC_model <- RNA_PairedTF_OPCOGC_model %>% add_root_using_expression(c("Pdgfra"), RNA_PairedTF_OPCOGC_dataset$expression)
RNA_PairedTF_OPCOGC_model <- label_milestones_markers(
  RNA_PairedTF_OPCOGC_model,
  markers = list(
    OPC = c("Pdgfra"),
    OGC = c("Plp1")),
  RNA_PairedTF_OPCOGC_dataset$expression
)

plot_heatmap(
  RNA_PairedTF_OPCOGC_model,
  expression_source = RNA_PairedTF_OPCOGC_dataset$expression,
  grouping = RNA_PairedTF_OPCOGC_meta_filter$celltype,
  features_oi = 50
)

branch_feature_importance <- calculate_branch_feature_importance(RNA_PairedTF_OPCOGC_model, expression_source=RNA_PairedTF_OPCOGC_dataset$expression)

OGC_features <- branch_feature_importance %>% 
  filter(to == names(RNA_PairedTF_OPCOGC_model$milestone_labelling)[which(RNA_PairedTF_OPCOGC_model$milestone_labelling =="OGC")]) %>% 
  top_n(50, importance) %>% 
  pull(feature_id)

plot_heatmap(
  RNA_PairedTF_OPCOGC_model, 
  expression_source = RNA_PairedTF_OPCOGC_dataset$expression, 
  features_oi = OGC_features
)

branching_milestone <- RNA_PairedTF_OPCOGC_model$milestone_network %>% group_by(from) %>% filter(n() > 1) %>% pull(from) %>% first()

branch_feature_importance <- calculate_branching_point_feature_importance(RNA_PairedTF_OPCOGC_model, expression_source=RNA_PairedTF_OPCOGC_dataset$expression, milestones_oi = branching_milestone)

branching_point_features <- branch_feature_importance %>% top_n(20, importance) %>% pull(feature_id)

plot_heatmap(
  RNA_PairedTF_OPCOGC_model,
  expression_source = RNA_PairedTF_OPCOGC_dataset$expression,
  features_oi = branching_point_features
)

space <- dyndimred::dimred_mds(RNA_PairedTF_OPCOGC_dataset$expression)

map(branching_point_features[1:12], function(feature_oi) {
  plot_dimred(RNA_PairedTF_OPCOGC_model, dimred = space, expression_source = RNA_PairedTF_OPCOGC_dataset$expression, feature_oi = feature_oi, label_milestones = FALSE) +
    theme(legend.position = "none") +
    ggtitle(feature_oi)
}) %>% patchwork::wrap_plots()

