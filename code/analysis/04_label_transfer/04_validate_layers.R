suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
})
########################################################################
# Perform spatial registration of the Xenium layers against the Visium
# layers using the spatialLIBD package
# Enrichment analysis reference: https://github.com/LieberInstitute/spatialDLPFC_SCZ/blob/main/code/analysis/04_SpD_marker_genes/PB_analysis_SpD_enrichment.r
#########################################################################

# Load the pseudobulked SPE
spe_pseudo <- readRDS(here("processed-data", "05_differential_expression", "spe_pseudo_donor_domain_spaTransfer_k50_smoothed_predictions.rds"))
spe_pseudo$Dx <- factor(spe_pseudo$Dx)
spe_pseudo$predictions_smooth <- factor(spe_pseudo$domain_annotations)
spe_pseudo$Age <- as.numeric(spe_pseudo$Age)


# ------------------------------------------------------------
# First, need to perform enrichment analysis with the spatial domain as the
# registration variable to get domain-speciifc t-statistics
# ------------------------------------------------------------

# Create the DE model
dx_mod <-
    registration_model(
      spe_pseudo,
      covars = c("Age", "Sex", "slide_id"),
      var_registration = "predictions_smooth"
    )

 layer_block_cor <-
    registration_block_cor(
      spe_pseudo,
      registration_model = dx_mod,
      var_sample_id = "BrNum"
    )
# Get the t-statistics for the spatial domain
rownames(spe_pseudo) <- spe_pseudo$ID
layer_res <- registration_stats_enrichment(
    spe_pseudo,
    block_cor = layer_block_cor,
    covars = c("Age", "Sex", "slide_id"),
    var_registration = "predictions_smooth",
    gene_ensembl = "ID",
    gene_name = "Symbol"
  )

rownames(layer_res) <- layer_res$ensembl

# Save the xenium enrichment results
saveRDS(layer_res, here("processed-data", 
                    "04_label_transfer", 
                    "spaTransfer_k50_smoothed_layer_enrichment_results.rds"))


# register to manual annotation
t_stats <- layer_res[, grep("^t_stat_", colnames(layer_res))]
colnames(t_stats) <- gsub("^t_stat_", "", colnames(t_stats))

# Get the manual modeling results
manual_modeling_results <- fetch_data(type = "modeling_results")

manual_cor <- layer_stat_cor(
    t_stats,
    manual_modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    top_n = NULL
  )
# reorder to get things to be diagonal
manual_cor <- manual_cor[c("L1", "L2/3", "L3/4", "L5", "L6", "WMtz", "WM"),]
cor_annot <- annotate_registered_clusters(manual_cor)

pdf(here("plots", "04_label_transfer", 
            "layer_enrichment_correlation_manual_allgenes.pdf"), width = 10, height = 10)
layer_stat_cor_plot(
    manual_cor,
    max = max(manual_cor)
  )
dev.off()


# ------------------------------------------------------------
# Try registration of spatialDLPFC layers to manaul annotations, using only
# the genes from the Xenium panel
# ------------------------------------------------------------

# Get the spatialDLPFC modeling results
spatialDLPFC_modeling_results <- fetch_data(
    type = "spatialDLPFC_Visium_modeling_results"
  )

spatialDLPFC_enrich <- spatialDLPFC_modeling_results$enrichment[, grep("^t_stat_", colnames(spatialDLPFC_modeling_results$enrichment))]
spatialDLPFC_enrich <- spatialDLPFC_enrich[rownames(spatialDLPFC_enrich) %in% rownames(layer_res),]

spatialDLPFC_cor <- layer_stat_cor(
    spatialDLPFC_enrich,
    manual_modeling_results,
    model_type = "enrichment",
    reverse = FALSE,
    top_n = NULL
  )

# reorder to be diagonal
spatialDLPFC_cor <- spatialDLPFC_cor[c("t_stat_Sp09D01", "t_stat_Sp09D02", "t_stat_Sp09D03", "t_stat_Sp09D05", "t_stat_Sp09D08", "t_stat_Sp09D04", "t_stat_Sp09D07", "t_stat_Sp09D06", "t_stat_Sp09D09"),]

cor_annot <- annotate_registered_clusters(spatialDLPFC_cor)
pdf(here("plots", "04_label_transfer", 
            "layer_enrichment_correlation_spatialDLPFC_allgenes.pdf"), width = 10, height = 10)
layer_stat_cor_plot(
    spatialDLPFC_cor,
    max = max(spatialDLPFC_cor)
  )
dev.off()


# ------------------------------------------------------------
# Register against the Visium layers from Boyi's dataset
# ------------------------------------------------------------

vis_enrich <- readRDS("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/rds/04_SpD_marker_genes/test_enrich_PRECAST_07.rds")
vis_use <- list()
vis_use$enrichment <- vis_enrich

vis_cor <- layer_stat_cor(
    t_stats,
    vis_use,
    model_type = "enrichment",
    reverse = FALSE,
    top_n = NULL
  )

# reorder to be diagonal
vis_cor <- vis_cor[colnames(vis_cor),]
pdf(here("plots", "04_label_transfer", 
            "layer_enrichment_correlation_visium_allgenes.pdf"), width = 10, height = 10)
layer_stat_cor_plot(  
    vis_cor,
    max = max(vis_cor)
  )
dev.off()
