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
# Perform spatial registration of the Xenium Banksy
# clusters against the Visium layers using the spatialLIBD package
# Enrichment analysis reference: https://github.com/LieberInstitute/spatialDLPFC_SCZ/blob/main/code/analysis/04_SpD_marker_genes/PB_analysis_SpD_enrichment.r
#########################################################################

# Create pseudobulked SPE
# read in the data
spe <- readRDS(here("processed-data", "01_build_spe", "raw_spe_N24.RDS"))
outlier_ids <- read.csv(here("processed-data", "02_xenium_qc", "outlier_ids.csv"))$x
spe <- spe[, -which(colnames(spe) %in% outlier_ids)]
rm (outlier_ids)
spe <- spe[which(rowData(spe)$Type=="Gene Expression"), ]

#### Banksy parameters ####
lambda <- 0.8
res <- 0.5
compute_agf <- FALSE
use_agf <- FALSE
cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(use_agf), lambda, res)

clusts <- read.csv(here("processed-data", "03_clustering", sprintf("banksy_clustering_lambda%s_res%s.csv", lambda, res)))
colData(spe)[[cnm]] <- factor(paste0("spd_",as.character(clusts$V1)))

# Add some slide and date related metadata
slide_id <- unlist(lapply(strsplit(spe$Sample, split="/"), "[", 9))
slide_id <- unlist(lapply(strsplit(slide_id, split="__Br"), "[", 1))
spe$slide_id <- slide_id

run_date <- unlist(lapply(strsplit(spe$Sample, split="/"), "[", 7))
spe$run_date <- run_date

spe_pseudo <- registration_pseudobulk(
      spe,
      var_registration = cnm,
      var_sample_id = "BrNum",
      covars = c("Dx", "Age", "Sex", "slide_id", "run_date"),
      min_ncells = 10,
      pseudobulk_rds_file = here(
        "processed-data", "03_clustering",
        paste0("spe_pseudo_donor_domain_", "Banksy", cnm, ".rds")
      )
    )

#spe_pseudo[[cnm]] <- droplevels(spe_pseudo[[cnm]])
spe_pseudo$Age <- as.numeric(spe_pseudo$Age)

# Create the DE model
dx_mod <-
    registration_model(
      spe_pseudo,
      covars = c("Age", "Sex", "slide_id"),
      var_registration = cnm
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
    var_registration = cnm,
    gene_ensembl = "ID",
    gene_name = "Symbol")

rownames(layer_res) <- layer_res$ensembl
# Save the xenium enrichment results
saveRDS(layer_res, here("processed-data", 
                    "04_label_transfer", 
                    "banksy_layer_enrichment_results.rds"))


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

manual_cor_ordered <- manual_cor[c("spd_7","spd_6", "spd_5",  "spd_2", "spd_1", "spd_12", "spd_4", "spd_8", "spd_11","spd_3", "spd_10", "spd_9" ),]
pdf(here("plots", "03_clustering", 
            "banksy_layer_enrichment_correlation_manual_allgenes.pdf"), width = 10, height = 10)
layer_stat_cor_plot(
    manual_cor_ordered,
    max = max(manual_cor)
  )
dev.off()
