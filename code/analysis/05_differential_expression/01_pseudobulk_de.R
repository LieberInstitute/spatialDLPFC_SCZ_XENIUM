suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
})

####################################################################################################
# Pseudobulk the SPE to the donor-domain level using the smoothed labels from label transfer.
# Then, follow Boyi's workflow to perform DE analysis.
# Reference: https://github.com/LieberInstitute/spatialDLPFC_SCZ/blob/main/code/analysis/pseudobulk_dx/PB_analysis_dx_gene.R
####################################################################################################

spe <- readRDS(here("processed-data", "01_build_spe", "raw_spe_N24.RDS"))
spds <- read.csv(here("processed-data", "04_label_transfer", "label_transfer_N24_k50_smoothed_labels.csv"))

rownames(spds) <- spds$X 
spds$X <- NULL
colnames(spds) <- c("predictions_smooth")
spe <- spe[,colnames(spe) %in% rownames(spds)]
colData(spe) <- merge(colData(spe), spds, by="row.names", all.x=TRUE)

# Add some slide and date related metadata
slide_id <- unlist(lapply(strsplit(spe$Sample, split="/"), "[", 9))
slide_id <- unlist(lapply(strsplit(slide_id, split="__Br"), "[", 1))
spe$slide_id <- slide_id

run_date <- unlist(lapply(strsplit(spe$Sample, split="/"), "[", 7))
spe$run_date <- run_date


## First, create pseudobulk data:
# Reference: https://github.com/LieberInstitute/spatialDLPFC_SCZ/blob/main/code/analysis/pseudobulk_dx/create_pb_data.R

spe_pseudo <- registration_pseudobulk(
      spe,
      var_registration = "predictions_smooth",
      var_sample_id = "BrNum",
      covars = c("Dx", "Age", "Sex", "slide_id", "slide_id", "run_date"),
      min_ncells = 10,
      pseudobulk_rds_file = here(
        "processed-data", "05_differential_expression",
        paste0("spe_pseudo_", "spaTransfer_k50_smoothed_predictions", ".rds")
      )
    )
