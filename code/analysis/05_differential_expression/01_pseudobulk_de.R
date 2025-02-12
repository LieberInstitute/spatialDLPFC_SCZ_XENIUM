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
spds <- readRDS(here("processed-data", "04_label_transfer", "label_transfer_N24_k50_smoothed_labels.csv"))

spe <- spe[,colnames(spe) %in% rownames(spds)]
colData(spe) <- merge(colData(spe), spds, by="row.names", all.x=TRUE)

## First, create pseudobulk data:
# Reference: https://github.com/LieberInstitute/spatialDLPFC_SCZ/blob/main/code/analysis/pseudobulk_dx/create_pb_data.R

spe_pseudo <-
    registration_pseudobulk(
      spe,
      var_registration = "predictions_smooth",
      var_sample_id = "sample_id",
      covars = c("Dx", "age", "sex", "lot_num", "slide_id"),
      min_ncells = 10,
      pseudobulk_rds_file = here(
        "processed-data", "rds", "layer_spd",
        paste0("test_spe_pseudo_", .var, ".rds")
      )
    )
}