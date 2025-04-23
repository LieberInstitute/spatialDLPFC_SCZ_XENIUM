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
# Perform DE analysis on the pseudobulked object
# Reference: https://github.com/LieberInstitute/spatialDLPFC_SCZ/blob/main/code/analysis/pseudobulk_dx/PB_analysis_dx_gene.R
####################################################################################################

# Load the pseudobulked SPE
spe_pseudo <- readRDS(here("processed-data", "07_cell_type_de", "spe_pseudo_donor_celltype_clust_M0_lam0.1_k50_res0.7.rds"))
spe_pseudo$Dx <- factor(spe_pseudo$Dx)
spe_pseudo$predictions_smooth <- factor(spe_pseudo$annots)
spe_pseudo$Age <- as.numeric(spe_pseudo$Age)


#### Banksy parameters ####
lambda <- 0.1
res <- 0.7 # higher = more clusters
use_agf <- FALSE
cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(use_agf), lambda, res)


# Create the DE model
dx_mod <-
    registration_model(
      spe_pseudo,
      covars = c("annots", "Age", "Sex", "slide_id"),
      var_registration = "Dx"
    )


  dx_block_cor <-
    registration_block_cor(
      spe_pseudo,
      registration_model = dx_mod,
      var_sample_id = "BrNum"
    )

  dx_res <- registration_stats_enrichment(
    spe_pseudo,
    block_cor = dx_block_cor,
    covars = c("annots", "Age", "Sex", "slide_id"),
    var_registration = "Dx",
    gene_ensembl = "ID",
    gene_name = "Symbol"
  )

# Save results
write.csv(dx_res, here("processed-data", "07_cell_type_de", sprintf("donor_cell_type_level_pseudobulk_Dx_DEGs_.csv", cnm)))


pdf(here("plots", "07_cell_type_de", sprintf("donor_cell_type_level_pseudobulk_Dx_DEGs_%s.pdf", cnm)))

## Volcano Plot ----
  out_gene_df <- dx_res |>
    arrange(fdr_SCZ) |>
    slice_head(n = 1)

  impl_gene_df <- dx_res |>
    filter(gene %in% c(
      "PVALB",
      "NOS1",
      "SST",
      "CHODL",
      "GRIN2A",
      "SV2A",
      "DLG4",
      "C4A",
      "C3"
    )) |>
    select(ensembl, gene, ends_with("SCZ"))
    
    sig_gene_df <- dx_res |>
        filter(fdr_SCZ <= 0.1) |>
        select(ensembl, gene, ends_with("SCZ"))

  n_sig_gene <- dx_res |>
    filter(fdr_SCZ <= 0.1) |>
    nrow()

  print(
    ggplot(
      dx_res,
      aes(
        x = logFC_SCZ, y = -log10(fdr_SCZ),
        color = fdr_SCZ <= 0.05
      )
    ) +
      geom_point(alpha = 0.8) +
      geom_label_repel(
        data = sig_gene_df,
        aes(label = gene),
        force = 2,
        nudge_y = 0.1
      ) +
      labs(
        title = paste0(
          "Donor-cell type level pseudobulk analysis by Dx - ", 
          " ( ", n_sig_gene, " sig genes)"
        )
      ) +
      theme_minimal()
  )
dev.off()


