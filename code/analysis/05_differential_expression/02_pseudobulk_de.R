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
spe_pseudo <- readRDS(here("processed-data", "05_differential_expression", "spe_pseudo_donor_domain_spaTransfer_k50_smoothed_predictions.rds"))
spe_pseudo$Dx <- factor(spe_pseudo$Dx)
spe_pseudo$predictions_smooth <- factor(spe_pseudo$predictions_smooth)
spe_pseudo$Age <- as.numeric(spe_pseudo$Age)
# Create the DE model
dx_mod <-
    registration_model(
      spe_pseudo,
      covars = c("predictions_smooth", "Age", "Sex", "slide_id"),
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
    covars = c("predictions_smooth", "Age", "Sex", "slide_id"),
    var_registration = "Dx",
    gene_ensembl = "ID",
    gene_name = "Symbol"
  )

# Save results
write.csv(dx_res, here("processed-data", "05_differential_expression", "donor_domain_level_pseudobulk_Dx_DEGs_spaTransfer_k50_smoothed_predictions.csv"))


pdf(here("plots", "05_differential_expression", "donor_domain_level_pseudobulk_Dx_DEGs_k50_volcano.pdf"))
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
        filter(gene %in% c("FKBP5","SERPINA3","GSTM5", "HIF3A", "CHI3L1",
          "MBP", "XRRA1", "IFITM3", "NPTX2", "TNFSF10", "KLF2", "SLC44A1")) %>%
        mutate(
          Direction = case_when(
            fdr_SCZ >= 0.1 ~ "N.S. (FDR > 0.1)",
            logFC_SCZ > 0 & fdr_SCZ < 0.1 ~ "Upregulated",
            logFC_SCZ < 0 & fdr_SCZ < 0.1 ~ "Downregulated"
        )) %>%
        select(ensembl, gene, ends_with("SCZ"), Direction)

  n_sig_gene <- dx_res |>
    filter(fdr_SCZ <= 0.1) |>
    nrow()
 dx_res <- dx_res %>%
  mutate(
    Direction = case_when(
      fdr_SCZ >= 0.1 ~ "N.S. (FDR > 0.1)",
      logFC_SCZ > 0 & fdr_SCZ < 0.1 ~ "Upregulated",
      logFC_SCZ < 0 & fdr_SCZ < 0.1 ~ "Downregulated"
    ))
  print(
    ggplot(
      dx_res,
      aes(
        x = logFC_SCZ, y = -log10(fdr_SCZ),
        color = Direction
      )
    ) +
      geom_point(alpha = 0.8, show.legend=TRUE) +
      geom_label_repel(size=6,
        data = sig_gene_df,
        aes(label = gene),
        force = 2,
        nudge_y = 0.1, show.legend=FALSE
      ) +
      scale_color_manual(values=c("Upregulated"="red",
               "Downregulated"="blue", "N.S. (FDR > 0.1)"="grey"),
               guide = guide_legend(
              override.aes = list(size = 8)  # legend points size
         )) +
     theme_minimal()+
     theme(axis.text.x=element_text(size=18),
            axis.text.y=element_text(size=18),
            axis.title.x=element_text(size=20),
            axis.title.y=element_text(size=20),
            legend.title=element_text(size=18),
            legend.position="bottom", legend.direction="horizontal")
    
  )
dev.off()


## check whether the sig genes are in SHK's gene list
panel_markers <- readxl::read_xlsx((here("raw-data", 
        "experiment_info", 
        "Xenium_SHK_celltype_REannot_2025-04-13.xlsx")), sheet=2)

panel_markers <- panel_markers %>%
    as.data.frame() %>%
    mutate(cell_type_updated=case_when(cell_type_updated=="NA" ~ NA,
                                        TRUE ~ cell_type_updated))


panel_markers <- panel_markers %>%
  filter(Gene %in% sig_gene_df$gene) %>%
  select(Gene, dx_deg,  cell_type_updated, layer_marker)%>%
  column_to_rownames("Gene")

sig_gene_df <- sig_gene_df %>%
  mutate(xenium_direction = case_when(logFC_SCZ > 0 ~ "up",
                               logFC_SCZ < 0 ~ "down"))%>%
  select(xenium_direction)

sig_df <- merge(panel_markers, sig_gene_df, by="row.names")

sig_df <- sig_df[c("Row.names", "dx_deg", "xenium_direction", "cell_type_updated", "layer_marker")]
print(sig_df)