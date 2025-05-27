suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
  library(ComplexHeatmap)
  library(spatialLIBD)
})

# Need to first create a new pseudobulked object for each layer

spe <- readRDS(here("processed-data", "07_cell_type_de", "cleaned_spe_N24_with_cell_type_and_spds.RDS"))
layers <- unique(spe$domain_annotations)
spe$cell_types <- as.factor(colData(spe)[["Banksy-clust_M0_lam0.1_k50_res0.7-cell-types"]])
spe$Age <- as.numeric(spe$Age)


# Read in Sangho's annotations
panel_markers <- readxl::read_xlsx((here("raw-data", 
        "experiment_info", 
        "Xenium_SHK_celltype_REannot_2025-04-13.xlsx")), sheet=2)
panel_markers$`...1` <- NULL

panel_markers <- panel_markers %>%
    as.data.frame() %>%
    mutate(cell_type_updated=case_when(cell_type_updated=="NA" ~ NA,
                                        TRUE ~ cell_type_updated))%>%
    mutate(dx_deg=case_when(dx_deg=="NA" ~ NA,
                                        TRUE ~ dx_deg))

dx_degs <- panel_markers[!is.na(panel_markers$dx_deg),]%>%
    select(Gene, dx_deg)%>%
    mutate(gene=Gene) %>%
    select(gene, dx_deg)


pdf(here("plots", "07_cell_type_de", "layer_stratified_cell_type_as_covariate_de.pdf"))
for (i in 1:length(layers)){
    layer_use <- layers[[i]]
    spe_sub <- spe[, spe$domain_annotations == layer_use]

    #--------------------------
    # Create pseudobulked data
    #--------------------------

    spe_pseudo <- registration_pseudobulk(
      spe_sub,
      var_registration = "cell_types",
      var_sample_id = "BrNum",
      covars = c("Dx", "Age", "Sex", "slide_id", "run_date"),
      min_ncells = 10,
      pseudobulk_rds_file = NULL
    )
    spe_pseudo$cell_types <- as.factor(spe_pseudo$cell_types)
    spe_pseudo$cell_types <- droplevels(spe_pseudo$cell_types)

    # Create the DE model
    dx_mod <- registration_model(
      spe_pseudo,
      covars = c("cell_types", "Age", "Sex", "slide_id"),
      var_registration = "Dx"
    )

    dx_block_cor <- registration_block_cor(
        spe_pseudo,
        registration_model = dx_mod,
        var_sample_id = "BrNum"
    )

  dx_res <- registration_stats_enrichment(
    spe_pseudo,
    block_cor = dx_block_cor,
    covars = c("cell_types", "Age", "Sex", "slide_id"),
    var_registration = "Dx",
    gene_ensembl = "ID",
    gene_name = "Symbol"
  )

  # merge with sangho's markers
  dx_res <- merge(dx_res, dx_degs, by="gene", all.x=TRUE) %>%
    mutate(dx_deg=case_when(
        dx_deg=="Dx_DEG_Up" ~ "up",
        dx_deg=="Dx_DEG_Down" ~ "down",
        TRUE ~ "none")) 
        
    sig_gene_df <- dx_res |>
        arrange(fdr_SCZ) |>
        slice_head(n=10)


    n_sig_gene <- dx_res |>
        filter(fdr_SCZ <= 0.1) |>
        nrow()

    print(
        ggplot(
        dx_res,
        aes(
            x = logFC_SCZ, y = -log10(fdr_SCZ),
            color = fdr_SCZ <= 0.05,
            shape=dx_deg
        )
        ) +
        geom_point(alpha = 0.8) +
        geom_label_repel(size=6,
            data = sig_gene_df,
            aes(label = gene),
            force = 2,
            nudge_y = 0.1
        ) +
        labs(
            title = paste0("Cell-type adjusted DE within ",
             layer_use)
        ) +theme_minimal()+
         theme(axis.text.x=element_text(size=16),
            axis.text.y=element_text(size=16),
            axis.title.x=element_text(size=16),
            axis.title.y=element_text(size=16),
            plot.title = element_text(hjust = 0.5, size=20))
    )

    layer_use <- gsub("\\/", "_", layer_use)
    write.csv(dx_res, here("processed-data", "07_cell_type_de", "layer_stratified_results",
        paste0("layer_stratified_cell_type_as_covariate_de_results_", layer_use, ".csv")), row.names=FALSE)
}
dev.off()



# dx_df <- do.call(rbind, dx_results)
# dx_mat <- dx_df %>%
#     as.data.frame() %>% 
#     pivot_wider(names_from=cell_type, values_from=logFC_SCZ)%>%
#     filter(gene %in% dx_degs$Gene) %>%
#     column_to_rownames(var="gene") %>%
#     as.matrix()