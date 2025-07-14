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

spe$domain_annotations <- gsub(spe$domain_annotations, pattern = "\\/", replacement = "_")
layers <- unique(spe$domain_annotations)

spe$cell_types <- as.factor(colData(spe)[["Banksy-clust_M0_lam0.1_k50_res0.7-cell-types"]])
spe$Age <- as.numeric(spe$Age)

cell_types <- unique(spe$cell_types)



for (i in 1:length(layers)){
    layer_use <- layers[[i]]
    pdf(here("plots", "07_cell_type_de", "layer_stratified_de",
        sprintf("%s_and_cell_type_stratified_de.pdf", layer_use)))
    for (j in 1:length(cell_types)){
        cell_type_use <- cell_types[[j]]
        print(sprintf("-------%s-------", cell_type_use))
        spe_sub <- spe[, spe$cell_types == cell_type_use & spe$domain_annotations == layer_use]
        
        #--------------------------
        # Create pseudobulked data
        #--------------------------
        spe_pseudo <- registration_pseudobulk(
                spe_sub,
                var_registration = "cell_types",
                var_sample_id = "BrNum",
                covars = c("Dx", "Age", "Sex", "slide_id", "run_date"),
                min_ncells = 20,
                pseudobulk_rds_file = NULL
            )
        spe_pseudo$cell_types <- as.factor(spe_pseudo$cell_types)
        spe_pseudo$cell_types <- droplevels(spe_pseudo$cell_types)

        # Create the DE model
        dx_mod <- registration_model(
            spe_pseudo,
            covars = c("Age", "Sex", "slide_id"),
            var_registration = "Dx"
        )
        dx_res <- registration_stats_enrichment(
            spe_pseudo,
            block_cor = NA,
            covars = c("Age", "Sex", "slide_id"),
            var_registration = "Dx",
            gene_ensembl = "ID",
            gene_name = "Symbol"
        )
        sig_gene_df <- dx_res |>
            arrange(fdr_SCZ) |>
            slice_head(n=10) 

        print(ggplot(
            dx_res,
            
            aes(
                x = logFC_SCZ, y = -log10(fdr_SCZ),
                color = fdr_SCZ <= 0.05
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
                title = paste0(cell_type_use ," DE within ",
                layer_use)
            ) +theme_minimal()+
            theme(axis.text.x=element_text(size=16),
                axis.text.y=element_text(size=16),
                axis.title.x=element_text(size=16),
                axis.title.y=element_text(size=16),
                plot.title = element_text(hjust = 0.5, size=20))
            )

    }
    dev.off()

}
