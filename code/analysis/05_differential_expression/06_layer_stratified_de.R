suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
  library(ComplexHeatmap)
})
####################################################################################################
# Perform DE analysis stratified by layer on the pseudobulked object
# Reference: https://github.com/LieberInstitute/spatialDLPFC_SCZ/blob/main/code/analysis/pseudobulk_dx/PB_analysis_dx_gene.R
####################################################################################################
# Load the pseudobulked SPE
spe_pseudo <- readRDS(here("processed-data", "05_differential_expression", "spe_pseudo_donor_domain_spaTransfer_k50_smoothed_predictions.rds"))
spe_pseudo$Dx <- factor(spe_pseudo$Dx)
spe_pseudo$predictions_smooth <- factor(spe_pseudo$predictions_smooth)
spe_pseudo$Age <- as.numeric(spe_pseudo$Age)


cnm <- "spaTransfer_k50_smoothed"
# Run PCA and plot PC plot to see if there is variation driven by Dx within the layers
set.seed(1108)
spe_pseudo <- scater::runPCA(spe_pseudo, ncomponents=50)
pdf(here("plots", "05_differential_expression", 
        sprintf("stratified_donor_layer_level_PCA_%s.pdf", cnm)),
        height=15, width=20)
scater::plotPCA(spe_pseudo, colour_by = "Dx", ncomponents=5)
scater::plotPCA(spe_pseudo, colour_by = "predictions_smooth", shape_by="Dx", ncomponents=5)
scater::plotPCA(spe_pseudo, colour_by = "BrNum", ncomponents=5)
dev.off()

pdf(here("plots", "05_differential_expression", 
        sprintf("stratified_donor_layer_level_pseudobulk_Dx_DEGs_%s.pdf", cnm)))
# Loop over the cell types and perform DE within each cell type

layers <- unique(spe_pseudo$predictions_smooth)
dx_results <- list()
for(i in 1:length(layers)){
    layer_use <- layers[[i]]

    spe_pseudo_sub <- spe_pseudo[,spe_pseudo$predictions_smooth==layer_use]

    # Create the DE model
    dx_mod <-
        registration_model(
        spe_pseudo_sub,
        covars = c("Age", "Sex", "slide_id"),
        var_registration = "Dx"
        )

    dx_res <- registration_stats_enrichment(
        spe_pseudo_sub,
        block_cor = NaN,
        covars = c("Age", "Sex", "slide_id"),
        var_registration = "Dx",
        var_sample_id="BrNum",
        gene_ensembl = "ID",
        gene_name = "Symbol"
    )
 


    # save the results
    layer_use_fname <- gsub("/", "-", layer_use)
    layer_use_fname <- gsub("\\:", "-", layer_use_fname)
    write.csv(dx_res, here("processed-data", "05_differential_expression", 
            "stratified_results", 
            sprintf("stratified_donor_layer_level_pseudobulk_Dx_DEGs_%s.csv", layer_use_fname)))


   
    ## Volcano Plot ----
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
        arrange(fdr_SCZ) |>
        slice_head(n=10)
        #filter(fdr_SCZ <= 0.1) |>
        #select(ensembl, gene, ends_with("SCZ"))

    n_sig_gene <- dx_res |>
        filter(fdr_SCZ <= 0.1) |>
        nrow()

    print(
        ggplot(
        dx_res,
        aes(
            x = logFC_SCZ, y = -log10(fdr_SCZ),
            color = fdr_SCZ <= 0.1
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
            title = paste0(layer_use, 
            " Donor-layer level pseudobulk analysis by Dx - ", 
            " ( ", n_sig_gene, " sig genes)"
            )
        ) +
        theme_minimal()
  )
    dx_res$layer<- layer_use
    dx_res_use <- dx_res %>%
        select(gene, logFC_SCZ, layer)
    dx_results[[i]] <- dx_res_use

}
dev.off()