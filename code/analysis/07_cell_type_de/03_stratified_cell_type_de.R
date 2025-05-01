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
# Perform DE analysis stratified by cell type on the pseudobulked object
# Reference: https://github.com/LieberInstitute/spatialDLPFC_SCZ/blob/main/code/analysis/pseudobulk_dx/PB_analysis_dx_gene.R
####################################################################################################


# Read in the pseudobulked SPE
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


# Run PCA and plot PC plot to see if there is variation driven by Dx within the cell types
set.seed(1108)
spe_pseudo <- scater::runPCA(spe_pseudo, ncomponents=50)
pdf(here("plots", "07_cell_type_de", 
        sprintf("stratified_donor_cell_type_level_PCA_%s.pdf", cnm)),
        height=15, width=20)
scater::plotPCA(spe_pseudo, colour_by = "Dx", ncomponents=5)
scater::plotPCA(spe_pseudo, colour_by = "annots", shape_by="Dx", ncomponents=5)
scater::plotPCA(spe_pseudo, colour_by = "BrNum", ncomponents=5)
dev.off()

pdf(here("plots", "07_cell_type_de", 
        sprintf("stratified_donor_cell_type_level_pseudobulk_Dx_DEGs_%s.pdf", cnm)))
# Loop over the cell types and perform DE within each cell type

cell_types <- unique(spe_pseudo$annots)
dx_results <- list()
for(i in 1:length(cell_types)){
    cell_type_use <- cell_types[[i]]

    spe_pseudo_sub <- spe_pseudo[,spe_pseudo$annots==cell_type_use]

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
    cell_type_use_fname <- gsub("/", "-", cell_type_use)
    cell_type_use_fname <- gsub("\\:", "-", cell_type_use_fname)
    write.csv(dx_res, here("processed-data", "07_cell_type_de", 
            "stratified_results", 
            sprintf("stratified_donor_cell_type_level_pseudobulk_Dx_DEGs_%s.csv", cell_type_use_fname)))


   
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
            title = paste0(cell_type_use, 
            " Donor-cell type level pseudobulk analysis by Dx - ", 
            " ( ", n_sig_gene, " sig genes)"
            )
        ) +
        theme_minimal()
  )
    dx_res$cell_type <- cell_type_use
    dx_res_use <- dx_res %>%
        select(gene, logFC_SCZ, cell_type)
    dx_results[[i]] <- dx_res_use

}
dev.off()

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
    select(Gene, dx_deg)

dx_df <- do.call(rbind, dx_results)
dx_mat <- dx_df %>%
    as.data.frame() %>% 
    pivot_wider(names_from=cell_type, values_from=logFC_SCZ)%>%
    filter(gene %in% dx_degs$Gene) %>%
    column_to_rownames(var="gene") %>%
    as.matrix()

# Create row annotation of directionality 
# Ensure gene names are rownames
dx_deg_vec <- dx_degs %>%
  filter(Gene %in% rownames(dx_mat)) %>%
  distinct(Gene, .keep_all = TRUE) %>% # avoid duplicated rows
  column_to_rownames("Gene")

# Reorder to match dx_mat
dx_deg_vec <- dx_deg_vec[rownames(dx_mat),]
dx_col_fun <- c("Dx_DEG_Up" = "red", "Dx_DEG_Down" = "blue")

row_anno <- rowAnnotation(
  DEG = dx_deg_vec,
  col = list(DEG = dx_col_fun),
  show_annotation_name = TRUE
)

pdf(here("plots", "07_cell_type_de", "stratified_cell_type_de_heatmap.pdf"),
    height=25, width=10)
ComplexHeatmap::Heatmap(dx_mat, right_annotation = row_anno)
dev.off()
