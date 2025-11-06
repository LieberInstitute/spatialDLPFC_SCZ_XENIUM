suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
  library(grid)
  library(ComplexHeatmap)
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
manual_cor <- manual_cor[c("L1/M", "L2/3", "L3/4", "L5", "L6", "WMtz", "WM"),]
cor_annot <- annotate_registered_clusters(manual_cor, confidence_threshold = 0.25)
cor_annot$layer_label <- factor(cor_annot$layer_label, levels = c("L1/M", "L2/3", "L3/4", "L5", "L6", "WMtz", "WM"))
# map the layer labels that start with L to "LayerX"
cor_annot$layer_label <- gsub("^L", "Layer", cor_annot$layer_label)
cor_annot$layer_label <- gsub("\\/", "/Layer", cor_annot$layer_label)


create_annotation_matrix <- function(annotation_df, cor_stats_layer) {
    if (!setequal(rownames(cor_stats_layer), annotation_df$cluster)) {
        stop(
            "Query cluster names do not match between rownames(cor_stats_layer) and annotation$cluster.\n
        Be sure the annotation data matches.",
            call. = FALSE
        )
    }

    anno_list <- lapply(
        rownames(cor_stats_layer),
        function(cluster) {
            # look up confidence
            confidence <- annotation_df[match(cluster, annotation_df$cluster), "layer_confidence"]
            print(confidence)
            sym <- ifelse(confidence == "good", "X", "*")
            print(sym)
            # match annotations
            anno <- gsub("\\*", "", annotation_df[match(cluster, annotation_df$cluster), "layer_label"])
            anno_lgl_row <- unlist(lapply(colnames(cor_stats_layer), "%in%", unlist(strsplit(anno, split = "/"))))
            print(anno_lgl_row)
            return(ifelse(anno_lgl_row, sym, ""))
        }
    )
    anno_matrix <- t(data.frame(anno_list))
    rownames(anno_matrix) <- rownames(cor_stats_layer)
    colnames(anno_matrix) <- colnames(cor_stats_layer)

    return(anno_matrix)
}

layer_stat_cor_plot <- function(cor_stats_layer,
    color_max = max(cor_stats_layer),
    color_min = min(cor_stats_layer),
    color_scale = RColorBrewer::brewer.pal(7, "PRGn"),
    query_colors = NULL,
    reference_colors = NULL,
    annotation = NULL,
    ...) {
    ## define color pallet
    stopifnot(color_min < color_max)
    stopifnot(color_min < 0)
    stopifnot(length(color_scale) >= 3)

    # create a sequence from color_min to color max centered around 0
    n.col <- length(color_scale)
    zero_center_seq <- unique(c(
        seq(color_min, 0, length.out = ceiling(n.col / 2)),
        seq(0, color_max, length.out = ceiling(n.col / 2))
    ))

    if (!length(color_scale) == length(zero_center_seq)) {
        warning(sprintf(
            "Using %d/%d colors to center zero, dropping %s",
            length(zero_center_seq),
            n.col,
            color_scale[n.col]
        ), call. = FALSE)
        color_scale <- color_scale[seq(length(zero_center_seq))]
    }

    my.col <- circlize::colorRamp2(
        breaks = zero_center_seq,
        colors = color_scale
    )

    # ## query annotations on row
    if (!is.null(query_colors)) {
        stopifnot(all(rownames(cor_stats_layer) %in% names(query_colors)))
        query_colors <- query_colors[rownames(cor_stats_layer)]


        query_row_annotation <- ComplexHeatmap::rowAnnotation(
            " " = rownames(cor_stats_layer),
            col = list(" " = query_colors),
            show_legend = FALSE
        )
    } else {
        query_row_annotation <- NULL
    }

    ## reference annotation on bottom
    if (!is.null(reference_colors)) {
        stopifnot(all(colnames(cor_stats_layer) %in% names(reference_colors)))
        reference_colors <- reference_colors[colnames(cor_stats_layer)]

        ref_col_annotation <- ComplexHeatmap::columnAnnotation(
            " " = colnames(cor_stats_layer),
            col = list(" " = reference_colors),
            show_legend = FALSE
        )
    } else {
        ref_col_annotation <- NULL
    }

    ## add annotation
    if (!is.null(annotation)) {
        anno_matrix <- create_annotation_matrix(annotation, cor_stats_layer)
        print(anno_matrix)

        ## plot heatmap
        return(
            ComplexHeatmap::Heatmap(
                matrix = cor_stats_layer,
                col = my.col,
                name = "Cor",
                bottom_annotation = ref_col_annotation,
                right_annotation = query_row_annotation,
                cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(anno_matrix[i, j], x, y, gp = gpar(fontsize = 18))
                },
                ...
            )
        )
    }

    ## plot heatmap
    return(
        ComplexHeatmap::Heatmap(
            matrix = cor_stats_layer,
            col = my.col,
            name = "Cor",
            bottom_annotation = ref_col_annotation,
            right_annotation = query_row_annotation,
            column_names_gp = grid::gpar(fontsize=16),
            row_names_gp = grid::gpar(fontsize=16),
            ...
        )
    )
}

pdf(here("plots", "04_label_transfer", 
            "layer_enrichment_correlation_manual_allgenes.pdf"), width = 10, height = 10)
ht <- layer_stat_cor_plot(
    manual_cor,
    annotation=cor_annot, cluster_rows=FALSE, cluster_columns=FALSE
  )

draw(ht, row_title = "Xenium predicted domains", 
    column_title = "Manual annotations (Maynard, Collado-Torres et al)",
    row_title_gp = grid::gpar(fontsize = 18, fontface = "bold"),
    column_title_gp = grid::gpar(fontsize = 18, fontface = "bold")
    )
dev.off()


# ------------------------------------------------------------
# Try registration of spatialDLPFC layers to manaul annotations, using only
# the genes from the Xenium panel
# ------------------------------------------------------------

# Get the spatialDLPFC modeling results
# spatialDLPFC_modeling_results <- fetch_data(
#     type = "spatialDLPFC_Visium_modeling_results"
#   )

# spatialDLPFC_enrich <- spatialDLPFC_modeling_results$enrichment[, grep("^t_stat_", colnames(spatialDLPFC_modeling_results$enrichment))]
# spatialDLPFC_enrich <- spatialDLPFC_enrich[rownames(spatialDLPFC_enrich) %in% rownames(layer_res),]

# spatialDLPFC_cor <- layer_stat_cor(
#     spatialDLPFC_enrich,
#     manual_modeling_results,
#     model_type = "enrichment",
#     reverse = FALSE,
#     top_n = NULL
#   )

# # reorder to be diagonal
# spatialDLPFC_cor <- spatialDLPFC_cor[c("t_stat_Sp09D01", "t_stat_Sp09D02", "t_stat_Sp09D03", "t_stat_Sp09D05", "t_stat_Sp09D08", "t_stat_Sp09D04", "t_stat_Sp09D07", "t_stat_Sp09D06", "t_stat_Sp09D09"),]

# cor_annot <- annotate_registered_clusters(spatialDLPFC_cor)
# pdf(here("plots", "04_label_transfer", 
#             "layer_enrichment_correlation_spatialDLPFC_allgenes.pdf"), width = 10, height = 10)
# layer_stat_cor_plot(
#     spatialDLPFC_cor,
#     max = max(spatialDLPFC_cor)
#   )
# dev.off()


# ------------------------------------------------------------
# Register against the Visium layers from Boyi's dataset
# ------------------------------------------------------------

# vis_enrich <- readRDS("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/rds/04_SpD_marker_genes/test_enrich_PRECAST_07.rds")
# vis_use <- list()
# vis_use$enrichment <- vis_enrich

# vis_cor <- layer_stat_cor(
#     t_stats,
#     vis_use,
#     model_type = "enrichment",
#     reverse = FALSE,
#     top_n = NULL
#   )

# # reorder to be diagonal
# vis_cor <- vis_cor[colnames(vis_cor),]
# pdf(here("plots", "04_label_transfer", 
#             "layer_enrichment_correlation_visium_allgenes.pdf"), width = 10, height = 10)
# layer_stat_cor_plot(  
#     vis_cor,
#     max = max(vis_cor)
#   )
# dev.off()
