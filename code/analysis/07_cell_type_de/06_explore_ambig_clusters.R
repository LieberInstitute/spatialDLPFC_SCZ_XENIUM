suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
  library(ComplexHeatmap)
  library(grid)
})


# Plot a heatmap of the neg ctrl features for each pseudobulk sample
spe <- readRDS(here("processed-data", "01_build_spe", "raw_spe_N24.RDS"))
outlier_ids <- read.csv(here("processed-data", "02_xenium_qc", "outlier_ids.csv"))$x
spe <- spe[, -which(colnames(spe) %in% outlier_ids)]
rm (outlier_ids)
print(spe)

#### Banksy parameters ####
lambda <- 0.1
res <- 0.7 # higher = more clusters
use_agf <- FALSE
cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(use_agf), lambda, res)



clusts <- read.csv(here("processed-data", "06_cell_type_clustering", sprintf("banksy_clustering_lambda%s_res%s.csv", lambda, res)))
clusts <- as.data.frame(clusts)
print(head(clusts))
colData(spe)[["Banksy"]] <- as.character(clusts$V1)


# Aggregate to broader groups

clusts <- clusts %>% as.data.frame() %>% 
  mutate(V1 = as.numeric(V1)) %>%
  mutate(annots=case_when(V1 %in% c(6, 8, 1) ~ "Oligo", # 10 and 15 i am not super sure about
                          V1 %in% c(10) ~ "Ambig/Oligo",
                          V1 %in% c(7) ~ "Mic",
                          V1 %in% c(15) ~ "Ambig/In/Endo",
                          V1 %in% c(18) ~ "L5 Ex",
                          V1 %in% c(14) ~ "L6 Ex",
                          V1 %in% c(11) ~ "L4/5 Ex",
                          V1 %in% c(2) ~ "L2/3 Ex",
                          V1 %in% c(17, 16, 4) ~ "Ast",
                          V1 %in% c(5, 13, 3) ~ "Endo",
                          V1 %in% c(12) ~ "In: VIP, LAMP5",
                          V1 %in% c(9) ~ "In: SST, PVALB",
                       TRUE ~ "NA")) %>%
    mutate(annots_combined = paste(annots, V1, sep="."))
  
colData(spe)[["annots"]] <- clusts$annots
colData(spe)[["annots_combined"]] <- clusts$annots


# Add some slide and date related metadata
slide_id <- unlist(lapply(strsplit(spe$Sample, split="/"), "[", 9))
slide_id <- unlist(lapply(strsplit(slide_id, split="__Br"), "[", 1))
spe$slide_id <- slide_id

run_date <- unlist(lapply(strsplit(spe$Sample, split="/"), "[", 7))
spe$run_date <- run_date


## First, create pseudobulk data:

# spe_pseudo <- registration_pseudobulk(
#       spe,
#       var_registration = "annots",
#       var_sample_id = "BrNum",
#       covars = c("Dx", "Age", "Sex", "slide_id", "run_date"),
#       min_ncells = 10,
#       pseudobulk_rds_file = here(
#         "processed-data", "07_cell_type_de",
#         paste0(sprintf("spe_pseudo_donor_celltype_with_negctrls_%s", cnm), "", ".rds")
#       )
#     )

spe_pseudo <- aggregateAcrossCells(spe, ids=DataFrame(cbind(spe$annots, spe$BrNum)))
spe_pseudo <- logNormCounts(spe_pseudo, size.factors=NULL)
colnames(spe_pseudo) <- spe_pseudo$annots


plot_counts <- as.matrix(logcounts(spe_pseudo)[grepl("NegControl",rownames(spe_pseudo)),])
plot_counts <- t(scale(t(plot_counts)))
colnames(plot_counts) <- colnames(spe_pseudo)

dend = cluster_between_groups(plot_counts, colnames(plot_counts))

plot_counts <- plot_counts[,order.dendrogram(dend)]

# Set colours of cell types
# Define custom colors
celltype_colors <- c(
  "Ast" = "#A6CEE3",
  "Endo" = "#1F78B4",
  "Mic" = "#B2DF8A",
  "Oligo" = "#33A02C",
  "L2/3 Ex" = "#FB9A99",
  "L4/5 Ex" = "#E31A1C",
  "L5 Ex" = "#FDBF6F",
  "L6 Ex" = "#FF7F00",
  "In: SST, PVALB" = "#CAB2D6",
  "In: VIP, LAMP5" = "#6A3D9A",
  "Ambig/In/Endo" = "#B15928",
  "Ambig/Oligo" = "#FF69B4"
)

# Create a column annotation using cell type information
cell_types <- colnames(plot_counts)
# Create the colour annotation
column_ha <- HeatmapAnnotation(
  CellType = cell_types,
  annotation_name_gp = gpar(fontsize = 10),
  col = list(CellType = celltype_colors)
)
# Plot the raw counts too
plot_counts_raw <- as.matrix(counts(spe_pseudo)[grepl("NegControl",rownames(spe_pseudo)),])
plot_counts_raw <- t(scale(t(plot_counts_raw)))
colnames(plot_counts_raw) <- colnames(spe_pseudo)

dend_raw = cluster_between_groups(plot_counts_raw, colnames(plot_counts_raw))

plot_counts_raw <- plot_counts[,order.dendrogram(dend_raw)]

# Create a column annotation using cell type information
cell_types_raw <- colnames(plot_counts_raw)
#names(cell_types) <- colnames(spe_pseudo)

# Create the annotation
column_ha_raw <- HeatmapAnnotation(
  CellType = cell_types_raw,
  annotation_name_gp = gpar(fontsize = 10),
  col = list(CellType = celltype_colors)
)

pdf(here("plots", "07_cell_type_de", "heatmap_ambig_clusters.pdf"))

ComplexHeatmap::Heatmap(plot_counts, cluster_columns=dend,
        name="logcounts",
         row_names_gp = gpar(fontsize = 10),
         column_names_gp=gpar(fontsize = 10), column_names_rot=45,
         bottom_annotation=column_ha)

ComplexHeatmap::Heatmap(plot_counts_raw, cluster_columns=dend_raw,
        name="counts", bottom_annotation = column_ha_raw, 
        row_names_gp = gpar(fontsize = 10), 
        show_column_names=TRUE)

dev.off()











