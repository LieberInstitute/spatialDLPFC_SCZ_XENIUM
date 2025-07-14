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
  library(scales)
  library(ggridges)
  library(scattermore)
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

# read in the clusters and add them to the SPE
clusts <- read.csv(here("processed-data", "06_cell_type_clustering", sprintf("banksy_clustering_lambda%s_res%s.csv", lambda, res)))
clusts <- as.data.frame(clusts)
print(head(clusts))
colData(spe)[["Banksy"]] <- as.character(clusts$V1)

# Aggregate to broader groups and annotate
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

# Create the pseudobulk SPE for plotting
spe_pseudo <- aggregateAcrossCells(spe, ids=DataFrame(cbind(spe$annots, spe$BrNum)),
                                  coldata.merge=list(
                                    cell_area=sum,
                                    nucleus_area=sum))
spe_pseudo <- logNormCounts(spe_pseudo, size.factors=NULL)
colnames(spe_pseudo) <- spe_pseudo$annots


plot_counts <- as.matrix(logcounts(spe_pseudo)[grepl("NegControl",rownames(spe_pseudo)),])
plot_counts <- t(scale(t(plot_counts)))
colnames(plot_counts) <- colnames(spe_pseudo)

dend = cluster_between_groups(plot_counts, colnames(plot_counts))

plot_counts <- plot_counts[,order.dendrogram(dend)]

# Set colours of cell types
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
plot_counts_raw <- as.matrix(counts(spe_pseudo)[grepl("NegControl",rownames(spe_pseudo)),])/spe_pseudo$ncells
plot_counts_raw <- t(scale(t(plot_counts_raw)))
colnames(plot_counts_raw) <- colnames(spe_pseudo)

dend_raw = cluster_between_groups(plot_counts_raw, colnames(plot_counts_raw))

plot_counts_raw <- plot_counts_raw[,order.dendrogram(dend_raw)]

# Create a column annotation using cell type information
cell_types_raw <- colnames(plot_counts_raw)
#names(cell_types) <- colnames(spe_pseudo)

# Create the annotation
column_ha_raw <- HeatmapAnnotation(
  CellType = cell_types_raw,
  annotation_name_gp = gpar(fontsize = 10),
  col = list(CellType = celltype_colors)
)

plot_counts_cell_size_norm <- as.matrix(counts(spe_pseudo)[grepl("NegControl",rownames(spe_pseudo)),])/spe_pseudo$cell_area
plot_counts_cell_size_norm <- t(scale(t(plot_counts_cell_size_norm)))
colnames(plot_counts_cell_size_norm) <- colnames(spe_pseudo)

dend_cell_size_norm = cluster_between_groups(plot_counts_cell_size_norm, 
                colnames(plot_counts_cell_size_norm))
plot_counts_cell_size_norm <- plot_counts_cell_size_norm[,order.dendrogram(dend_cell_size_norm)]

cell_types_cell_size_norm <- colnames(plot_counts_cell_size_norm)
#names(cell_types_cell_size_norm) <- colnames(spe_pseudo)

# Create the annotation
column_ha_cell_size_norm <- HeatmapAnnotation(
  CellType = cell_types_cell_size_norm,
  annotation_name_gp = gpar(fontsize = 10),
  col = list(CellType = celltype_colors)
)

pdf(here("plots", "07_cell_type_de", "heatmap_ambig_clusters.pdf"))

ComplexHeatmap::Heatmap(plot_counts, cluster_columns=TRUE,
        name="logcounts",
         row_names_gp = gpar(fontsize = 10),
         column_names_gp=gpar(fontsize = 10), column_names_rot=45,
         bottom_annotation=column_ha)

ComplexHeatmap::Heatmap(plot_counts_raw, cluster_columns=dend_raw,
        name="counts divided\nby ncells", bottom_annotation = column_ha_raw, 
        row_names_gp = gpar(fontsize = 10), 
        show_column_names=TRUE)

ComplexHeatmap::Heatmap(plot_counts_cell_size_norm, cluster_columns=dend_cell_size_norm,
        name="counts divided\nby cell size", bottom_annotation = column_ha_cell_size_norm, 
        row_names_gp = gpar(fontsize = 10), 
        show_column_names=TRUE)
dev.off()

#---------------------------------------------------
# Plot the nucleus and cell area boxplots for 
# each of the cell types
#---------------------------------------------------

spe <- scuttle::addPerCellQCMetrics(spe, subsets=list(
  control_codeword_counts = grepl("NegControlCodeword", rownames(spe)),
  control_probe_counts = grepl("NegControlProbe", rownames(spe))
))
spe$detected_genes <- colSums(counts(spe)[rowData(spe)$Type == "Gene Expression",] > 0)
spe$gene_counts <-  colSums(counts(spe)[rowData(spe)$Type == "Gene Expression",])

area_df <- colData(spe) %>%
  as.data.frame() %>%
  select(cell_area, nucleus_area, annots, 
  subsets_control_probe_counts_percent, subsets_control_codeword_counts_percent, 
  detected, detected_genes, total_counts, gene_counts,
  subsets_control_codeword_counts_sum, subsets_control_probe_counts_sum)%>%
  mutate(is_neuron = ifelse(annots %in% c("L2/3 Ex", "L4/5 Ex", "L5 Ex", "L6 Ex", "In: SST, PVALB", "In: VIP, LAMP5"), TRUE, FALSE))


pdf(here("plots", "07_cell_type_de", "cell_nucleus_area_boxplots.pdf"), width=10, height=10)
# ggplot(area_df, aes(x=annots, y=cell_area, colour=annots))+
#   geom_violin()+
#   geom_boxplot(width=.1)+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   ylab("Cell area")+
#   xlab("Cell type")+
#   scale_color_manual(values=celltype_colors)+
#   ggtitle("Cell area by cell type")

# ggplot(area_df, aes(x=annots, y=nucleus_area, colour=annots))+
#   geom_violin()+
#   geom_boxplot(width=.1)+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   ylab("Nucleus area")+
#   xlab("Cell type")+
#   scale_color_manual(values=celltype_colors)+
#   ggtitle("Nucleus area by cell type")

ggplot(area_df, aes(x=nucleus_area, y=annots, group=annots))+
    geom_density_ridges(aes(fill=is_neuron), alpha=0.8, scale=1.5)+
    scale_fill_discrete()+
    ggtitle("Cell type specific nucleus area")

ggplot(area_df, aes(x=cell_area, y=annots, group=annots))+
    geom_density_ridges(aes(fill=is_neuron), alpha=0.8, scale=1.5)+
    scale_fill_discrete()+
    ggtitle("Cell type specific cell area")

ggplot(area_df, aes(x=total_counts, y=annots, group=annots))+
    geom_density_ridges(aes(fill=is_neuron), alpha=0.8, scale=1.5)+
    scale_fill_discrete()+
    ggtitle("Cell type specific total counts")

ggplot(area_df, aes(x=gene_counts, y=annots, group=annots))+
    geom_density_ridges(aes(fill=is_neuron), alpha=0.8, scale=1.5)+
    scale_fill_discrete()+
    ggtitle("Cell type specific gene counts")


ggplot(area_df, aes(x=detected_genes, y=annots, group=annots))+
    geom_density_ridges(aes(fill=is_neuron), alpha=0.8, scale=1.5)+
    scale_fill_discrete()+
    ggtitle("Cell type specific number of detected genes")

ggplot(area_df, aes(x=subsets_control_codeword_counts_sum, y=annots, group=annots))+
    geom_density_ridges(aes(fill=is_neuron), alpha=0.8, scale=1.5)+
    scale_fill_discrete()+
    ggtitle("Cell type specific number of control codewords")

ggplot(area_df, aes(x=subsets_control_probe_counts_sum, y=annots, group=annots))+
    geom_density_ridges(aes(fill=is_neuron), alpha=0.8, scale=1.5)+
    scale_fill_discrete()+
    ggtitle("Cell type specific number of control probes")


ggplot(area_df, aes(x=subsets_control_codeword_counts_percent, y=annots, group=annots))+
    geom_density_ridges(aes(fill=is_neuron), alpha=0.8, scale=1.5)+
    scale_fill_discrete()+
    ggtitle("Cell type specific percentage of control codewords")

ggplot(area_df, aes(x=subsets_control_probe_counts_percent, y=annots, group=annots))+
    geom_density_ridges(aes(fill=is_neuron), alpha=0.8, scale=1.5)+
    scale_fill_discrete()+
    ggtitle("Cell type specific percentage of control probes")

  ggplot(area_df, aes(x=subsets_control_probe_counts_sum/cell_area, y=annots, group=annots))+
      geom_density_ridges(aes(fill=is_neuron), alpha=0.8, scale=1.5)+
      scale_fill_discrete()+
      ggtitle("Cell type specific density of control probes")
  
  ggplot(area_df, aes(x=subsets_control_codeword_counts_sum/cell_area, y=annots, group=annots))+
      geom_density_ridges(aes(fill=is_neuron), alpha=0.8, scale=1.5)+
      scale_fill_discrete()+
      ggtitle("Cell type specific density of control codewords")
  
  ggplot(area_df, aes(x=detected_genes/cell_area, y=annots, group=annots))+
      geom_density_ridges(aes(fill=is_neuron), alpha=0.8, scale=1.5)+
      scale_fill_discrete()+
      ggtitle("Number of genes detected per micron squared")

# ggplot(area_df, aes(x=total_counts, y=subsets_control_probe_counts_percent))+
#   geom_point(aes(colour=is_neuron), alpha=0.5)+
#   #geom_smooth(method = "lm", se=FALSE)+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   xlab("Total counts")+
#   ylab("Percentage of control probe counts")+
#   facet_wrap(~annots)+
#   geom_scattermore()

# ggplot(area_df, aes(x=total_counts, y=subsets_control_codeword_counts_percent))+
#   geom_point(aes(colour=is_neuron), alpha=0.5)+
#   #geom_smooth(method = "lm", se=FALSE)+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   xlab("Total counts")+
#   ylab("Percentage of control probe counts")+
#   facet_wrap(~annots)+
#   geom_scattermore()

#   ggplot(area_df, aes(x=detected_genes, y=subsets_control_codeword_counts_sum/cell_area))+
#     geom_point(aes(colour=is_neuron), alpha=0.5)+
#     #geom_smooth(method = "lm", se=FALSE)+
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#     xlab("Detected")+
#     ylab("Control codeword counts divided by cell area")+
#     facet_wrap(~annots)+
#     geom_scattermore()
 
  ggplot(area_df, aes(x=detected_genes, y=subsets_control_probe_counts_sum/cell_area))+
    geom_point(aes(colour=is_neuron), alpha=0.5)+
    #geom_smooth(method = "lm", se=FALSE)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("Detected")+
    ylab("Control probe counts divided by cell area")+
    facet_wrap(~annots)+
    geom_scattermore()

  ggplot(area_df, aes(x=cell_area, y=subsets_control_probe_counts_sum))+
    geom_point(aes(colour=is_neuron), alpha=0.5)+
    #geom_smooth(method = "lm", se=FALSE)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("Cell area")+
    ylab("Control probe counts")+
    geom_scattermore()+
    geom_smooth(method="lm", se=FALSE,
                  data = subset(area_df, cell_area<600))

  ggplot(area_df, aes(x=cell_area, y=subsets_control_codeword_counts_sum))+
    geom_point(aes(colour=is_neuron), alpha=0.5)+
    #geom_smooth(method = "lm", se=FALSE)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("Cell area")+
    ylab("Control codeword counts")+
    geom_scattermore()+
    geom_smooth(method="lm", se=FALSE,
                  data = subset(area_df, cell_area<600))
# ggplot(area_df, aes(x=annots, y=subsets_control_codeword_counts_percent, colour=annots))+
#   geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.01))+
#   ylab("Control codeword counts")+
#   xlab("Cell type")+
#   scale_color_manual(values=celltype_colors)+
#   ggtitle("Control codeword counts by cell type")

# ggplot(area_df, aes(x=annots, y=subsets_control_probe_counts_percent, colour=annots))+
#   geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.01))+
#   ylab("Control probe counts")+
#   xlab("Cell type")+
#   scale_color_manual(values=celltype_colors)+
#   ggtitle("Control probe counts by cell type")

# ggplot(area_df, aes(x=nucleus_area, y=subsets_control_codeword_counts_percent, colour=is_neuron))+
#   geom_point()+
#   geom_smooth(method = "lm", se=FALSE)+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   xlab("Nucleus area")+
#   ylab("Percentage of control codeword counts")+
#   facet_wrap(~is_neuron)

# ggplot(area_df, aes(x=nucleus_area, y=subsets_control_probe_counts_percent, colour=is_neuron))+
#   geom_point()+
#   geom_smooth(method = "lm", se=FALSE)+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   xlab("Nucleus area")+
#   ylab("Percentage of control probe counts")+
#   facet_wrap(~is_neuron)

# ggplot(area_df, aes(x=nucleus_area, y=subsets_control_probe_counts_percent, colour=is_neuron))+
#   geom_point()+
#   geom_smooth(method = "lm", se=FALSE)+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   xlab("Nucleus area")+
#   ylab("Percentage of control probe counts")+
#   facet_wrap(~annots)

# ggplot(area_df, aes(x=nucleus_area, y=subsets_control_codeword_counts_percent, colour=is_neuron))+
#   geom_point()+
#   geom_smooth(method = "lm", se=FALSE)+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   xlab("Nucleus area")+
#   ylab("Percentage of control codeword counts")+
#   facet_wrap(~annots)

# ggplot(area_df, aes(x=detected_genes, y=subsets_control_probe_counts_percent))+
#   geom_point()+
#   facet_wrap(~annots)+
#   geom_smooth(method = "lm", se=FALSE)+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   xlab("Number of detected genes")+
#   ylab("Percentage of control probe counts")

# ggplot(area_df, aes(x=detected_genes, y=subsets_control_codeword_counts_percent))+
#   geom_point()+
#   facet_wrap(~annots)+
#   geom_smooth(method = "lm", se=FALSE)+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   xlab("Number of detected genes")+
#   ylab("Percentage of control codeword counts")
dev.off()

# proportion of neurons with 0% control codeword counts
area_df %>%
  filter(is_neuron) %>%
  summarise(prop_zero = sum(subsets_control_codeword_counts_percent == 0)/n()) %>%
  print()

# proportion of non-neurons with 0% control codeword counts
area_df %>%
  filter(!is_neuron) %>%
  summarise(prop_zero = sum(subsets_control_codeword_counts_percent == 0)/n()) %>%
  print()

# proportion of neurons with 0% control probe counts
area_df %>%
  filter(is_neuron) %>%
  summarise(prop_zero = sum(subsets_control_probe_counts_percent == 0)/n()) %>%
  print()
# proportion of non-neurons with 0% control probe counts
area_df %>%
  filter(!is_neuron) %>%
  summarise(prop_zero = sum(subsets_control_probe_counts_percent == 0)/n()) %>%
  print()





