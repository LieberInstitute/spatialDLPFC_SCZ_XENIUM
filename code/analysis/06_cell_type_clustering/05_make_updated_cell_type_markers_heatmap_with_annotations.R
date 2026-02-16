library(SummarizedExperiment)
library(SpatialExperiment)
library(Seurat)
library(scater)
library(cowplot)
library(ggplot2)
library(here)
library(tidyverse)
library(escheR)
library(scran)
library(scater)
library(ComplexHeatmap)
library(circlize)

##########################################################################################
# Run marker gene detection on BANKSY clusters
##########################################################################################

# Read in the data
spe <- readRDS(here("processed-data", "01_build_spe", "raw_spe_N24.RDS"))
outlier_ids <- read.csv(here("processed-data", "02_xenium_qc", "outlier_ids.csv"))$x
spe <- spe[, -which(colnames(spe) %in% outlier_ids)]
rm (outlier_ids)
spe <- spe[which(rowData(spe)$Type=="Gene Expression"), ]
print(spe)

lambda <- 0.1
res <- 0.7

clusts <- read.csv(here("processed-data", "06_cell_type_clustering", sprintf("banksy_clustering_lambda%s_res%s.csv", lambda, res)))
clusts <- as.data.frame(clusts)
print(head(clusts))
colData(spe)[["Banksy"]] <- as.character(clusts$V1)

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
                          V1 %in% c(12) ~ "MGE",
                          V1 %in% c(9) ~ "CGE",
                       TRUE ~ "NA")) %>%
    mutate(annots_combined = paste(annots, V1, sep="."))
  
colData(spe)[["annots"]] <- clusts$annots

spe <- logNormCounts(spe)

markers <- findMarkers(logcounts(spe), groups=colData(spe)[["Banksy"]],
                        test.type="t", block=colData(spe)[["BrNum"]])

top_markers <- lapply(markers, function(x){return(rownames(x)[1:20])})
top_markers <- unique(do.call(c, top_markers))

spe_pseudo <- aggregateAcrossCells(spe, ids=DataFrame(cbind(spe$annots, spe$BrNum)))
spe_pseudo <- logNormCounts(spe_pseudo, size.factors=NULL)
spe_pseudo <- scater::runPCA(spe_pseudo, ncomponents=50)

# Sangho's markers
panel_markers <- readxl::read_xlsx((here("raw-data", 
        "experiment_info", 
        "Xenium_SHK_celltype_REannot_2025-04-13.xlsx")), sheet=2)

panel_markers <- panel_markers %>%
    as.data.frame() %>%
    mutate(cell_type_updated=case_when(cell_type_updated=="NA" ~ NA,
                                        TRUE ~ cell_type_updated))

markers_plot <- panel_markers %>%
    filter(!is.na(cell_type_updated))

markers_plot <- markers_plot[order(markers_plot$cell_type_updated),]

pdf(here::here("plots", "06_cell_type_clustering", sprintf("annotated_banksy_lambda%s_res%s_cell_types_markers.pdf", lambda, res)),
            height=15, width=28)
plot_counts <- as.matrix(logcounts(spe_pseudo)[rownames(spe_pseudo) %in% markers_plot$Gene,])
plot_counts <- plot_counts[markers_plot$Gene,]
plot_counts <- t(scale(t(plot_counts)))
dend = cluster_between_groups(plot_counts, spe_pseudo$annots)

#row_annotation = rowAnnotation(cell_type = markers_plot$cell_type)
row_annotation = rowAnnotation(foo = anno_text(markers_plot$cell_type_updated, location = 0.5, just = "center",
    gp = gpar(fill = plyr::mapvalues(markers_plot$cell_type_updated, 
                from=unique(markers_plot$cell_type_updated), to=1:length(unique(markers_plot$cell_type_updated))), col="white"),
    width = max_text_width(markers_plot$cell_type_updated)*1.2))

ha = HeatmapAnnotation(Banksy_label = anno_text(spe_pseudo$annots))
ComplexHeatmap::Heatmap(plot_counts, name="counts", bottom_annotation=ha, 
            cluster_columns = dend, row_names_gp = gpar(fontsize = 16), 
            right_annotation=row_annotation,
            cluster_rows=FALSE)

plot_counts <- as.matrix(logcounts(spe_pseudo)[rownames(spe_pseudo) %in% markers_plot$Gene,])
colnames(plot_counts) <- paste(spe_pseudo$BrNum, spe_pseudo$Dx, sep="_")
plot_counts <- plot_counts[markers_plot$Gene, ]
plot_counts <- t(scale(t(plot_counts)))
dend = cluster_between_groups(plot_counts, spe_pseudo$annots)
ha = HeatmapAnnotation(Banksy_label = anno_text(spe_pseudo$annots, rot=45),
            gp = gpar(fontsize = 6))



# 1. Extract Dx group from column names
Dx_group <- ifelse(grepl("NTC$", colnames(plot_counts)), "NTC", "SCZ")

# 2. Create a factor to preserve order in legend
Dx_group <- factor(Dx_group, levels = c("NTC", "SCZ"))

# 3. Define colors
dx_colors <- c("NTC" ="steelblue", "SCZ" = "firebrick")

# 4. Create bottom annotation
col_ha_bottom <- HeatmapAnnotation(
  Dx = Dx_group,
  col = list(Dx = dx_colors),
  annotation_name_side = "left"
)

# 5. Draw heatmap with bottom annotation
ComplexHeatmap::Heatmap(plot_counts, name="counts", bottom_annotation=col_ha_bottom, 
            cluster_columns = dend, row_names_gp = gpar(fontsize = 16),
            right_annotation = row_annotation, cluster_rows=FALSE)

dev.off()








