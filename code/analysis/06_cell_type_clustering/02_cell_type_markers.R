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
spe <- logNormCounts(spe)

markers <- findMarkers(logcounts(spe), groups=colData(spe)[["Banksy"]],
                        test.type="t", block=colData(spe)[["BrNum"]])

top_markers <- lapply(markers, function(x){return(rownames(x)[1:20])})
top_markers <- unique(do.call(c, top_markers))

spe_pseudo <- aggregateAcrossCells(spe, ids=DataFrame(cbind(spe$Banksy, spe$BrNum)))
spe_pseudo <- logNormCounts(spe_pseudo, size.factors=NULL)
spe_pseudo <- scater::runPCA(spe_pseudo, ncomponents=50)

pdf(here::here("plots", "06_cell_type_clustering", sprintf("banksy_lambda%s_res%s_PCA.pdf", lambda, res)))
scater::plotPCA(spe_pseudo, colour_by = "Dx", ncomponents=5)
scater::plotPCA(spe_pseudo, colour_by = "Banksy", shape_by="Dx", ncomponents=5)
scater::plotPCA(spe_pseudo, colour_by = "BrNum", ncomponents=5)
dev.off()

# Sangho's markers
panel_markers <- readxl::read_xlsx((here("raw-data", 
        "experiment_info", 
        "Xenium_SHK_celltype_REannot_2025-04-13.xlsx")), sheet=2)

panel_markers <- panel_markers %>%
    as.data.frame() %>%
    mutate(cell_type_updated=case_when(cell_type_updated=="NA" ~ NA,
                                        TRUE ~ cell_type_updated))
    #  mutate(cell_type = case_when(cell_type_updated == 'BDNF' ~ NA,
    #                        cell_type == "In (EI balance)" ~ 'In',
    #                        TRUE ~ cell_type))

markers_plot <- panel_markers %>%
    filter(!is.na(cell_type_updated))

markers_plot <- markers_plot[order(markers_plot$cell_type_updated),]

pdf(here::here("plots", "06_cell_type_clustering", sprintf("banksy_lambda%s_res%s_cell_types_markers.pdf", lambda, res)),
            height=15, width=28)
plot_counts <- as.matrix(logcounts(spe_pseudo)[rownames(spe_pseudo) %in% top_markers,])
plot_counts <- plot_counts[top_markers,]
plot_counts <- t(scale(t(plot_counts)))
dend = cluster_between_groups(plot_counts, spe_pseudo$Banksy)

ha = HeatmapAnnotation(Banksy_label = anno_text(spe_pseudo$Banksy))
ComplexHeatmap::Heatmap(plot_counts, name="counts", bottom_annotation=ha, 
            cluster_columns = dend, row_names_gp = gpar(fontsize = 16))

plot_counts <- as.matrix(logcounts(spe_pseudo)[rownames(spe_pseudo) %in% markers_plot$Gene,])
plot_counts <- plot_counts[markers_plot$Gene, ]
plot_counts <- t(scale(t(plot_counts)))
dend = cluster_between_groups(plot_counts, spe_pseudo$Banksy)
ha = HeatmapAnnotation(Banksy_label = anno_text(spe_pseudo$Banksy, rot=45),
            gp = gpar(fontsize = 6))


#row_annotation = rowAnnotation(cell_type = markers_plot$cell_type)
row_annotation = rowAnnotation(foo = anno_text(markers_plot$cell_type_updated, location = 0.5, just = "center",
    gp = gpar(fill = plyr::mapvalues(markers_plot$cell_type_updated, 
                from=unique(markers_plot$cell_type_updated), to=1:length(unique(markers_plot$cell_type_updated))), col="white"),
    width = max_text_width(markers_plot$cell_type_updated)*1.2))

ComplexHeatmap::Heatmap(plot_counts, name="counts", bottom_annotation=ha, 
            cluster_columns = dend, row_names_gp = gpar(fontsize = 16), 
            right_annotation = row_annotation, cluster_rows=FALSE)


dev.off()







