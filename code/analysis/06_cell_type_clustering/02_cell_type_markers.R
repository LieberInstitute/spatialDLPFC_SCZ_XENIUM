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


pdf(here("plots", "06_cell_type_clustering", 
            sprintf("banksy_lambda%s_res%s_cell_types_markers.pdf", lambda, res)),
            height=15, width=25)
plot_counts <- as.matrix(counts(spe_pseudo)[rownames(spe_pseudo) %in% top_markers,])
dend = cluster_between_groups(plot_counts, spe_pseudo$Banksy)
ha = HeatmapAnnotation(Banksy_label = spe_pseudo$Banksy)
ComplexHeatmap::Heatmap(plot_counts, name="counts", bottom_annotation=ha, 
            cluster_columns = dend, row_names_gp = gpar(fontsize = 16))
dev.off()

