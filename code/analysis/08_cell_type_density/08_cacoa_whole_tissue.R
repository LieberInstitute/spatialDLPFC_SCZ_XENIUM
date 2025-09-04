suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
  library(escheR)
  library(readxl)
  library(edgeR)
  library(cacoa)
  library(cowplot)

})

####################################################################################################
# Try compositional data analysis of cell types using the cacoa package
# Reference: https://pklab.med.harvard.edu/viktor/cacoa/walkthrough_short.html
####################################################################################################


spe <- readRDS(here("processed-data", "07_cell_type_de", 
            "cleaned_spe_N24_with_cell_type_and_spds.RDS"))

spe$cell_types <- colData(spe)[["Banksy-clust_M0_lam0.1_k50_res0.7-cell-types"]]

brnums <- unique(spe$BrNum)

layers <- unique(spe$domain_annotations)

# Reference for creating new Cacoa object: https://github.com/kharchenkolab/cacoa
# Create a Cacoa object
pdf(here("plots", "08_cell_type_density", "whole_tissue_cacoa_cell_loadings.png"))
sample.groups <- colData(spe)[["Dx"]]
names(sample.groups) <- colData(spe)[["BrNum"]]

cell.groups <- colData(spe)[["cell_types"]]
names(cell.groups) <- spe$cell_id

sample.per.cell <- colData(spe)[["BrNum"]]
names(sample.per.cell) <- spe$cell_id

ref.level <- "NTC"
target.level <- "SCZ"

colnames(spe) <- spe$cell_id

cao <- Cacoa$new(
    as.matrix(counts(spe)), sample.groups=sample.groups, cell.groups=cell.groups, sample.per.cell=sample.per.cell, 
    ref.level=ref.level, target.level=target.level)


cao$estimateCellLoadings()
p <- cao$plotCellLoadings(show.pvals=TRUE)#+
#   theme(
#     axis.title.x = element_text(size = 20),  # X-axis label
#     axis.title.y = element_text(size = 20),  # Y-axis label
#     axis.text.x  = element_text(size = 18),  # tick labels
#     axis.text.y  = element_text(size = 18)
#   )

title <- ggdraw() + 
    draw_label("Whole tissue abundance analysis", fontface = 'bold')
print(plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)))



dev.off()
