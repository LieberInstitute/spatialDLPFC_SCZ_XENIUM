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
# Compute and plot the proportion of each cell type in the spatial domains
####################################################################################################


spe <- readRDS(here("processed-data", "07_cell_type_de", 
            "cleaned_spe_N24_with_cell_type_and_spds.RDS"))

spe$cell_types <- colData(spe)[["Banksy-clust_M0_lam0.1_k50_res0.7-cell-types"]]

brnums <- unique(spe$BrNum)
diags <- unique(spe$Dx)

all_tabs <- list()
for (i in 1:length(diags)){
    dx_use <- diags[[i]]
    spe_use <- spe[, spe$Dx == dx_use]


    # number of cell types in each spatial domain
    tab <- table(spe_use$spatransfer_k50_predictions_smooth, spe_use$cell_types)
    tab <- tab/rowSums(tab) # proportions in each spatial domain, sums to 1 across each row
    tab <- as.data.frame(tab)
    tab$Dx <- unique(spe_use$Dx)

    all_tabs[[i]] <- tab
    rm(tab)
    gc()

}

rm(spe)
gc()


# ntc_tab <- all_tabs[[1]]
# colnames(ntc_tab)  <- c("spatial_domain", "cell_type", "proportion", "Dx")
# ntc_tab <- ntc_tab %>%
#     pivot_wider(names_from = spatial_domain, values_from = proportion)

tabs_df <- do.call(rbind, all_tabs)
colnames(tabs_df) <- c("spatial_domain", "cell_type", "proportion", "Dx")
tabs_df_wide <- tabs_df %>%
    pivot_wider(names_from = cell_type, values_from = proportion)
tabs_df_wider <- tabs_df_wide %>%
    pivot_wider(values_from=contains("`"),
                names_from = Dx) 

# plot a heatmap of the proportions
pdf(here("plots", "07_cell_type_de", "cell_type_props_in_spds_heatmap.pdf"), height=10, width=15)
# Reorder BrNum so NTC bars are on the left
tabs_df_wide <- tabs_df_wide %>%
  mutate(Dx_ordered = factor(Dx, levels = unique(Dx[order(Dx)])))

plot_mat <- tabs_df_wide %>%
    select(starts_with("spd")) %>%
    as.data.frame()%>%
    as.matrix()
Heatmap(plot_mat)
dev.off()
