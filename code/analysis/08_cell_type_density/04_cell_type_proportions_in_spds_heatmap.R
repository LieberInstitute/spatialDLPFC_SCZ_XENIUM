suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
  library(ComplexHeatmap)
  library(circlize)
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
    tab <- table(spe_use$domain_annotations, spe_use$cell_types)
    tab <- tab/rowSums(tab) # proportions in each spatial domain, sums to 1 across each row
    tab <- as.data.frame(tab)
    tab$Dx <- unique(spe_use$Dx)

    all_tabs[[i]] <- tab
    rm(tab)
    gc()

}

rm(spe)
gc()

ntc_tab <- all_tabs[[1]]
colnames(ntc_tab)  <- c("spatial_domain", "cell_type", "proportion", "Dx")

cell_types_unique<- unique(ntc_tab$cell_type)

ntc_tab <- ntc_tab %>%
    pivot_wider(names_from = cell_type, values_from = proportion)%>%
    rename_at(vars(all_of(cell_types_unique)),function(x) paste0(x,"_NTC"))


scz_tab <- all_tabs[[2]]
colnames(scz_tab)  <- c("spatial_domain", "cell_type", "proportion", "Dx")
scz_tab <- scz_tab %>%
    pivot_wider(names_from = cell_type, values_from = proportion)%>%
    rename_at(vars(all_of(cell_types_unique)),function(x) paste0(x,"_SCZ"))


all_tab <- cbind(ntc_tab, scz_tab)

# plot a heatmap of the proportions

all_tab_plot <- all_tab %>%
  select(-c(spatial_domain, Dx)) %>%
  as.matrix()
rownames(all_tab_plot) <- all_tab$spatial_domain

all_tab_plot <- all_tab_plot[, order(colnames(all_tab_plot))] # sort the columns by name

column_condition <- ifelse(grepl("NTC$", colnames(all_tab_plot)), "NTC", "SCZ")
col_ha <- HeatmapAnnotation(
  Dx = column_condition,
  col = list(Dx = c("NTC" = "skyblue", "SCZ" = "tomato"))
)

# set colour limits
col_fun <- colorRamp2(c(0, 0.5, 1), c("blue","white", "red"))

pdf(here("plots", "07_cell_type_de", "cell_type_props_in_spds_heatmap.pdf"), height=10, width=15)


Heatmap(all_tab_plot, 
          bottom_annotation = col_ha,
          cluster_columns = FALSE,
          cluster_rows=FALSE,
          col = col_fun,
          name="Proportion")
dev.off()