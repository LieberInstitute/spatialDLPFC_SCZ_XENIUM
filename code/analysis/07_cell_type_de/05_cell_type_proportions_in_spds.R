suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
})

####################################################################################################
# Compute and plot the proportion of each cell type in the spatial domains
####################################################################################################


spe <- readRDS(here("processed-data", "07_cell_type_de", 
            "cleaned_spe_N24_with_cell_type_and_spds.RDS"))

spe$cell_types <- colData(spe)[["Banksy-clust_M0_lam0.1_k50_res0.7-cell-types"]]

brnums <- unique(spe$BrNum)

all_tabs <- list()
for (i in 1:length(brnums)){
    br_use <- brnums[[i]]
    spe_use <- spe[, spe$BrNum == br_use]


    # number of cell types in each spatial domain
    tab <- table(spe_use$spatransfer_k50_predictions_smooth, spe_use$cell_types)
    tab <- tab/rowSums(tab) # proportions in each spatial domain, sums to 1 across each row
    tab <- as.data.frame(tab)
    tab$BrNum <- br_use
    tab$Dx <- unique(spe_use$Dx)

    all_tabs[[i]] <- tab
    rm(tab)
    gc()

}

rm(spe)
gc()

tabs_df <- do.call(rbind, all_tabs)
pdf(here("plots", "07_cell_type_de", "cell_type_props_in_spds.pdf"), height=10, width=15)
# Reorder BrNum so NTC bars are on the left
tabs_df <- tabs_df %>%
  mutate(BrNum_ordered = factor(BrNum, levels = unique(BrNum[order(Dx)])))

# Create rectangle data: one row per BrNum
rect_data <- tabs_df %>%
  distinct(BrNum_ordered, Dx) %>%
  mutate(
    xmin = as.numeric(BrNum_ordered) - 0.5,
    xmax = as.numeric(BrNum_ordered) + 0.5,
    ymin = -0.08,
    ymax = -0.02  # rectangle height below x-axis
  )

# Plot with rectangles
ggplot(tabs_df, aes(x = BrNum_ordered, y = Freq, fill = Var2)) +
  geom_bar(position = "fill", stat = "identity", aes(fill=Var2)) +
  facet_wrap(~Var1) +
  # Dx annotation rectangles
  geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Dx),
            inherit.aes = FALSE, colour = NA) +
  scale_fill_manual(values = c("NTC" = "skyblue", "SCZ" = "salmon")) +  # tweak colors here
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    plot.margin = margin(10, 10, 40, 10),
    legend.position = "right"
  ) +
  coord_cartesian(ylim = c(-0.1, 1))  # makes space for rectangles below
dev.off()
