suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
  library(escheR)
  library(patchwork)
  library(cowplot)
  library(ggplot2)
  library(ggpubr)
})

####################################################################################################
# Pseudobulk the SPE to the donor-domain level using the smoothed labels from label transfer.
# Reference: https://github.com/LieberInstitute/spatialDLPFC_SCZ/blob/main/code/analysis/pseudobulk_dx/create_pb_data.R
####################################################################################################

spe <- readRDS(here("processed-data", "01_build_spe", "raw_spe_N24.RDS"))
outlier_ids <- read.csv(here("processed-data", "02_xenium_qc", "outlier_ids.csv"))$x
spe <- spe[, -which(colnames(spe) %in% outlier_ids)]
rm(outlier_ids)
spds <- read.csv(here("processed-data", "04_label_transfer", "label_transfer_N24_k50_smoothed_labels.csv"))

rownames(spds) <- spds$X 
spds$X <- NULL
colnames(spds) <- c("predictions_smooth")
stopifnot(all(colnames(spe)==rownames(spds)))
colData(spe)$predictions_smooth <- factor(spds$predictions_smooth)

brs_use <- c("Br8667", "Br5973")

spe <- spe[,spe$BrNum %in% brs_use]


# map spd style labels to actual annotations
domain_annotations <- colData(spe) %>%
  as.data.frame() %>%
  mutate(domain_annotations = case_when(predictions_smooth == "spd07" ~ "L1/M",
                                        predictions_smooth == "spd06" ~ "L2/3",
                                        predictions_smooth == "spd02" ~ "L3/4",
                                        predictions_smooth == "spd05" ~ "L5",
                                        predictions_smooth == "spd03" ~ "L6",
                                        predictions_smooth == "spd01" ~ "WMtz",
                                        predictions_smooth == "spd04" ~ "WM",
                                          TRUE ~ NA))

pal <- set_names(Polychrome::palette36.colors(7)[seq.int(7)],
c("WMtz", "L3/4", "L6", "WM", "L5", "L2/3", "L1/M"))

colData(spe)$domain_annotations <- domain_annotations$domain_annotations

plist <- list()
for (i in 1:length(brs_use)){
    spe_use <- spe[, spe$BrNum == brs_use[i]]

    p <- make_escheR(spe_use, y_reverse=FALSE) %>%
        add_fill("domain_annotations") +
        ggtitle(paste(brs_use[i], unique(spe_use$Dx), sep=" "))+
        scale_fill_manual(
        name = "Spatial Domain",
        values = pal,
        guide = guide_legend(override.aes = list(size = 7)))+
        theme(legend.direction = "horizontal")
         

    plist[[i]] <- p

}

png(here("plots", "04_label_transfer", "transferred_spds_on_tissue.png"),
     width=2000, height=3500, res=300)
wrap_plots(plist, ncol=1)+
    plot_layout(guides="collect")&
    theme(legend.position="none")
dev.off()


# plot the legend separately
p_legend <- get_legend(plist[[1]])
legend_plot <- as_ggplot(p_legend)

pdf(here("plots", "04_label_transfer", "transferred_spds_on_tissue_legend.pdf"),
     width=8, height=4)
print(legend_plot)
dev.off()







