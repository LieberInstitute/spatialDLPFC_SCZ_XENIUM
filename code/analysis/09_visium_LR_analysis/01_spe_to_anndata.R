library(here)
library(SpatialExperiment)
library(ggplot2)
library(escheR)
library(patchwork)
library(tidyverse)
#library(zellkonverter)

spe <- readRDS(
  here::here(
    "/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/",
    "processed-data/rds/01_build_spe",
    "fnl_spe_kept_spots_only.rds"
  )
)

dim(spe)
#[1]  36601 279806

spe_neurons <- spe[,spe$neun_pos == TRUE]
dim(spe_neurons)
# [1] 36601 55017

spe_neuropil <- spe[,spe$neuropil_pos == TRUE]
dim(spe_neuropil)
# [1]  36601 147965

spe_pnn <- spe[,spe$pnn_pos == TRUE]
dim(spe_pnn)
# [1] 36601 12333

spe_vasc <- spe[,spe$vasc_pos == TRUE]
dim(spe_vasc)
# [1] 36601 21906

# check overlaps between the spot types
sum(colnames(spe_neurons) %in% colnames(spe_neuropil))
# [1] 13717

sum(colnames(spe_neurons) %in% colnames(spe_pnn))
# [1] 4365
sum(colnames(spe_neurons) %in% colnames(spe_vasc))
# [1] 3478
sum(colnames(spe_neuropil) %in% colnames(spe_pnn))
# [1] 5382
sum(colnames(spe_neuropil) %in% colnames(spe_vasc))
# [1] 5598
sum(colnames(spe_pnn) %in% colnames(spe_vasc))
# [1] 765

# Make spot plots of each type of spe to see where the spots are located
brnums <- unique(spe$BrNumbr)

neur_plots <- list()
neuropil_plots <- list()
pnn_plots <- list()
vasc_plots <- list()
for(i in 1:length(brnums)){
    br_use <- brnums[[i]]
    spe_sub <- spe[,spe$BrNumbr == br_use]

    p_neurons <- make_escheR(spe_sub) %>%
        add_fill("neun_pos") +
        scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "lightgrey")) +
        ggtitle(paste0(br_use, " Neurons"))
    
    p_neuropil <- make_escheR(spe_sub) %>%
        add_fill("neuropil_pos") +
        scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "lightgrey")) +
        ggtitle(paste0(br_use, " Neuropil"))
    
    p_pnn <- make_escheR(spe_sub) %>%
        add_fill("pnn_pos") +
        scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "lightgrey")) +
        ggtitle(paste0(br_use, " PNN"))
    
    p_vasc <- make_escheR(spe_sub) %>%
        add_fill("vasc_pos") +
        scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "lightgrey")) +
        ggtitle(paste0(br_use, " Vasc"))

    neur_plots[[i]] <- p_neurons
    neuropil_plots[[i]] <- p_neuropil
    pnn_plots[[i]] <- p_pnn
    vasc_plots[[i]] <- p_vasc

}

pdf(here("plots", "09_visium_LR_analysis", "spot_types_by_brnum.pdf"), 
    width = 20, height = 15)
print(wrap_plots(neur_plots, ncol=8) + plot_annotation(title = "Neurons"))
print(wrap_plots(neuropil_plots, ncol=8) + plot_annotation(title = "Neuropil"))
print(wrap_plots(pnn_plots, ncol=8) + plot_annotation(title = "PNN"))
print(wrap_plots(vasc_plots, ncol=8) + plot_annotation(title = "Vasc"))
dev.off()

# Turn the SPE objects into SCEs so that zellkonverter can convert them to AnnData
colData(spe_neurons)$spot_type <- "neurons"
colData(spe_neurons)$coord_x <- spatialCoords(spe_neurons)[,1]
colData(spe_neurons)$coord_y <- spatialCoords(spe_neurons)[,2]
sce_neurons <- as(spe_neurons, "SingleCellExperiment")
print(sce_neurons)
writeH5AD(sce_neurons, 
          here("processed-data", "09_visium_LR_analysis", "sce_neurons.h5ad"))

colData(spe_neuropil)$spot_type <- "neuropil"
colData(spe_neuropil)$coord_x <- spatialCoords(spe_neuropil)[,1]
colData(spe_neuropil)$coord_y <- spatialCoords(spe_neuropil)[,2]
sce_neuropil <- as(spe_neuropil, "SingleCellExperiment")
print(sce_neuropil)
writeH5AD(sce_neuropil, 
          here("processed-data", "09_visium_LR_analysis", "sce_neuropil.h5ad"))

colData(spe_pnn)$spot_type <- "pnn"
colData(spe_pnn)$coord_x <- spatialCoords(spe_pnn)[,1]
colData(spe_pnn)$coord_y <- spatialCoords(spe_pnn)[,2]
sce_pnn <- as(spe_pnn, "SingleCellExperiment")
print(sce_pnn)
writeH5AD(sce_pnn, 
          here("processed-data", "09_visium_LR_analysis", "sce_pnn.h5ad"))

colData(spe_vasc)$spot_type <- "vasc"
colData(spe_vasc)$coord_x <- spatialCoords(spe_vasc)[,1]
colData(spe_vasc)$coord_y <- spatialCoords(spe_vasc)[,2]
sce_vasc <- as(spe_vasc, "SingleCellExperiment")
print(sce_vasc)
writeH5AD(sce_vasc, 
          here("processed-data", "09_visium_LR_analysis", "sce_vasc.h5ad"))



