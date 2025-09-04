# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(spatialLIBD)
  # library(limma)
  library(sessioninfo)
  library(here)
  library(SpatialExperiment)
  library(ggplot2)
})

spe_pseudo <- readRDS(here("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/", 
                "processed-data/rds/PB_dx_spg",
                "pseudo_neuropil_pos_donor_spd.rds"))


gene_plot <- "BDNF"

counts_plt <- logcounts(spe_pseudo)[rowData(spe_pseudo)$gene_name==gene_plot,]
spe_pseudo$counts_plt <- counts_plt 

df <- colData(spe_pseudo) %>%
    as.data.frame() %>%
    select(counts_plt, spd_label, DX, sex, age, slide_id)


pdf(here("plots", "00_misc", "check_visium_boxplots", 
         "neuropil_boxplot_BDNF.pdf"), width = 15, height = 10)

ggplot(df, aes(x=factor(DX), y=counts_plt))+
    geom_violin()+
    geom_point(aes(colour=spd_label), position = position_jitter(seed = 1, width = 0.2))+
    ggtitle(gene_plot)


ggplot(df, aes(x=factor(DX), y=counts_plt))+
    geom_violin()+
    geom_point(aes(colour=spd_label), position = position_jitter(seed = 1, width = 0.2))+
    ggtitle(gene_plot)+
    facet_wrap(~spd_label)

ggplot(df, aes(x=factor(DX), y=counts_plt))+
    geom_violin()+
    geom_point(aes(colour=spd_label), position = position_jitter(seed = 1, width = 0.2))+
    ggtitle(gene_plot)+
    facet_wrap(~sex)
dev.off()





