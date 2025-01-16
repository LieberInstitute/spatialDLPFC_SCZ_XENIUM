library(here)
library(SpatialExperiment)
library(scran)
library(tidyverse)
library(escheR)
library(scater)
library(scattermore)

spe <- readRDS(here("processed-data", "01_build_spe", "raw_spe_N24.RDS"))

plot_coldata_on_tissue <- function(x, column_name){
    plist <- list()
    brnums <- unique(x$BrNum)
    for (i in 1:length(brnums)){
        x_sub <- x[, x$BrNum == brnums[i]]
        p <- make_escheR(x_sub) %>%
            add_fill(column_name)+
            ggtitle(brnums[i])+
            geom_scattermore()
        plist[[i]] <- p
    }
    return(plist)
}



pdf(here("plots", "02_xenium_qc", "01_total_counts.pdf"))
plotColData(spe, y="total_counts", x="BrNum", colour_by="Dx")+
    geom_scattermore()
plotColData(spe, y="total_counts", x="Dx")+
    geom_scattermore()
total_counts_plots <- plot_coldata_on_tissue(spe, "total_counts")
lapply(total_counts_plots, print)
dev.off()


pdf(here("plots", "02_xenium_qc", "01_control_probe_counts.pdf"))
plotColData(spe, y="control_probe_counts", x="BrNum", colour_by="Dx")+
    geom_scattermore()
plotColData(spe, y="control_probe_counts", x="Dx")+
    geom_scattermore()
total_counts_plots <- plot_coldata_on_tissue(spe, "control_probe_counts")
lapply(total_counts_plots, print)
dev.off()

pdf(here("plots", "02_xenium_qc", "01_unassigned_codewords.pdf"))
plotColData(spe, y="unassigned_codeword_counts", x="BrNum", colour_by="Dx")+
    geom_scattermore()
plotColData(spe, y="unassigned_codeword_counts", x="Dx")+
    geom_scattermore()
total_counts_plots <- plot_coldata_on_tissue(spe, "unassigned_codeword_counts")
lapply(total_counts_plots, print)
dev.off()

pdf(here("plots", "02_xenium_qc", "01_cell_area.pdf"))
plotColData(spe, y="cell_area", x="BrNum", colour_by="Dx")+
    geom_scattermore()
plotColData(spe, y="cell_area", x="Dx")+
    geom_scattermore()
total_counts_plots <- plot_coldata_on_tissue(spe, "cell_area")
lapply(total_counts_plots, print)
dev.off()

pdf(here("plots", "02_xenium_qc", "01_nucleus_area.pdf"))
plotColData(spe, y="nucleus_area", x="BrNum", colour_by="Dx")+
    geom_scattermore()
plotColData(spe, y="nucleus_area", x="Dx")+
    geom_scattermore()
total_counts_plots <- plot_coldata_on_tissue(spe, "nucleus_area")
lapply(total_counts_plots, print)
dev.off()

pdf(here("plots", "02_xenium_qc", "01_transcript_counts.pdf"))
plotColData(spe, y="transcript_counts", x="BrNum", colour_by="Dx")+
    geom_scattermore()
plotColData(spe, y="transcript_counts", x="Dx")+
    geom_scattermore()
total_counts_plots <- plot_coldata_on_tissue(spe, "transcript_counts")
lapply(total_counts_plots, print)
dev.off()
