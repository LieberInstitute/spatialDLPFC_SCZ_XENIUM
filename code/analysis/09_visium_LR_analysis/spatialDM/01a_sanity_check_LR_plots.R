library(here)
library(SpatialExperiment)
library(ggplot2)
library(escheR)
library(patchwork)
library(tidyverse)
library(zellkonverter)

spe <- readRDS(
  here::here(
    "/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/",
    "processed-data/rds/01_build_spe",
    "fnl_spe_kept_spots_only.rds"
  )
)

spe$BDNF_counts <- assay(spe, "counts")[rowData(spe)$gene_name=="BDNF", ]
spe$NTRK2_counts <- assay(spe, "counts")[rowData(spe)$gene_name=="NTRK2", ]
spe$CADM1_counts <- assay(spe, "counts")[rowData(spe)$gene_name=="CADM1", ]

brnums <- unique(spe$BrNumbr)
pdf(here("plots", "09_visium_LR_analysis", "01a_sanity_check_BDNF_NTRK2_in_neuropil.pdf"), width=15, height=5)
for(i in 1:length(brnums)){
    br_use <- brnums[[i]]
    spe_sub <- spe[,spe$BrNumbr == br_use]
    #spe_sub <- spe_sub[,spe_sub$neuropil_pos==TRUE]

    p_bdnf <- make_escheR(spe_sub) %>%
        add_fill("BDNF_counts") %>%
        add_ground("neuropil_pos") +
        ggtitle(paste0("BDNF in neuropil - BrNum ", br_use))+
        scale_fill_gradient(low = "white", high = "black")+
        scale_color_manual(
        name = "", # turn off legend name for ground
        values = c("TRUE" = "red", "FALSE" = "transparent")
      )
        #scale_color_manual(c(TRUE = "red", FALSE = "grey"))
        #scale_color_manual(values=c(TRUE="red", FALSE="black"))

    p_ntrk2 <- make_escheR(spe_sub) %>%
        add_fill("NTRK2_counts")%>%
        add_ground("neuropil_pos")+
        ggtitle(paste0("NTRK2 in neuropil - BrNum ", br_use))+
        scale_fill_gradient(low = "white", high = "black")+
        scale_color_manual(
        name = "", # turn off legend name for ground
        values = c("TRUE" = "red", "FALSE" = "transparent")
      )
    p_cadm1 <- make_escheR(spe_sub) %>%
        add_fill("CADM1_counts")%>%
        add_ground("neuropil_pos")+
        ggtitle(paste0("CADM1 in neuropil - BrNum ", br_use))+
        scale_fill_gradient(low = "white", high = "black")+
        scale_color_manual(
        name = "", # turn off legend name for ground
        values = c("TRUE" = "red", "FALSE" = "transparent")
      )
    print(p_bdnf + p_ntrk2 + p_cadm1)

}
dev.off()
