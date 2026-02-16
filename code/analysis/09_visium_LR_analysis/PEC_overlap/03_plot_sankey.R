library(here)
library(tidyverse)
library(SpatialExperiment)
library(escheR)
library(patchwork)
library(ggsankey)
library(patchwork)




microenvs <- c("neuropil", "pnn", "neun", "vasc")


p_list <- list()
titles <- list("Neuropil", "PNN", "Neuronal", "Vasculature")
colours <- list("grey", "magenta", "yellow", "green")
for(i in 1:length(microenvs)){
    #### read in LR pairs ####
    microenv <- microenvs[i]
    lr_pairs <- read.csv(here("processed-data", 
            "09_visium_LR_analysis", 
            sprintf("PEC_LR_overlap_with_%s_DEGs_SHK.csv", microenv)))
    lr_pairs <- lr_pairs %>%
        dplyr::rename_with(~"SHK_Tier", starts_with("SHK"))%>%
        dplyr::select(genesymbol_intercell_source, genesymbol_intercell_target) %>%
        dplyr::rename(Ligand=genesymbol_intercell_source,
                      Receptor=genesymbol_intercell_target)
    #### read in DE results #####
    if(microenv=="vasc"){
        degs <- read.csv("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spg_pb_de/test_SPD_pseudo_vasc_pos.csv")
    }else{
        degs <- read.csv(here("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/",
            "code/analysis/",
            sprintf("dx_deg_spg_%s/", microenv), 
            sprintf("%s-dx_DEG-GM.csv", microenv)))
    }

    degs <- degs %>%
        dplyr::filter(p_value_scz < 0.05)
    
    lr_pairs <- make_long(lr_pairs, Ligand, Receptor)
    p <- ggplot(lr_pairs, aes(x=x, next_x=next_x, node=node, next_node=next_node, label=node))+
        geom_alluvial(flow.alpha=0.5, node.color="lightgrey", fill=colours[[i]], color="grey50")+
        geom_alluvial_label(size=6, color = "white", fill = "grey40")+
        theme_alluvial(base_size=16)+
        ggtitle(titles[[i]])+
        theme(plot.title = element_text(hjust = 0.5, size=28),axis.text.y=element_blank(),
            axis.ticks=element_blank(), axis.title.x=element_blank(),
            axis.text.x=element_text(size=24))
        
    p_list[[microenv]] <- p
}
pdf(here("plots", 
        "09_visium_LR_analysis", "PEC_overlap",
        "PEC_LR_overlap_expression_patterns_sankey.pdf"),
        height=14, width=65.8)

wrap_plots(p_list, ncol=4) 
dev.off()
