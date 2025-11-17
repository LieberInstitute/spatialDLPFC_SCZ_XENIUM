library(here)
library(tidyverse)
library(SpatialExperiment)
library(escheR)
library(patchwork)
library(ggsankey)




microenvs <- c("neuropil", "neun", "pnn", "vasc")

pdf(here("plots", 
        "09_visium_LR_analysis", "PEC_overlap",
        "PEC_LR_overlap_expression_patterns_sankey.pdf"
        ),height=8, width=12)
for(microenv in microenvs){


    #### read in LR pairs ####
    lr_pairs <- read.csv(here("processed-data", 
            "09_visium_LR_analysis", 
            sprintf("PEC_LR_overlap_with_%s_DEGs_SHK.csv", microenv)))
    lr_pairs <- lr_pairs %>%
        dplyr::rename_with(~"SHK_Tier", starts_with("SHK"))%>%
        dplyr::select(genesymbol_intercell_source, genesymbol_intercell_target)

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
    

    # lr_pairs <- lr_pairs %>%
    #     dplyr::mutate(logFC_SCZ = ifelse(lr_pairs$genesymbol_intercell_source %in% degs$gene,
    #                                 degs$logFC_scz[match(genesymbol_intercell_source, degs$gene)],
    #                                 degs$logFC_scz[match(genesymbol_intercell_target, degs$gene)]
    #                                 )) %>%
    #     dplyr::mutate(which_DE = ifelse(lr_pairs$genesymbol_intercell_source %in% degs$gene & lr_pairs$genesymbol_intercell_target %in% degs$gene,
    #                             "Both",
    #                             ifelse(lr_pairs$genesymbol_intercell_source %in% degs$gene,
    #                                     "Source",
    #                                     "Target"
    #                                     )
    #                             )
    #             ) 


    lr_pairs <- make_long(lr_pairs, genesymbol_intercell_source, genesymbol_intercell_target)
    p <- ggplot(lr_pairs, aes(x=x, next_x=next_x, node=node, next_node=next_node, label=node))+
        geom_alluvial(flow.alpha=0.5, node.color="lightgrey")+
        geom_alluvial_label(size=3,color = "white", fill = "gray40")+
        theme_alluvial(base_size=16)+
        ggtitle(sprintf("PEC LR pairs overlapping with %s DEGs", microenv))+
        theme(plot.title = element_text(hjust = 0.5))
    print(p)
}
dev.off()
