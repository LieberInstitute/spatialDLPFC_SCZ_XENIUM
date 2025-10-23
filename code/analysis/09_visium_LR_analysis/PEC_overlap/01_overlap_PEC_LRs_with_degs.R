library(here)
library(tidyverse)


microenvs <- c("neuropil", "neun", "pnn", "vasc")

for(microenv in microenvs){
    if(microenv=="vasc"){
        degs <- read.csv("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spg_pb_de/test_SPD_pseudo_vasc_pos.csv")
    }else{
        degs <- read.csv(here("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/",
            "code/analysis/",
            sprintf("dx_deg_spg_%s/", microenv), 
            sprintf("%s-dx_DEG-GM.csv", microenv)))
    }

    degs <- degs %>%
        filter(p_value_scz < 0.05)

    pec_lrs <- read.csv(here("processed-data", 
            "09_visium_LR_analysis", 
            "TableS5_SCZ_LorR_intersect_LIANA_op_annotation.csv"))

    pec_lrs <- pec_lrs %>%
        select(genesymbol_intercell_source, genesymbol_intercell_target, receptor_SCZ_risk, ligand_SCZ_risk) %>%
        filter(receptor_SCZ_risk==TRUE| ligand_SCZ_risk==TRUE)


    lrs_overlap <- pec_lrs %>%
        filter(genesymbol_intercell_source %in% degs$gene | genesymbol_intercell_target %in% degs$gene)

    write.csv(lrs_overlap, here("processed-data", 
            "09_visium_LR_analysis", 
            sprintf("PEC_LR_overlap_with_%s_DEGs.csv", microenv)), row.names=FALSE)
}