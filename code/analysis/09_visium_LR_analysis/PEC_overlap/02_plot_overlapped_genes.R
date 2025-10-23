library(here)
library(tidyverse)
library(SpatialExperiment)
library(escheR)
library(patchwork)

spe <- readRDS(
  here::here(
    "/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/",
    "processed-data/rds/01_build_spe",
    "fnl_spe_kept_spots_only.rds"
  )
)

microenvs <- c("neuropil", "neun", "pnn", "vasc")
brs_use <- c("Br8667", "Br5373")

spe_ntc <- spe[,spe$BrNumbr == "Br8667"]
spe_scz <- spe[,spe$BrNumbr == "Br5973"]



for(microenv in microenvs){
    lr_pairs <- read.csv(here("processed-data", 
            "09_visium_LR_analysis", 
            sprintf("PEC_LR_overlap_with_%s_DEGs_SHK.csv", microenv)))
    lr_pairs <- lr_pairs %>%
        rename_with(~"SHK_Tier", starts_with("SHK"))%>%
        select(genesymbol_intercell_source, genesymbol_intercell_target)
    
    pdf(here("plots", 
            "09_visium_LR_analysis", 
            sprintf("PEC_LR_overlap_with_%s_DEGs_expression_patterns.pdf", microenv)
            ),height=6, width=15)
    
    for(i in 1:nrow(lr_pairs)){
        receptor <- lr_pairs$genesymbol_intercell_target[i]
        ligand <- lr_pairs$genesymbol_intercell_source[i]

        spe_ntc$receptor_counts <- counts(spe_ntc)[rowData(spe_ntc)$gene_name == receptor]
        spe_ntc$ligand_counts <- counts(spe_ntc)[rowData(spe_ntc)$gene_name == ligand]

        spe_scz$receptor_counts <- counts(spe_scz)[rowData(spe_scz)$gene_name == receptor]
        spe_scz$ligand_counts <- counts(spe_scz)[rowData(spe_scz)$gene_name == ligand]

        # determine whether each spot has expression of receptor, ligand, or both
        spe_ntc$expression_pattern <- case_when(
            spe_ntc$receptor_counts > 0 & spe_ntc$ligand_counts > 0 ~ "Both",
            spe_ntc$receptor_counts > 0 & spe_ntc$ligand_counts == 0 ~ "Receptor only",
            spe_ntc$receptor_counts == 0 & spe_ntc$ligand_counts > 0 ~ "Ligand only",
            TRUE ~ "Neither"
        )

        spe_scz$expression_pattern <- case_when(
            spe_scz$receptor_counts > 0 & spe_scz$ligand_counts > 0 ~ "Both",
            spe_scz$receptor_counts > 0 & spe_scz$ligand_counts == 0 ~ "Receptor only",
            spe_scz$receptor_counts == 0 & spe_scz$ligand_counts > 0 ~ "Ligand only",
            TRUE ~ "Neither"
        )

        p1 <- make_escheR(spe_ntc) %>%
            add_fill("expression_pattern") %>%
            add_ground(sprintf("%s_pos", microenv))+
            ggtitle(sprintf("%s: %s - %s", unique(spe_ntc$BrNumbr), ligand, receptor))+
            scale_fill_manual(values=c("Both"="purple", "Receptor only"="lightblue", 
                        "Ligand only"="lightpink", "Neither"="lightgrey"))+
            scale_color_manual(
                    name = "", # turn off legend name for ground
                    values = c("TRUE" = "red", "FALSE" = "transparent")
                )

        p2 <- make_escheR(spe_scz) %>%
            add_fill("expression_pattern")%>%
            add_ground(sprintf("%s_pos", microenv))+
            ggtitle(sprintf("%s: %s - %s", unique(spe_scz$BrNumbr), ligand, receptor))+
            scale_fill_manual(values=c("Both"="purple", "Receptor only"="lightblue", 
                        "Ligand only"="lightpink", "Neither"="lightgrey"))+
            scale_color_manual(
                    name = "", # turn off legend name for ground
                    values = c("TRUE" = "red", "FALSE" = "transparent")
                )

        print(p1 + p2)
    }
    dev.off()
    

}
