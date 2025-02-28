library(tidyverse)
library(here)
library(ggplot2)


vis_de <- read.csv(here("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100", 
            "processed-data/PB_dx_genes",
            "test_PRECAST_07.csv"))

xen_de <- read.csv(here("processed-data", "05_differential_expression", 
    "donor_domain_level_pseudobulk_Dx_DEGs_spaTransfer_k50_smoothed_predictions.csv"))

probe_list <- readxl::read_excel(here("raw-data",
                "experiment_info", 
                "Xenium_SCZ_ProbeSelection5_SHK_v1.xlsx"),
                sheet="SHK_Final_300_v1", col_names=FALSE)

colnames(probe_list) <- c("Gene", "cell_type", "dx_deg", "layer_marker")

probe_list_updated <- readxl::read_excel(here("raw-data",
                "experiment_info", 
                "Xenium_SCZ_ProbeSelection5_SHK_v4_15_w_ensemblgene_list.xlsx"),
                col_names=TRUE)

probe_list_final <- merge(probe_list_updated, probe_list, by="Gene", all.x=TRUE) 

xen_de_sig <- xen_de %>%
  filter(fdr_SCZ <= 0.1)

probe_list[probe_list$gene %in% xen_de_sig$gene,]

vis_de <- vis_de %>% 
    mutate(logFC_SCZ = logFC_scz) %>%
    mutate(p_value_SCZ = p_value_scz) %>% 
    select(gene, logFC_SCZ, p_value_SCZ) %>%
    filter(gene %in% xen_de$gene)

xen_de <- xen_de %>%
    select(gene, logFC_SCZ, p_value_SCZ) %>%
    filter(gene %in% vis_de$gene)

all_de <- merge(vis_de, xen_de, by="gene", all=TRUE, suffixes=c("_vis", "_xen"))

pdf(here("plots", "05_differential_expression", "visium_xenium_de_logFC.pdf"))
ggplot(all_de, aes(x=logFC_SCZ_vis, y=logFC_SCZ_xen)) +
    geom_point() +
    geom_abline(intercept=0, slope=1, color="red") +
    labs(x="Visium logFC", y="Xenium logFC")
dev.off()

# Some sanity checks
sum(!(is.na(probe_list_final$dx_deg)))
probe_list_final[!is.na(probe_list_final$dx_deg),]

xen_de %>% 
    filter(fdr_SCZ <= 0.05) %>%
    nrow()


