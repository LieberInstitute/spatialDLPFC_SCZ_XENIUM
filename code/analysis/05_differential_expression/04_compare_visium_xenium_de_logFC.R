library(tidyverse)
library(here)
library(ggplot2)
library(ggrepel)


vis_de <- read.csv(here("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100", 
            "processed-data/rds/10_dx_deg_adjust_spd/dx-deg_PRECAST07.csv"))

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
write.csv(probe_list_final, here("raw-data",
                "experiment_info", 
                "combined_Xenium_SCZ_ProbeSelection5_SHK_v4_15_w_ensemblgene_list.xlsx"))

#xen_de_sig <- xen_de %>%
 # filter(fdr_SCZ <= 0.1)

#probe_list[probe_list$gene %in% xen_de_sig$gene,]
colnames(probe_list_final)[[ncol(probe_list_final)]] <- "notes"
dx_deg <- probe_list_final %>%
   as.data.frame() %>%
    filter(!is.na(dx_deg)) %>%
    select(Gene, dx_deg)

vis_de <- vis_de %>% 
    mutate(logFC_SCZ = logFC_scz) %>%
    mutate(fdr_SCZ = fdr_scz) %>% 
    select(gene, logFC_SCZ, fdr_SCZ) %>%
    filter(gene %in% xen_de$gene) %>%
    filter(gene %in% dx_deg$Gene)

xen_de <- xen_de %>%
    select(gene, logFC_SCZ, fdr_SCZ) %>%
    filter(gene %in% vis_de$gene)

all_de <- merge(vis_de, xen_de, by="gene", all=TRUE, suffixes=c("_vis", "_xen"))
all_de <- all_de %>%
  mutate(
    direction = case_when(
      logFC_SCZ_vis > 0 & logFC_SCZ_xen > 0 ~ "both upregulated",
      logFC_SCZ_vis < 0 & logFC_SCZ_xen  < 0 ~ "both downregulated",
      TRUE ~ "mixed"
    )
  )

# create a dataframe of the genes that are significant in xenium
xen_de_sig <- all_de %>%
    filter(fdr_SCZ_xen <= 0.1)
print(xen_de_sig)

pdf(here("plots", "05_differential_expression", "visium_xenium_de_logFC.pdf"))
ggplot(all_de, aes(x=logFC_SCZ_vis, y=logFC_SCZ_xen, colour=direction)) +
    geom_point() +
    geom_abline(intercept=0, slope=1, color="red") +
    scale_color_manual(values = c("both upregulated" = "red", "both downregulated" = "blue", "mixed" = "grey"))+
    labs(x="Visium logFC", y="Xenium logFC")+
    geom_label_repel(size=3,
        data = xen_de_sig,
        aes(label = gene),
        force = 4,
        nudge_y = 0.2
      )+
      theme_minimal() 
dev.off()

# Some sanity checks
sum(!(is.na(probe_list_final$dx_deg)))
probe_list_final[!is.na(probe_list_final$dx_deg),]

# xen_de %>% 
#     filter(fdr_SCZ <= 0.05) %>%
#     nrow()


