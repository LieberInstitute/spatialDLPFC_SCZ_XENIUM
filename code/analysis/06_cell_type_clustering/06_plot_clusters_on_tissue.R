suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
  library(escheR)
  library(patchwork)
})

spe <- readRDS(here("processed-data", "01_build_spe", "raw_spe_N24.RDS"))
outlier_ids <- read.csv(here("processed-data", "02_xenium_qc", "outlier_ids.csv"))$x
spe <- spe[, -which(colnames(spe) %in% outlier_ids)]
rm (outlier_ids)
spe <- spe[which(rowData(spe)$Type=="Gene Expression"), ]
print(spe)

#### Banksy parameters ####
lambda <- 0.1
res <- 0.7 # higher = more clusters
use_agf <- FALSE
cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(use_agf), lambda, res)



clusts <- read.csv(here("processed-data", "06_cell_type_clustering", sprintf("banksy_clustering_lambda%s_res%s.csv", lambda, res)))
clusts <- as.data.frame(clusts)
print(head(clusts))
colData(spe)[["Banksy"]] <- as.character(clusts$V1)


# Aggregate to broader groups

clusts <- clusts %>% as.data.frame() %>% 
  mutate(V1 = as.numeric(V1)) %>%
  mutate(annots=case_when(V1 %in% c(6, 8, 1) ~ "Oligo", # 10 and 15 i am not super sure about
                          V1 %in% c(10) ~ "Ambig/Oligo",
                          V1 %in% c(7) ~ "Mic",
                          V1 %in% c(15) ~ "Ambig/In/Endo",
                          V1 %in% c(18) ~ "L5 Ex",
                          V1 %in% c(14) ~ "L6 Ex",
                          V1 %in% c(11) ~ "L4/5 Ex",
                          V1 %in% c(2) ~ "L2/3 Ex",
                          V1 %in% c(17, 16, 4) ~ "Ast",
                          V1 %in% c(5, 13, 3) ~ "Endo",
                          V1 %in% c(12) ~ "In: VIP, LAMP5",
                          V1 %in% c(9) ~ "In: SST, PVALB",
                       TRUE ~ "NA")) %>%
    mutate(annots_combined = paste(annots, V1, sep="."))
  
colData(spe)[["annots"]] <- clusts$annots
colData(spe)[["annots_combined"]] <- clusts$annots


brs_use <- c("Br5973", "Br8667")
plist <- list()
for(i in 1:length(brs_use)){
    br_plt <- brs_use[i]
    spe_use <- spe[,spe$BrNum==br_plt]

    p <- make_escheR(spe_use, y_reverse=FALSE) %>%
        add_fill("annots")+ 
        ggtitle(paste(brs_use[i], unique(spe_use$Dx), sep=" "))+
        scale_fill_discrete()+
        labs(fill = "Cell Type")

    plist[[i]] <- p
}

png(here("plots", "06_cell_type_clustering", "cell_types_on_tissue.png"),
     width=3500, height=2000, res=300)
wrap_plots(plist, ncol=2)+
    plot_layout(guides="collect")
dev.off()