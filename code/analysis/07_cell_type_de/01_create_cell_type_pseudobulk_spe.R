suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
  library(escheR)
})

####################################################################################################
# Pseudobulk the SPE to the donor-domain level using the smoothed labels from label transfer.
# Reference: https://github.com/LieberInstitute/spatialDLPFC_SCZ/blob/main/code/analysis/pseudobulk_dx/create_pb_data.R
####################################################################################################

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


# Add some slide and date related metadata
slide_id <- unlist(lapply(strsplit(spe$Sample, split="/"), "[", 9))
slide_id <- unlist(lapply(strsplit(slide_id, split="__Br"), "[", 1))
spe$slide_id <- slide_id

run_date <- unlist(lapply(strsplit(spe$Sample, split="/"), "[", 7))
spe$run_date <- run_date


## First, create pseudobulk data:

spe_pseudo <- registration_pseudobulk(
      spe,
      var_registration = "annots",
      var_sample_id = "BrNum",
      covars = c("Dx", "Age", "Sex", "slide_id", "run_date"),
      min_ncells = 10,
      pseudobulk_rds_file = here(
        "processed-data", "07_cell_type_de",
        paste0(sprintf("spe_pseudo_donor_celltype_%s", cnm), "", ".rds")
      )
    )

# sanity check plots of spds
pdf(here("plots", "07_cell_type_de", "pseudobulk_sanity_check_plots.pdf"), width=10, height=10)
brnums <- unique(spe$BrNum)

for(i in 1:length(brnums)){
  spe_use <- spe[,spe$BrNum == brnums[i]]
  p <- make_escheR(spe_use) %>%
    add_fill("annots")+
    ggtitle(unique(spe_use$BrNum))+
    scale_fill_discrete()
  print(p)
}

dev.off()
