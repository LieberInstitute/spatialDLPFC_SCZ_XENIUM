suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
})

####################################################################################################
# Create a QCed and normalized SPE object with cell type and spatial domain info in the colData
####################################################################################################

spe <- readRDS(here("processed-data", "01_build_spe", "raw_spe_N24.RDS"))
outlier_ids <- read.csv(here("processed-data", "02_xenium_qc", "outlier_ids.csv"))$x
spe <- spe[, -which(colnames(spe) %in% outlier_ids)]
rm(outlier_ids)
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
                       TRUE ~ "NA"))

cell_types_name <- sprintf("Banksy-%s-cell-types", cnm)
colData(spe)[[cell_types_name]] <- clusts$annots


# Add some slide and date related metadata
slide_id <- unlist(lapply(strsplit(spe$Sample, split="/"), "[", 9))
slide_id <- unlist(lapply(strsplit(slide_id, split="__Br"), "[", 1))
spe$slide_id <- slide_id

run_date <- unlist(lapply(strsplit(spe$Sample, split="/"), "[", 7))
spe$run_date <- run_date

### Spatial domain information
spds <- read.csv(here("processed-data", "04_label_transfer", "label_transfer_N24_k50_smoothed_labels.csv"))
print(head(spds))
rownames(spds) <- spds$X
spds$X <- NULL

colnames(spds) <- c("spatransfer_k50_predictions_smooth")
print(head(spds))

print(dim(spe))
print(dim(spds))

print(head(colnames(spe)))
print(head(rownames(spds)))
stopifnot(all(colnames(spe)==rownames(spds)))

colData(spe)$spatransfer_k50_predictions_smooth <- factor(spds$spatransfer_k50_predictions_smooth)
print(spe)

# map spd style labels to actual annotations
domain_annotations <- colData(spe) %>%
  as.data.frame() %>%
  mutate(domain_annotations = case_when(spatransfer_k50_predictions_smooth == "spd07" ~ "L1",
                                          spatransfer_k50_predictions_smooth == "spd06" ~ "L2/3",
                                          spatransfer_k50_predictions_smooth == "spd02" ~ "L3/4",
                                          spatransfer_k50_predictions_smooth == "spd05" ~ "L5",
                                          spatransfer_k50_predictions_smooth == "spd03" ~ "L6",
                                          spatransfer_k50_predictions_smooth == "spd01" ~ "WMtz",
                                          spatransfer_k50_predictions_smooth == "spd04" ~ "WM",
                                          TRUE ~ NA))



colData(spe)$domain_annotations <- domain_annotations$domain_annotations
print(spe$domain_annotations)

# Also add Sang Ho's gene annotations to the rowData for ease of access later
panel_markers <- readxl::read_xlsx((here("raw-data", 
        "experiment_info", 
        "Xenium_SHK_celltype_REannot_2025-04-13.xlsx")), sheet=2)
panel_markers$`...1` <- NULL

panel_markers <- panel_markers %>%
    as.data.frame() %>%
    mutate(cell_type_updated=case_when(cell_type_updated=="NA" ~ NA,
                                        TRUE ~ cell_type_updated))%>%
    mutate(Symbol=Gene)

df <- merge(rowData(spe), panel_markers, by="Symbol")
rowData(spe) <- df
                        
saveRDS(spe, here("processed-data", "07_cell_type_de", "cleaned_spe_N24_with_cell_type_and_spds.RDS"))