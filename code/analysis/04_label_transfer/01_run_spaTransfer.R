library(nmfLabelTransfer)
library(SpatialExperiment)
library(SingleCellExperiment)
library(here)
library(tidyverse)
library(escheR)

# Read in the source dataset and merge in the PRECAST clusters based on code from Boyi
## Load Spe ----
raw_spe <- readRDS(
  here("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100",
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "spe_wo_spg_N63_PRECAST.rds"
  )
)

# error prevention
stopifnot(all(raw_spe$in_tissue))

# create sample_label
raw_spe$sample_label <- paste0(
  raw_spe$brnum, "_", toupper(raw_spe$dx)
)

## Load PRECAST df ----
PRECAST_df <- readRDS(
  here("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100",
    "processed-data/rds/spatial_cluster",
    "PRECAST",
    "test_clus_label_df_semi_inform_k_2-16.rds"
  )
)

## Merge PRECAST df ----
precast_vars <- grep(
  "^PRECAST_", colnames(PRECAST_df),
  value = TRUE
)

source_spe <- raw_spe[, raw_spe$key %in% PRECAST_df$key]
# raw_spe[, precast_vars] <- PRECAST_df[raw_spe$key, precast_vars]
col_data_df <- PRECAST_df |>
  right_join(
    colData(source_spe) |> data.frame(),
    by = c("key"),
    relationship = "one-to-one"
  )
rownames(col_data_df) <- colnames(source_spe)
colData(source_spe) <- DataFrame(col_data_df)

# error prevention
stopifnot(is.character(source_spe$PRECAST_07))
rm(raw_spe)
rm(PRECAST_df)
gc()

target_spe <- readRDS(here("processed-data", "01_build_spe", "raw_spe_N24.RDS"))
outlier_ids <- read.csv(here("processed-data", "02_xenium_qc", "outlier_ids.csv"))$x
target_spe <- target_spe[, -which(colnames(target_spe) %in% outlier_ids)]
rm (outlier_ids)
target_spe <- target_spe[which(rowData(target_spe)$Type=="Gene Expression"), ]
print(target_spe)

rowData(target_spe)$gene_name <- rownames(target_spe)

target_spe_list <- list()
brnums <- unique(target_spe$BrNum)
for (i in 1:length(brnums)) {
  target_spe_list[[i]] <- target_spe[, target_spe$BrNum == brnums[i]]
}
target_spe_list <- lapply(target_spe_list, scuttle::logNormCounts)

res <- transfer_labels(source=source_spe, target=target_spe_list, 
                    assay="logcounts",
                    annotationsName="PRECAST_07",
                    seed=0,
                    k=50,
                    technicalVarName="sample_label")

print(res$targets)
saveRDS(res, here("processed-data", "04_label_transfer", "label_transfer_N24_k50.rds"))


pdf(here("plots", "04_label_transfer", "label_transfer_N24_k50.pdf"), height=15, width=15)
for (i in 1:length(res$targets)){
    spe <- res$targets[[i]]
    p <- make_escheR(spe, y_reverse=FALSE)%>%
        add_fill("nmf_preds")+
        scale_fill_discrete()
    print(p)
}
dev.off()