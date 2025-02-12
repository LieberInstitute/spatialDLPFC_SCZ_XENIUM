library(nmfLabelTransfer)
library(SpatialExperiment)
library(SingleCellExperiment)
library(here)
library(tidyverse)
library(escheR)
library(ggplot2)
library(scattermore)

# Read in the source dataset and merge in the PRECAST clusters based on code from Boyi
if (!file.exists(here("processed-data", "04_label_transfer", "label_transfer_N24_k50.rds"))){
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

}else{
  res <- readRDS(here("processed-data", "04_label_transfer", "label_transfer_N24_k50.rds"))
}
pdf(here("plots", "04_label_transfer", "label_transfer_N24_k50.pdf"), height=15, width=25)
  for (i in 1:length(res$targets)){
      spe <- res$targets[[i]]
      spe$predictions <- as.character(spe$nmf_preds)

      #print(unlist(strsplit(spe$predictions, split="0")))
      spe$predictions_num <- unlist(lapply(strsplit(spe$predictions, split="0"), "[", 2))
      print(head(spe$predictions_num))
      spe$predictions_num <- as.numeric(spe$predictions_num)


      #props <- c(0.5, 0.8, 0.6, 0.5, 0.2, 0.4, 0.5)
      #names(props) <- 1:7
      spe$predictions_smooth <- nmfLabelTransfer::smoother(labels_curr=colData(spe)[,"predictions_num"],
                                                        locs=spatialCoords(spe),
                                                        props=NULL,
                                                        k=10)
      spe$predictions_smooth <- paste0("spd0", as.character(spe$predictions_smooth))


      p <- make_escheR(spe, y_reverse=FALSE) %>%
          add_fill("predictions")+
           scale_fill_manual(
              name = "Spatial Domain",
              values = set_names(
                Polychrome::palette36.colors(7)[seq.int(7)],
                unique(spe$predictions) |> sort()
              ),
              guide = guide_legend(override.aes = list(size = 7))
            ) +
          ggtitle(paste(unique(colData(spe)$BrNum), unique(colData(spe)$Dx)))

      p_smooth <- make_escheR(spe, y_reverse=FALSE) %>%
          add_fill("predictions_smooth")+
           scale_fill_manual(
              name = "Smoothed Spatial Domain",
              values = set_names(
                Polychrome::palette36.colors(7)[seq.int(7)],
                unique(spe$predictions_smooth) |> sort()
              ),
              guide = guide_legend(override.aes = list(size = 7))
            ) +
          ggtitle(paste(unique(colData(spe)$BrNum), unique(colData(spe)$Dx)))

        gridExtra::grid.arrange(p, p_smooth, ncol=2)

      # Make faceted plots for Br8667 and Br5973:
      if (unique(colData(spe)$BrNum) %in% c("Br8667", "Br5973")){


        plist_faceted <- nmfLabelTransfer::plot_faceted_clusters(spe, "predictions_smooth")
        do.call(gridExtra::grid.arrange, c(plist_faceted, ncol=2))
        
      }
}
dev.off()

