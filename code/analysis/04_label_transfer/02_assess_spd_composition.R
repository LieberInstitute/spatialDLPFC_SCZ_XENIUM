library(here)
library(tidyverse)
library(SpatialExperiment)


####################################################################################################
# Read in the label transfer results, smooth the labels, and make pseudobulk PCA plots to
# see what's driving the variation in the pseudobulked samples
####################################################################################################

#res <- readRDS(here("processed-data", "04_label_transfer", "label_transfer_N24_k50.rds"))

# pdf(here("plots", "04_label_transfer", "label_transfer_faceted_plots.pdf"), width = 10, height = 10)
# for (i in 1:length(res$targets)){
#   spe <- res$targets[[i]]
#   spe <- spe[which(rowData(spe)$Type=="Gene Expression"), ]
#   spe$predictions <- as.character(spe$nmf_preds)
#   #print(unlist(strsplit(spe$predictions, split="0")))
#   spe$predictions_num <- unlist(lapply(strsplit(spe$predictions, split="0"), "[", 2))
#   print(head(spe$predictions_num))
#   spe$predictions_num <- as.numeric(spe$predictions_num)

#   spe$predictions_smooth <- nmfLabelTransfer::smoother(labels_curr=colData(spe)[,"predictions_num"],
#                                                           locs=spatialCoords(spe),
#                                                           props=NULL,
#                                                           k=10)
#   spe$predictions_smooth <- paste0("spd0", as.character(spe$predictions_smooth))
#   res$targets[[i]] <- spe

#   if (unique(colData(spe)$BrNum) %in% c("Br8667", "Br5973")){
#         plist_faceted <- nmfLabelTransfer::plot_faceted_clusters(spe, "predictions_smooth")
#         for (j in 1:length(plist_faceted)){
#           print(plist_faceted[[j]])
#         }
        
#       }
#   rm(spe)
#   gc()
# }
# dev.off()
# spe_all <- do.call(cbind, res$targets)
# rm(res)
# gc()
# print(spe_all)

# labs <- as.data.frame(colData(spe_all)[,"predictions_smooth"])
# rownames(labs) <- colnames(spe_all)
# print(head(labs))
# write.csv(labs,
#   here("processed-data", "04_label_transfer", "label_transfer_N24_k50_smoothed_labels.csv"))


# spe_all$sample_id <- paste(spe_all$BrNum, spe_all$Dx, sep="_")
# domain_annotations <- colData(spe_all) %>%
#   as.data.frame() %>%
#   mutate(domain_annotations = case_when(predictions_smooth == "spd07" ~ "L1",
#                                         predictions_smooth == "spd06" ~ "L2/3",
#                                         predictions_smooth == "spd02" ~ "L3/4",
#                                         predictions_smooth == "spd05" ~ "L5",
#                                         predictions_smooth == "spd03" ~ "L6",
#                                         predictions_smooth == "spd01" ~ "WMtz",
#                                         predictions_smooth == "spd04" ~ "WM",
#                                           TRUE ~ NA))


spe_all <- readRDS(here("processed-data", "01_build_spe", "raw_spe_N24.RDS"))
outlier_ids <- read.csv(here("processed-data", "02_xenium_qc", "outlier_ids.csv"))$x
spe_all <- spe_all[, -which(colnames(spe_all) %in% outlier_ids)]
rm(outlier_ids)
spds <- read.csv(here("processed-data", "04_label_transfer", "label_transfer_N24_k50_smoothed_labels.csv"))

rownames(spds) <- spds$X 
spds$X <- NULL
colnames(spds) <- c("predictions_smooth")
stopifnot(all(colnames(spe_all)==rownames(spds)))
colData(spe_all)$predictions_smooth <- factor(spds$predictions_smooth)

# map spd style labels to actual annotations
domain_annotations <- colData(spe_all) %>%
  as.data.frame() %>%
  mutate(domain_annotations = case_when(predictions_smooth == "spd07" ~ "L1/M",
                                        predictions_smooth == "spd06" ~ "L2/3",
                                        predictions_smooth == "spd02" ~ "L3/4",
                                        predictions_smooth == "spd05" ~ "L5",
                                        predictions_smooth == "spd03" ~ "L6",
                                        predictions_smooth == "spd01" ~ "WMtz",
                                        predictions_smooth == "spd04" ~ "WM",
                                          TRUE ~ NA))

colData(spe_all)$domain_annotations <- domain_annotations$domain_annotations
print(head(spe_all$domain_annotations))
colData(spe_all)$domain_annotations <- domain_annotations$domain_annotations

domain_levels <- c("L1/M", "L2/3", "L3/4", "L5", "L6", "WMtz", "WM")
colData(spe_all)$domain_annotations <- factor(colData(spe_all)$domain_annotations, levels = domain_levels)



spe_pseudo <- scuttle::aggregateAcrossCells(spe_all, 
      ids=DataFrame(domain_annotations=spe_all$domain_annotations, Dx=spe_all$sample_id))

spe_pseudo <- scuttle::logNormCounts(spe_pseudo, size.factors = NULL)
set.seed(1515)
spe_pseudo <- scater::runPCA(spe_pseudo, ncomponents=10)


domain_colors <- c(
    "L1/M" = "#FEAF16",
    "L2/3" = "#3283FE",
    "L3/4" = "#E4E1E3",
    "L5"   = "#16FF32",
    "L6"   = "#F6222E",
    "WMtz" = "#5A5156",
    "WM"   = "#FE00FA"
)
pdf(here("plots", "04_label_transfer", "pseudobulk_label_transfer_N24_k50.pdf"), width = 10, height = 10)
scater::plotPCA(spe_pseudo, colour_by = "Dx", ncomponents=4)
scater::plotPCA(spe_pseudo, colour_by = "domain_annotations", ncomponents=4, shape_by="Dx")+
           scale_colour_manual(
              name = "Spatial Domain",
              values = domain_colors,
              guide = guide_legend(override.aes = list(size = 7))
            ) 
dev.off()
 




