library(here)
library(tidyverse)
library(SpatialExperiment)


####################################################################################################
# Read in the label transfer results, smooth the labels, and make pseudobulk PCA plots to
# see what's driving the variation in the pseudobulked samples
####################################################################################################

res <- readRDS(here("processed-data", "04_label_transfer", "label_transfer_N24_k50.rds"))

pdf(here("plots", "04_label_transfer", "label_transfer_faceted_plots.pdf"), width = 10, height = 10)
for (i in 1:length(res$targets)){
  spe <- res$targets[[i]]
  spe$predictions <- as.character(spe$nmf_preds)
  #print(unlist(strsplit(spe$predictions, split="0")))
  spe$predictions_num <- unlist(lapply(strsplit(spe$predictions, split="0"), "[", 2))
  print(head(spe$predictions_num))
  spe$predictions_num <- as.numeric(spe$predictions_num)

  spe$predictions_smooth <- nmfLabelTransfer::smoother(labels_curr=colData(spe)[,"predictions_num"],
                                                          locs=spatialCoords(spe),
                                                          props=NULL,
                                                          k=10)
  spe$predictions_smooth <- paste0("spd0", as.character(spe$predictions_smooth))
  res$targets[[i]] <- spe

  if (unique(colData(spe)$BrNum) %in% c("Br8667", "Br5973")){
        plist_faceted <- nmfLabelTransfer::plot_faceted_clusters(spe, "predictions_smooth")
        for (j in 1:length(plist_faceted)){
          print(plist_faceted[[j]])
        }
        
      }
  rm(spe)
  gc()
}
dev.off()
spe_all <- do.call(cbind, res$targets)
rm(res)
gc()
print(spe_all)

write.csv(here("processed_data", "04_label_transfer", "label_transfer_N24_k50_smoothed_labels.csv"), colData(spe_all)[,"predictions_smooth"])

spe_all$sample_id <- paste(spe_all$BrNum, spe_all$Dx, sep="_")
spe_pseudo <- scuttle::aggregateAcrossCells(spe_all, 
      ids=DataFrame(predictions_smooth=spe_all$predictions_smooth, Dx=spe_all$sample_id))

spe_pseudo <- scuttle::logNormCounts(spe_pseudo, size.factors = NULL)
set.seed(1515)
spe_pseudo <- scater::runPCA(spe_pseudo, ncomponents=10)

pdf(here("plots", "04_label_transfer", "pseudobulk_label_transfer_N24_k50.pdf"), width = 10, height = 10)
scater::plotPCA(spe_pseudo, colour_by = "Dx", ncomponents=4)
scater::plotPCA(spe_pseudo, colour_by = "predictions_smooth", ncomponents=4)+
           scale_colour_manual(
              name = "Spatial Domain",
              values = set_names(
                Polychrome::palette36.colors(7)[seq.int(7)],
                unique(spe_pseudo$predictions_smooth) |> sort()
              ),
              guide = guide_legend(override.aes = list(size = 7))
            ) 
dev.off()
 


