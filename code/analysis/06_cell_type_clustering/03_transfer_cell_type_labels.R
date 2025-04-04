library(SingleCellExperiment)
library(SpatialExperiment)
library(tidyverse)
library(here)
library(nmfLabelTransfer)
library(scuttle)
library(escheR)




######################################################
# Transfer labels from DLPFC snRNA-seq data into 
# the Xenium datasets 
######################################################

sce_path <- "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata"
load(sce_path)
spe <- readRDS(here("processed-data", "01_build_spe", "raw_spe_N24.RDS"))
outlier_ids <- read.csv(here("processed-data", "02_xenium_qc", "outlier_ids.csv"))$x
spe <- spe[, -which(colnames(spe) %in% outlier_ids)]
rm (outlier_ids)
spe <- spe[which(rowData(spe)$Type=="Gene Expression"), ]
print(spe)


spes_list <- list()
brnums <- unique(spe$BrNum)
for(i in 1:length(brnums)){
    br_use <- brnums[[i]]
    spe_use <- logNormCounts(spe[,spe$BrNum == br_use])
    rowData(spe_use)$gene_name <- rowData(spe_use)$Symbol 
    spes_list[[i]] <- spe_use
    print(head(rowData(spe_use)))
    rm(spe_use)
    gc()
}

k <- 50
annots <- "cellType_broad_hc"

res_fname <- here("processed-data", "06_cell_type_clustering", sprintf("cell_type_label_transfer_N24_k%s_%s.rds", k, annots))
if(!file.exists(res_fname)){
  res <- transfer_labels(source=sce, target=spes_list, 
                      assay="logcounts",
                      annotationsName=annots,
                      seed=0,
                      k=k,
                      technicalVarName="Sample")
  saveRDS(res, res_fname)
}else{
  res <- readRDS(res_fname)
}
pdf(here("plots", "06_cell_type_clustering", sprintf("cell_type_label_transfer_N24_k%s_%s.pdf", k, annots)), height=15, width=25)
  for (i in 1:length(res$targets)){
      spe <- res$targets[[i]]
      spe$predictions <- as.character(spe$nmf_preds)


      p <- make_escheR(spe, y_reverse=FALSE) %>%
          add_fill("predictions")+
          scale_fill_discrete()+
          ggtitle(paste(unique(colData(spe)$BrNum), unique(colData(spe)$Dx)))


    print(p)
      # Make faceted plots for Br8667 and Br5973:
      if (unique(colData(spe)$BrNum) %in% c("Br8667", "Br5973")){
        plist_faceted <- nmfLabelTransfer::plot_faceted_clusters(spe, "predictions")
        do.call(gridExtra::grid.arrange, c(plist_faceted, ncol=2))
        
      }
}
dev.off()

