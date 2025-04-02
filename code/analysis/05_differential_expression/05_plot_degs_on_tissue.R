library(tidyverse)
library(here)
library(ggplot2)
library(SpatialExperiment)
library(escheR)

xen_de <- read.csv(here("processed-data", "05_differential_expression", 
    "donor_domain_level_pseudobulk_Dx_DEGs_spaTransfer_k50_smoothed_predictions.csv"))


xen_de_sig <- xen_de %>%
  filter(fdr_SCZ <= 0.1)

degs <- xen_de_sig$X

spe <- readRDS(here("processed-data", "01_build_spe", "raw_spe_N24.RDS"))
spds <- read.csv(here("processed-data", "04_label_transfer", "label_transfer_N24_k50_smoothed_labels.csv"))

rownames(spds) <- spds$X 
spds$X <- NULL
colnames(spds) <- c("predictions_smooth")
spe <- spe[,colnames(spe) %in% rownames(spds)]
colData(spe) <- merge(colData(spe), spds, by="row.names", all.x=TRUE)

# remove non-gene expression counts
spe <- spe[rowData(spe)$Type == "Gene Expression",]


# Plot the DEGs on the tissue


brnums <- unique(spe$BrNum)
for (i in 1:length(brnums)){
    br_use <- brnums[[i]]
    pdf(here("plots", "05_differential_expression", "donor_domain_level_degs", sprintf("%s_donor_domain_level_degs_on_tissue.pdf", br_use)), height=5, width=15)
    sub_spe <- spe[, spe$BrNum == brnums[i]]
    sub_spe <- scuttle::logNormCounts(sub_spe)
    
    p_list <- list()
    for(j in 1:length(degs)){
        deg_use <- degs[[j]]
        sub_spe$counts_DEG <- counts(sub_spe)[which(rownames(sub_spe)== deg_use),]  
        sub_spe$logcounts_DEG <- logcounts(sub_spe)[which(rownames(sub_spe)== deg_use),]  

        p_list[[j]] <- make_escheR(sub_spe) %>%
          add_fill("counts_DEG")+
          ggtitle(paste(brnums[[i]], unique(sub_spe$Dx)[[1]], deg_use))+
          scattermore::geom_scattermore() 
        p <- make_escheR(sub_spe) %>%
          add_fill("logcounts_DEG")+
          ggtitle(paste(brnums[[i]], unique(sub_spe$Dx)[[1]], deg_use))+
          scattermore::geom_scattermore() 
        gridExtra::grid.arrange(p_list[[j]], p, ncol=2)
    }
    
    #do.call(gridExtra::grid.arrange, c(p_list, ncol=4))

    dev.off()

    }



