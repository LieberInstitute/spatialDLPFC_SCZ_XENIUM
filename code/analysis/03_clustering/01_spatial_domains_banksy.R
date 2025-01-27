library(Banksy)
library(SummarizedExperiment)
library(SpatialExperiment)
library(Seurat)
library(scater)
library(cowplot)
library(ggplot2)
library(here)


##########################################################################################
# Use BANKSY to do integration-free spatial domain analysis as a first pass
##########################################################################################

# Read in the data
spe <- readRDS(here("processed-data", "01_build_spe", "raw_spe_N24.RDS"))
outlier_ids <- read.csv(here("processed-data", "02_xenium_qc", "outlier_ids.csv"))$x
spe <- spe[, -which(colnames(spe) %in% outlier_ids)]
rm (outlier_ids)

spe <- scuttle::logNormCounts(spe)

# Create a list of each of the SPEs
spe_list <- list()
brnums <- unique(spe$BrNum)
for (i in 1:length(brnums)) {
  spe_list[[i]] <- spe[, spe$BrNum == brnums[i]]
}
compute_agf <- FALSE
k_geom <- 10
spe_list <- lapply(spe_list, computeBanksy, assay_name = "logcounts", 
                   compute_agf = compute_agf, k_geom = k_geom)

# Merge the SPEs back 
spe_joint <- do.call(cbind, spe_list)
rm(spe_list)
gc()

# Run BANKSY PCA
lambda <- 0.8
use_agf <- FALSE
spe_joint <- runBanksyPCA(spe_joint, use_agf = use_agf, 
                    lambda = lambda, group = "BrNum", seed = 1000)

# Run UMAP on the BANKSY PCA embedding
spe_joint <- runBanksyUMAP(spe_joint, use_agf = use_agf,  
            lambda = lambda, seed = 1000)

# BANKSY clustering
res <- 0.7
spe_joint <- clusterBanksy(spe_joint, use_agf = use_agf, lambda = lambda, resolution = res, seed = 1000)
cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(use_agf), lambda, res)


# split the spe again and plot each one 
pdf(here(sprintf("plots", "03_clustering", "banksy_clustering_lambda%s_res%s.pdf", lambda, res)))
for (i in 1:length(brnums)){
    sub_spe <- spe[, spe$BrNum == brnums[i]]

    make_escheR(sub_spe) %>%
        add_fill(cnm)+
        ggtitle(paste(br_use, unique(spe_sub$Dx)[[1]]))
    
}
dev.off()