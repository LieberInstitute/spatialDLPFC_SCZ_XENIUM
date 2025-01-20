library(tidyverse)
library(escheR)
library(here)
library(SpatialExperiment)
library(scuttle)
library(SpotSweeper)


spe <- readRDS(here("processed-data", "01_build_spe", "raw_spe_N24.RDS"))

# compute number of genes detected
colData(spe)$detected <- Matrix::colSums(counts(spe) > 0)

is_neg <- stringr::str_detect(rownames(spe), "^NegControlProbe")
is_neg2 <- stringr::str_detect(rownames(spe), "^NegControlCodeword")
is_unassigned <- stringr::str_detect(rownames(spe), "^Unassigned")
is_mt <- stringr::str_detect(rownames(spe), "MT-")

# spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(negProbe = is_neg,
#                                                         negCodeword = is_neg2,
#                                                         unassigned = is_unassigned,
#                                                         mito=is_mt))

# Identify library size outliers
# spe <- localOutliers(spe,
#     metric = "total_counts",
#     direction = "lower",
#     log = TRUE,
#     n_neighbors=50
# )

brnums <- unique(spe$BrNum)
for(i in 1:length(brnums)){
    br_use <- brnums[[i]]
    spe_sub <- spe[, colData(spe)$BrNum == br_use]
    # Compute within-sample QC metrics based on negative controls and mitochondrial genes
    spe_sub <- scuttle::addPerCellQCMetrics(spe_sub, subsets = list(negProbe = is_neg,
                                                        negCodeword = is_neg2,
                                                        unassigned = is_unassigned,
                                                        mito=is_mt))
                                
    # Make some plots based on these metrics 
    p <- make_escheR(spe_sub)%>%
        add_fill("subsets_negProbe_percent")#+
        #ggtitle(paste(br_use, unique(spe_sub$Dx)[[1]]))
                                                        
}
