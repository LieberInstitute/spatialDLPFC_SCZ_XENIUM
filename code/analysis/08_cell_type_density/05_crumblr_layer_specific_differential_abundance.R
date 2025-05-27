suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(sessioninfo)
  library(escheR)
  library(crumblr)
  library(variancePartition)
  library(dreamlet)
})

####################################################################################################
# Perform layer-specific differential abundance analysis of cell types
# Reference: https://github.com/LieberInstitute/BLA_crossSpecies/blob/devel/code/07_annotation_and_characterization/14_crumblr_inhib_BLAsubregions.R
# https://diseaseneurogenomics.github.io/crumblr/articles/crumblr.html
####################################################################################################

## using conda_R/4.5
spe <- readRDS(here("processed-data", "07_cell_type_de", 
            "cleaned_spe_N24_with_cell_type_and_spds.RDS"))

spe$cell_types <- colData(spe)[["Banksy-clust_M0_lam0.1_k50_res0.7-cell-types"]]

brnums <- unique(spe$BrNum)

layers <- unique(spe$domain_annotations)
spe$Dx <- as.factor(spe$Dx)
spe$BrNum <- as.factor(spe$BrNum)
spe$domain_annotations <- as.factor(spe$domain_annotations)
spe$cell_types <- as.factor(spe$cell_types)

# create unique identifier based on Dx and BrNum
spe$id <- as.factor(paste(spe$Dx, spe$BrNum, spe$cell_types, spe$domain_annotations, sep = "_"))


#--------------------------
# Create pseudobulked data
#--------------------------
pb <- aggregateToPseudoBulk(spe,
    assay = "counts",
    cluster_id = "cell_types",
    sample_id = "id",
    verbose = TRUE
)
cobj <- crumblr(cellCounts(pb)) # crumblr object
# =========
# Dream DE
# =========
form = ~ + Dx + (1|BrNum) + domain_annotations
fit <- dream(cobj, form, colData(pb))
fit <- eBayes(fit)

# Extract results for each cell type
topTable(fit, coef = "DxSCZ", number = Inf)

vp <- fitExtractVarPartModel(cobj, form, colData(pb))
# layer adjusted analysis?
# pseudobulk to donor-domain-cell type level?



