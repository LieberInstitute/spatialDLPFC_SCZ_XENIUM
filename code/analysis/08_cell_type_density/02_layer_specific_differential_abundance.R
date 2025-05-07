suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
  library(escheR)
  library(readxl)
  library(edgeR)

})

####################################################################################################
# Perform layer-specific differential abundance analysis of cell types
# Reference: https://bioconductor.org/books/3.13/OSCA.multisample/differential-abundance.html
####################################################################################################


spe <- readRDS(here("processed-data", "07_cell_type_de", 
            "cleaned_spe_N24_with_cell_type_and_spds.RDS"))

spe$cell_types <- colData(spe)[["Banksy-clust_M0_lam0.1_k50_res0.7-cell-types"]]

brnums <- unique(spe$BrNum)

layers <- unique(spe$spatransfer_k50_predictions_smooth)

# Loop through each layer and perform differential cell type abundance analysis
pdf(here("plots", "08_cell_type_density", "layer_specific_differential_abundance.pdf"))
for (i in 1:length(layers)){
    layer_use <- layers[i]
    print(layer_use)
    spe_sub <- spe[,spe$spatransfer_k50_predictions_smooth == layer_use]

    abundances <- table(spe_sub$cell_types, spe_sub$BrNum)
    print(abundances)
    # Attach metadata and create DGEList object
    extra.info <- colData(spe_sub)[match(colnames(abundances), spe_sub$BrNum),]
    y.ab <- DGEList(abundances, samples=extra.info)
    y.ab

    # Filter out low abundance cell types
    # keep <- filterByExpr(y.ab, group=y.ab$samples$BrNum)
    # y.ab <- y.ab[keep,]
    # summary(keep)
    #y.ab <- calcNormFactors(y.ab)
    # Create design matrix
    design <- model.matrix(~factor(Dx) + factor(slide_id) 
                 + factor(Sex), y.ab$samples)

    # Estimate dispersion
    y.ab <- estimateDisp(y.ab, design, trend="none")
    summary(y.ab$common.dispersion)

    plotBCV(y.ab, cex=1)

    fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
    summary(fit.ab$var.prior)

    plotQLDisp(fit.ab, cex=1)

    res <- glmQLFTest(fit.ab, coef=2)
    print(summary(decideTests(res)))

    print(topTags(res))
}
dev.off()
