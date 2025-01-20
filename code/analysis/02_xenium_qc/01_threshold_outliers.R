library(tidyverse)
library(escheR)
library(here)
library(SpatialExperiment)
library(scuttle)
library(scattermore)
library(gridExtra)


spe <- readRDS(here("processed-data", "01_build_spe", "raw_spe_N24.RDS"))

# compute number of genes detected
#colData(spe)$detected <- Matrix::colSums(counts(spe) > 0)

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
                                


    # Identify the outliers based on these QC metrics
    spe_sub$neg_probe_out <- isOutlier(spe_sub$subsets_negProbe_percent)
    spe_sub$neg_codeword_out <- isOutlier(spe_sub$subsets_negCodeword_percent)
    spe_sub$unassigned_out <- isOutlier(spe_sub$subsets_unassigned_percent)
    spe_sub$mito_out <- isOutlier(spe_sub$subsets_mito_percent)
    spe_sub$detected_out <- isOutlier(spe_sub$detected)
    spe_sub$total_counts_out <- isOutlier(spe_sub$total_counts)


    print("Proportion of negative probe outliers:")
    print(sum(spe_sub$neg_probe_out, na.rm=TRUE))


    pdf(here("plots", "02_xenium_qc", paste0("01_", br_use, ".pdf")))
        # Make some plots based on these metrics 
    make_escheR(spe_sub)%>%
        add_fill("subsets_negProbe_percent")%>%
        add_ground("neg_probe_out")+
        ggtitle(paste(br_use, unique(spe_sub$Dx)[[1]]))+
        geom_scattermore()

    make_escheR(spe_sub)%>%
        add_fill("subsets_negCodeword_percent") %>%
        add_ground("neg_codeword_out")+
        ggtitle(paste(br_use, unique(spe_sub$Dx)[[1]]))+
        geom_scattermore()
    
    make_escheR(spe_sub)%>%
        add_fill("subsets_unassigned_percent")%>%
        add_ground("unassigned_out")+
        ggtitle(paste(br_use, unique(spe_sub$Dx)[[1]]))+
        geom_scattermore()
    
    make_escheR(spe_sub)%>%
        add_fill("subsets_mito_percent")%>%
        add_ground("mito_out")+
        ggtitle(paste(br_use, unique(spe_sub$Dx)[[1]]))+
        geom_scattermore()
    
    make_escheR(spe_sub)%>%
        add_fill("detected")%>%
        add_ground("detected_out")+
        ggtitle(paste(br_use, unique(spe_sub$Dx)[[1]]))+
        geom_scattermore()
    
    make_escheR(spe_sub)%>%
        add_fill("total_counts")%>%
        add_ground("total_counts_out")+
        ggtitle(paste(br_use, unique(spe_sub$Dx)[[1]]))+
        geom_scattermore()
    
    
    # Make some scatterplots of the QC metrics to see how they are related
    mito_1 <- ggplot(colData(spe_sub), aes(x=detected, y=subsets_mito_percent))+
        geom_point(aes(colour=mito_out))+
        geom_smooth(method="lm")+
        ggtitle(paste(br_use, unique(spe_sub$Dx)[[1]]))
    
    mito_2 <- ggplot(colData(spe_sub), aes(x=total_counts, y=subsets_mito_percent))+
        geom_point(aes(colour=mito_out))+
        geom_smooth(method="lm")+
        ggtitle(paste(br_use, unique(spe_sub$Dx)[[1]]))

    grid.arrange(mito_1, mito_2, ncol=2)

    negProbe_1 <- ggplot(colData(spe_sub), aes(x=detected, y=subsets_negProbe_percent))+
        geom_point(aes(colour=neg_probe_out))+
        geom_smooth(method="lm")+
        ggtitle(paste(br_use, unique(spe_sub$Dx)[[1]]))
    
    negProbe_2 <- ggplot(colData(spe_sub), aes(x=total_counts, y=subsets_negProbe_percent))+
        geom_point(aes(colour=neg_probe_out))+
        geom_smooth(method="lm")+
        ggtitle(paste(br_use, unique(spe_sub$Dx)[[1]]))
    
    grid.arrange(negProbe_1, negProbe_2, ncol=2)

    negCode_1 <- ggplot(colData(spe_sub), aes(x=detected, y=subsets_negCodeword_percent))+
        geom_point(aes(colour=neg_codeword_out))+
        geom_smooth(method="lm")+
        ggtitle(paste(br_use, unique(spe_sub$Dx)[[1]]))
    
    negCode_2 <- ggplot(colData(spe_sub), aes(x=total_counts, y=subsets_negCodeword_percent))+
        geom_point(aes(colour=neg_codeword_out))+
        geom_smooth(method="lm")+
        ggtitle(paste(br_use, unique(spe_sub$Dx)[[1]]))
    
    grid.arrange(negCode_1, negCode_2, ncol=2)

    unassigned_1 <- ggplot(colData(spe_sub), aes(x=detected, y=subsets_unassigned_percent))+
        geom_point(aes(colour=unassigned_out))+
        geom_smooth(method="lm")+
        ggtitle(paste(br_use, unique(spe_sub$Dx)[[1]])) 
    
    unassigned_2 <- ggplot(colData(spe_sub), aes(x=total_counts, y=subsets_unassigned_percent))+
        geom_point(aes(colour=unassigned_out))+
        geom_smooth(method="lm")+
        ggtitle(paste(br_use, unique(spe_sub$Dx)[[1]])) 

    grid.arrange(unassigned_1, unassigned_2, ncol=2)
    
    dev.off()                                                    
}
