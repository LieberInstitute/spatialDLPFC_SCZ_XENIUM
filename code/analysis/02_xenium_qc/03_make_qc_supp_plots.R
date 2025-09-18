library(Banksy)
library(SummarizedExperiment)
library(SpatialExperiment)
library(Seurat)
library(scater)
library(cowplot)
library(ggplot2)
library(here)
library(tidyverse)
library(escheR)
library(patchwork)


##########################################################################################
# Use BANKSY to do integration-free spatial domain analysis as a first pass
##########################################################################################

# Read in the data
spe <- readRDS(here("processed-data", "01_build_spe", "raw_spe_N24.RDS"))
outlier_ids <- read.csv(here("processed-data", "02_xenium_qc", "outlier_ids.csv"))$x

spe$is_outlier <- ifelse(colnames(spe) %in% outlier_ids, "Discarded", "Kept")
brnums <- unique(spe$BrNum)

is_neg <- stringr::str_detect(rownames(spe), "^NegControlProbe")
is_neg2 <- stringr::str_detect(rownames(spe), "^NegControlCodeword")
is_unassigned <- stringr::str_detect(rownames(spe), "^Unassigned")

# Add QC metrics for plotting
spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(negProbe = is_neg,
                                                        negCodeword = is_neg2,
                                                        unassigned = is_unassigned))
cols_use <- names(colData(spe))[str_detect(names(colData(spe)), "_percent$")]
cols_use <- cols_use[!grepl("subsets_blank_percent", cols_use)]
cols_use <- cols_use[!grepl("subsets_anti_percent", cols_use)]
cols_use <- cols_use[!grepl("depr", cols_use)]
cols_use <- cols_use[!grepl("log10_subsets_unassigned_percent", cols_use)]
cols_use <- cols_use[!grepl("unassigned", cols_use)]

cols_use <- c(cols_use, "detected", "total_counts")


plist <- list()
for (i in 1:length(brnums)){
    br_use <- brnums[[i]]
    spe_use <- spe[,spe$BrNum==br_use]
    p <- make_escheR(spe_use, y_reverse=FALSE) %>%
        add_fill("is_outlier", point_size=1)+
        ggtitle(paste(br_use, unique(spe_use$Dx))) +
        theme(plot.title = element_text(hjust = 0.5))+
        scale_fill_manual(values=c("Kept"="grey80", "Discarded"="deeppink"))
    plist[[i]] <- p
                                                        

     
}

png(here("plots", "02_xenium_qc", "qc_outlier_spots.png"), width=1000, height=1200)
wrap_plots(plist, ncol=4, guides="collect") & 
            guides(fill=guide_legend(title="Cells", override.aes = list(size=5))) &
            theme(legend.position="bottom", 
                legend.title=element_text(size=24), 
                legend.text=element_text(size=24),
                title=element_text(size=18))
dev.off()

spe$BrNum <- paste(spe$BrNum, spe$Dx, sep="_")
colData(spe)$suffix <- sub(".*(_NTC|_SCZ)$", "\\1", colData(spe)$BrNum)
colData(spe)$suffix <- factor(colData(spe)$suffix, levels = c("_NTC", "_SCZ"))

colData(spe)$BrNum<- factor(
  colData(spe)$BrNum,
  levels = unique(colData(spe)$BrNum[order(colData(spe)$suffix)])
)

p_detected <- plotColData(spe, x="BrNum", y="detected", colour_by="is_outlier")+
        ggtitle("Detected")+
        scale_colour_manual(values=c("Kept"="grey80", "Discarded"="deeppink"))

p_total_counts <- plotColData(spe, x="BrNum", y="total_counts", colour_by="is_outlier")+
        ggtitle("Total counts")+
        scale_colour_manual(values=c("Kept"="grey80", "Discarded"="deeppink"))

p_negProbe <- plotColData(spe, x="BrNum", y="subsets_negProbe_percent", colour_by="is_outlier")+
        ggtitle("Negative control probe percentage")+
        scale_colour_manual(values=c("Kept"="grey80", "Discarded"="deeppink"))

p_negCodeword <- plotColData(spe, x="BrNum", y="subsets_negCodeword_percent", colour_by="is_outlier")+
        ggtitle("Negative control codeword percentage")+
        scale_colour_manual(values=c("Kept"="grey80", "Discarded"="deeppink"))

violin_ps <- list(p_detected, p_total_counts, p_negProbe, p_negCodeword)
# Make violin plots for the QC metrics
png(here("plots", "02_xenium_qc", "qc_outlier_violins.png"), width=1200, height=800)
wrap_plots(violin_ps, ncol=2, guides="collect") & 
    guides(colour=guide_legend(title="Cells", override.aes = list(size=5))) &
    theme(legend.position="bottom", 
          legend.title=element_text(size=24), 
          legend.text=element_text(size=24),
          title=element_text(size=18),
          axis.text.x=element_text(angle=45, hjust=1)) 
dev.off()
