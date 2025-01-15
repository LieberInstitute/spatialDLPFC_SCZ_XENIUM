library(here)
library(SpatialExperiment)
library(scran)
library(tidyverse)
library(escheR)

spe <- readRDS(here("processed-data", "01_build_spe", "raw_spe_N24.RDS"))

pdf(here("plots", "02_xenium_qc", "01_total_counts.pdf"))
brnums <- unique(spe$BrNum)
for (i in 1:length(brnums)){
    spe_sub <- spe[, spe$BrNum == brnums[i]]
    p <- make_escheR(spe_sub) %>%
        add_fill("total_counts")
    print(p)
}
dev.off()
