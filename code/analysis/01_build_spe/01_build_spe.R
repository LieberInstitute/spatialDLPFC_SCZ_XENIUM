library(here)
library(SpatialExperiment)
library(SingleCellExperiment)
library(readxl)
library(tidyverse)


########################################################
# Read in the raw Xenium data into SPE objects and save.
#######################################################

sample_paths <- read.csv(here("raw-data", "experiment_info", "sample_paths.csv"))
sample_info <- read_excel(here("raw-data", "experiment_info", "PNN_64_collection_Summary_final.xlsx"))


sample_info <- sample_info %>%
    filter(Brnumbr %in% sample_paths$BrNum) %>%
    rename(BrNum=Brnumbr)

sample_paths <- merge(sample_paths, sample_info, by = "BrNum")

all_spes <- list()
for(i in 1:nrow(sample_paths)){
    sample_path <- sample_paths$path[i]
    brnum <- sample_paths$BrNum[i]

    print(sprintf("----------%s----------", brnum))

    counts_path <- here(sample_path, "cell_feature_matrix.h5")
    cell_info_path <- here(sample_path, "cells.csv.gz")

    sce <- DropletUtils::read10xCounts(counts_path)
    counts(sce) <- methods::as(DelayedArray::realize(counts(sce)), "dgCMatrix") # Convert to delayed array

    cell_info <- vroom::vroom(cell_info_path)

    colData(sce) <- cbind(colData(sce), cell_info)
    spe <- toSpatialExperiment(sce, spatialCoordsNames = c("x_centroid", "y_centroid"))
    rownames(spe) <- rowData(spe)$Symbol # change rownames to gene symbol

    colnames(spe) <- paste(brnum, rownames(cell_info), sep = "_") # assign donor-specific colnames
    
    # Add in some relevant metadata from the sample info excel file
    spe$BrNum <- brnum
    spe$Dx <- sample_paths$Dx
    spe$CaptureArea<- sample_paths$CaptureArea
    spe$PNN <- sample_paths$PNN
    spe$tear <- sample_paths$tear
    spe$Age <- sample_paths$AGE
    spe$Sex <- sample_paths$SEX

    all_spes[[brnum]] <- spe
    
    # clear memory
    rm(sce)
    rm(cell_info)
    gc()
}

all_spes <- do.call(cbind, all_spes)
saveRDS(all_spes, here("processed-data", "01_build_spe", "raw_spe_N24.RDS"))