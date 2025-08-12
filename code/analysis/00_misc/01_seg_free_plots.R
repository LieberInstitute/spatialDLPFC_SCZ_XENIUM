suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(scater)
  library(tidyverse)
  library(escheR)
  library(MoleculeExperiment)
  library(readxl)
  #library(Rarr)
})
Sys.setenv(TMPDIR = here())
##################################################
# Use the MoleculeExperiment package to explore the 
# segmentation-free transcript data from Xenium
##################################################
print("##################### HYP data #####################")
hyp_dir <- "/dcs04/lieber/marmaypag/spatialHYP_LIBD4195/spatial_HYP/xenium_HYP/processed-data/01_xeniumranger1.7-resegment/"
hyp_files <- list.files(hyp_dir)
hyp_files <- hyp_files[grepl("wnuclei", hyp_files)]
hyp_files <- hyp_files[!grepl(".tar.gz", hyp_files)]
hyp_files <- paste0(hyp_dir, hyp_files, "/outs")

for(i in 1:length(hyp_files)){
    sample_path <- hyp_files[[i]]
    #cell_info_path <- here(sample_path, "cells.csv.gz")
    transcripts_path <- here(sample_path, "transcripts.csv")

    #transcripts <- zarr_overview(transcripts_path)

    #cell_info <- vroom::vroom(cell_info_path)
    #transcripts <- vroom::vroom(transcripts_path)
    #transcripts <- arrow::read_parquet(transcripts_path)
    
    me <-  MoleculeExperiment::readXenium(sample_path, keepCols="essential",addBoundaries=c("nucleus","cell"))
    me@molecules$detected$outs <- lapply(me@molecules$detected$outs, FUN=function(g){return(g[which(g$qv>20),])})


}
