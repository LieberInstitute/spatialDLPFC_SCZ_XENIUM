library(here)
library(biomaRt)
library(tidyverse)
library(SpatialExperiment)
library(escheR)
library(ComplexHeatmap)
library(pgenlibr)
library(here)
library(data.table)
library(edgeR)
library(jaffelab)


source(here("code", "analysis", "10_eQTL", "00_utils.R"))
eGenes <- read.table(gzfile(paste0("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/eQTL/tqtl_out/", "neuropil", ".gene.map_independent.txt.gz")))


xen_samples <- read.table(here("processed-data", "10_eQTL", "xenium_brnums.keep"))
# read in the genotype info
plink_path <- "/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/genotypes"
pgen_prefix <- here(plink_path, "plink2", "merged_maf05")

genotypes <- prepGenotypes_subset(
  pgen_prefix = pgen_prefix,
  variant_ids = unique(eGenes$variant_id)
)

