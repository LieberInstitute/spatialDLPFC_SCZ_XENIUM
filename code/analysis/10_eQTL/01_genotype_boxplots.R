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

# read in the independent eqtls for each spd
spds <- paste0("spd0", 1:7)
spds <- lapply(spds, function(spd){
  eGenes <- read.table(gzfile(paste0("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/eQTL/tqtl_out/", spd, ".gene.map_independent.txt.gz")))
  eGenes$spd <- spd
  return(eGenes)
})

eGenes <- do.call(rbind, spds)

# read in the xenium object
spe_pseudo <- readRDS(here("processed-data/05_differential_expression/spe_pseudo_donor_domain_spaTransfer_k50_smoothed_predictions.rds"))


# subset eGenes to just the genes that are in the Xenium data
eGenes <- eGenes[eGenes$phenotype_id %in% rowData(spe_pseudo)$ID,]
# read in the genotype info
plink_path <- "/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/genotypes"

# Define the file path (without the extension)
pgen_prefix <- here(plink_path, "plink2", "merged_maf05")

prepGenotypes_subset <- function(pgen_prefix, variant_ids, gwas_vars = NULL) {
  # Load pvar
  pvar_file <- paste0(pgen_prefix, ".pvar")
  pvar <- fread(pvar_file)
  setnames(pvar, "#CHROM", "CHROM")
  pvar[, var_idx := .I]  # 1-based index
  
  # Keep only variants in variant_ids
  pvar_sub <- pvar[ID %in% variant_ids]
  if (nrow(pvar_sub)==0) stop("No matching variants in PLINK2 .pvar!")
  
  # Load psam
  psam <- fread(paste0(pgen_prefix, ".psam"))
  samples <- if ("#IID" %in% names(psam)) psam$`#IID` else psam$IID


  # Open pgen
  pgen <- NewPgen(paste0(pgen_prefix, ".pgen"))
  dosage <- ReadList(pgen, pvar_sub$var_idx, meanimpute = FALSE)
  rm(pgen); gc()
  
  rownames(dosage) <- samples
  colnames(dosage) <- pvar_sub$ID
  
  list(
    dosage = dosage,
    pvar = pvar_sub,
    samples = samples
  )
}

genotypes <- prepGenotypes_subset(
  pgen_prefix = pgen_prefix,
  variant_ids = unique(eGenes$variant_id)
)

# subset the genotypes matrix to the samples in Xenium
xen_samples <- unique(spe_pseudo$BrNum)

genotypes_xen <- list(dosage = genotypes$dosage[genotypes$samples %in% xen_samples, ],
                      pvar = genotypes$pvar,
                      samples = genotypes$samples[genotypes$samples %in% xen_samples])


# make tmm assay for spe_pseudo
norm <- calcNormFactors(DGEList(counts(spe_pseudo)), method = "TMM")
log2cpm <- cpm(norm, log = TRUE, prior.count = 1)
assays(spe_pseudo)$tmm <- log2cpm


# residualize MAPK3 expression 

# read in covariates
ds <- "spd07"
tqtl_in <- here("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100", "processed-data", "eQTL/tqtl_in")
covar_file <- file.path(tqtl_in, sprintf("%s.gene.covars.txt", ds))

covars <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
covars <- as.data.frame(t(covars))
colnames(covars) <- covars[1,]
covars <- covars[-1,]
covars <- covars[rownames(covars) %in% xen_samples, ] # subset to xenium samples

# remove the slide information because it's irrelevant for Xenium
covars <- covars %>% select(-grep("slide", colnames(covars)))

# add in the Xenium slide information as dummy variables
slide_df <- model.matrix(~slide_id-1, colData(spe_pseudo_ds))

covars <- merge(covars, slide_df, by = "row.names")


# subset spe_pseudo to just the samples from ds
spe_pseudo_ds <- spe_pseudo[,spe_pseudo$predictions_smooth ==ds]
colnames(spe_pseudo_ds) <- spe_pseudo_ds$BrNum

gene_expr <- assays(spe_pseudo_ds)$tmm[rownames(spe_pseudo)=="KANSL1-AS1",]

# use jaffelab::cleaningY to residualize the gene expression for the covariates
resid_expr <- jaffelab::cleaningY(gene_expr, mod = covars, P = 1)

# actually, we can't do neuropil because we don't have things pseudobulked for the IF domains
# since we don't have IF domain labels in Xenium....