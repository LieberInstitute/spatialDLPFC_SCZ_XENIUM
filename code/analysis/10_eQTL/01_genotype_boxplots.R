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



ds <- "spd07"
# subset spe_pseudo to just the samples from ds
spe_pseudo_ds <- spe_pseudo[,spe_pseudo$predictions_smooth == ds]
colnames(spe_pseudo_ds) <- spe_pseudo_ds$BrNum
# make tmm assay for spe_pseudo
norm <- calcNormFactors(DGEList(counts(spe_pseudo_ds)), method = "TMM")
log2cpm <- cpm(norm, log = TRUE, prior.count = 1)
assays(spe_pseudo_ds)$tmm <- log2cpm


# residualize MAPK3 expression 

# read in covariates

tqtl_in <- here("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100", "processed-data", "eQTL/tqtl_in")
covar_file <- file.path(tqtl_in, sprintf("%s.gene.covars.txt", ds))

covars <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
covars <- as.data.frame(t(covars))
colnames(covars) <- covars[1,]
covars <- covars[-1,]
covars <- covars[rownames(covars) %in% xen_samples, ] # subset to xenium samples





gene_expr <- assays(spe_pseudo_ds)$tmm[rownames(spe_pseudo)=="KANSL1-AS1",]

# add in the Xenium slide information as dummy variables
slide_df <- model.matrix(~slide_id-1, colData(spe_pseudo_ds))
slide_df <- as.data.frame(slide_df)


covars <- covars %>%
  select(-grep("expr|slide", colnames(covars)))%>% # remove the gene expression PCs from visium% # remove the slide information because it's irrelevant for Xenium
  mutate(DXSCZ=as.numeric(DXSCZ)) %>%
  mutate(sexM = as.numeric(sexM)) %>%
  mutate(across(starts_with("snp"), as.numeric)) %>%
  mutate("(Intercept)"= 1) %>%
  mutate(age = as.numeric(age))

#covars <- merge(covars, slide_df, by = "row.names")# %>%
#   mutate(across(starts_with("slide"), factor))

# Compute our own gene expression PCs for the Xenium data and add those to the covariates
mat_z <- t(scale(t(log2cpm), center = TRUE, scale = TRUE))
mat_z[!is.finite(mat_z)] <- 0
assays(spe_pseudo_ds)$zscore <- mat_z
set.seed(1301)
pc <- prcomp(t(mat_z), center = FALSE, scale. = FALSE)
expr_pcs <- pc$x
colnames(expr_pcs) <- paste0("exprPC", seq_len(ncol(expr_pcs)))

covars <- rownames_to_column(covars, var = "Row.names")

# Use the first 8 PCs similar to the visium analysis
expr_pcs_use <- as.data.frame(expr_pcs[,1:8])
expr_pcs_use <- rownames_to_column(expr_pcs_use, var = "Row.names")

# Merge the expr PCs into the covariate df
covars <- merge(covars, expr_pcs_use, by="Row.names")
gene_expr <- gene_expr[covars$Row.names]
gene_expr <- as.vector(gene_expr)
covars <- column_to_rownames(covars, var = "Row.names")

# Get rid of one of the slide variables to avoid rank issues
covars <- covars[, -grep("slide", colnames(covars))[1]]
covars <- covars %>%
  mutate(across(starts_with("slide"), as.numeric))

# Convert the vector to a 1-row matrix
gene_expr_mat <- matrix(gene_expr, nrow = 1)
colnames(gene_expr_mat) <- rownames(covars)
# use jaffelab::cleaningY to residualize the gene expression for the covariates
resid_expr <- cleaningY(gene_expr_mat, mod = as.matrix(covars), P = 1)

# Now get the dosage info
gene_id <- rowData(spe_pseudo_ds)$ID[rownames(spe_pseudo_ds)=="KANSL1-AS1"]
variant_id <- eGenes$variant_id[eGenes$phenotype_id == gene_id & eGenes$spd == ds]
dosage <- genotypes_xen$dosage[, variant_id]

dosage <- dosage[colnames(resid_expr)]

plot_df <- cbind(dosage, resid_expr=as.vector(resid_expr))
plot_df <- as.data.frame(plot_df)
stopifnot(all.equal(rownames(plot_df), rownames(covars)))

plot_df$Dx <- ifelse(covars$DXSCZ == 1, "SCZ", "Control")

pdf(here("plots", "10_eQTL", "spd07_KANSL1-AS1_dosage_boxplot.pdf"), width = 4, height = 4)
ggplot(plot_df, aes(x = factor(dosage), y = resid_expr)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, aes(colour=Dx)) +
  labs(x = "Dosage", y = "Residualized Expression of KANSL1-AS1") +
  theme_bw() +
  theme(legend.position = "top")
dev.off()
# actually, we can't do neuropil because we don't have things pseudobulked for the IF domains
# since we don't have IF domain labels in 

