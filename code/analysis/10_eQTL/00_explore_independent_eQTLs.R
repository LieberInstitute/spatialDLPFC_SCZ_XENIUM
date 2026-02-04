library(here)
library(biomaRt)
library(tidyverse)
library(SpatialExperiment)
library(escheR)
library(ComplexHeatmap)

spds <- paste0("spd0", 1:7)

spds <- lapply(spds, function(spd){
  eGenes <- read.table(gzfile(paste0("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/eQTL/tqtl_out/", spd, ".gene.map_independent.txt.gz")))
  eGenes$spd <- spd
  return(eGenes)
})

eGenes <- do.call(rbind, spds)

spe <- readRDS(here("processed-data", "01_build_spe", "raw_spe_N24.RDS"))
outlier_ids <- read.csv(here("processed-data", "02_xenium_qc", "outlier_ids.csv"))$x
spe <- spe[, -which(colnames(spe) %in% outlier_ids)]
rm (outlier_ids)
spe <- spe[which(rowData(spe)$Type=="Gene Expression"), ]
print(spe)

#### Banksy parameters ####
lambda <- 0.1
res <- 0.7 # higher = more clusters
use_agf <- FALSE
cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(use_agf), lambda, res)

clusts <- read.csv(here("processed-data", "06_cell_type_clustering", sprintf("banksy_clustering_lambda%s_res%s.csv", lambda, res)))
clusts <- as.data.frame(clusts)
print(head(clusts))
colData(spe)[["Banksy"]] <- as.character(clusts$V1)


# Aggregate to broader groups

clusts <- clusts %>% as.data.frame() %>% 
  mutate(V1 = as.numeric(V1)) %>%
  mutate(annots=case_when(V1 %in% c(6, 8, 1) ~ "Oligo", # 10 and 15 i am not super sure about
                          V1 %in% c(10) ~ "Ambig/Oligo",
                          V1 %in% c(7) ~ "Mic",
                          V1 %in% c(15) ~ "Ambig/In/Endo",
                          V1 %in% c(18) ~ "L5 Ex",
                          V1 %in% c(14) ~ "L6 Ex",
                          V1 %in% c(11) ~ "L4/5 Ex",
                          V1 %in% c(2) ~ "L2/3 Ex",
                          V1 %in% c(17, 16, 4) ~ "Ast",
                          V1 %in% c(5, 13, 3) ~ "Endo",
                          V1 %in% c(12) ~ "CGE",
                          V1 %in% c(9) ~ "MGE",
                       TRUE ~ "NA")) %>%
    mutate(annots_combined = paste(annots, V1, sep="."))
  
colData(spe)[["annots"]] <- clusts$annots
colData(spe)[["annots_combined"]] <- clusts$annots


spds <- read.csv(here("processed-data", "04_label_transfer", "label_transfer_N24_k50_smoothed_labels.csv"))

rownames(spds) <- spds$X 
spds$X <- NULL
colnames(spds) <- c("predictions_smooth")
stopifnot(all(colnames(spe)==rownames(spds)))
colData(spe)$predictions_smooth <- factor(spds$predictions_smooth)

# map spd style labels to actual annotations
domain_annotations <- colData(spe) %>%
  as.data.frame() %>%
  mutate(domain_annotations = case_when(predictions_smooth == "spd07" ~ "L1/M",
                                        predictions_smooth == "spd06" ~ "L2/3",
                                        predictions_smooth == "spd02" ~ "L3/4",
                                        predictions_smooth == "spd05" ~ "L5",
                                        predictions_smooth == "spd03" ~ "L6",
                                        predictions_smooth == "spd01" ~ "WMtz",
                                        predictions_smooth == "spd04" ~ "WM",
                                          TRUE ~ NA))

spe$domain_annotations <- domain_annotations$domain_annotations


#check which of the eGenes are in the Xenium data
eGenes_in_spe <- eGenes[eGenes$phenotype_id %in% rowData(spe)$ID, ]
eGenes_in_spe$gene_name <- rowData(spe)$Symbol[match(eGenes_in_spe$phenotype_id, rowData(spe)$ID)]


panel_markers <- readxl::read_xlsx((here("raw-data", 
        "experiment_info", 
        "Xenium_SHK_celltype_REannot_2025-04-13.xlsx")), sheet=2)

panel_markers <- panel_markers %>%
    as.data.frame() %>%
    mutate(cell_type_updated=case_when(cell_type_updated=="NA" ~ NA,
                                        TRUE ~ cell_type_updated))

eGenes_in_spe <- merge(eGenes_in_spe, panel_markers, by.x="gene_name", by.y="Gene", all.x=TRUE) %>%
  as.data.frame()

# plot eGenes on two example tissue sections
brs_use <- c("Br8667", "Br5973")

# pdf(here("plots", "10_eQTL", "eGenes_in_Xenium_example_sections.pdf"), width=25, height=15)
# for(i in 1:nrow(eGenes_in_spe)){
#   egene <- eGenes_in_spe[i, ]
#   egene_name <- egene$gene_name

#   spe1 <- spe[,spe$BrNum == brs_use[1]]
#   spe1$egene_expr <- counts(spe1)[which(rowData(spe1)$Symbol == egene_name), ]
#   p1 <- make_escheR(spe1) %>%
#         add_fill("egene_expr")%>%
#         add_ground("annots")+
#         scale_fill_gradient(low = "white", high = "black")+
#         ggtitle(paste0(egene_name, " in ", brs_use[1]))
  
#   spe2 <- spe[,spe$BrNum == brs_use[2]]
#   spe2$egene_expr <- counts(spe2)[which(rowData(spe2)$Symbol == egene_name), ]
#   p2 <- make_escheR(spe2) %>%
#         add_fill("egene_expr")%>%
#         add_ground("annots")+
#         scale_fill_gradient(low = "white", high = "black")+
#         ggtitle(paste0(egene_name, " in ", brs_use[2]))
  
#   gridExtra::grid.arrange(p1, p2, ncol=2)

# }


#dev.off()

spe <- scuttle::logNormCounts(spe)
# Make boxplots of expression by cell type for all eGenes in spe
pdf(here("plots", "10_eQTL", "eGenes_in_Xenium_boxplots_by_celltype.pdf"), width=10, height=7)
for(i in 1:nrow(eGenes_in_spe)){
  egene <- eGenes_in_spe[i, ]
  egene_name <- egene$gene_name

  spe$egene_expr <- logcounts(spe)[which(rowData(spe)$Symbol == egene_name), ]
  df_plot <- as.data.frame(colData(spe)) %>%
    dplyr::select(egene_expr, annots, annots_combined, BrNum, domain_annotations)

  p <- ggplot(df_plot, aes(x=annots, y=egene_expr))+
        geom_violin(aes(fill=annots))+
        theme_bw()+
        facet_wrap(~domain_annotations, scales="free_x")+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
        ggtitle(paste0(egene_name, " expression by cell type"))+
        xlab("Cell Type")+
        ylab("Expression (logcounts)")+
        guides(fill=FALSE)

  print(p)
}
dev.off()


plot_heatmap <- function(spe, gene, domain){
  spe$libsize <- colSums(counts(spe))
  col_use <- sprintf("%s_expr", gene)
  col_use <- make.names(col_use)
  colData(spe)[[col_use]] <- counts(spe)[gene, ]
  spe$BrNum <- paste(spe$BrNum, spe$Dx, sep="_")

  agg_df <- colData(spe) %>%
    as.data.frame(check.names=FALSE) %>%
    group_by(BrNum, domain_annotations, annots) %>%
    summarise(pb_libsize=sum(get("libsize")), 
              sum_expr = sum(get(col_use)), 
              .groups="drop")%>%
      mutate(
      norm_expr = log((sum_expr/ pb_libsize)+1)           # normalize by total library size
    )

  agg_df <- agg_df %>%
    mutate(col_name = paste(BrNum, domain_annotations, sep = "_"))

  # Pivot to wide format: rows = cell types, columns = sample_domain
  heatmap_df <- agg_df %>%
    select(annots, col_name, norm_expr) %>%
    pivot_wider(
      names_from = col_name,
      values_from = norm_expr,
      values_fill = 0
    )
  #counts_plt <- logcounts(spe_pseudo)[c("ZNF804A", "XRRA1"), ]
  heatmap_matrix <- as.matrix(heatmap_df[,-1]) # remove 'cell_type' column
  rownames(heatmap_matrix) <- heatmap_df$annots

  parts <- strsplit(colnames(heatmap_matrix), "_")

  diagnosis_vec <- sapply(parts, `[`, 2)      # second field = diagnosis
  domain_vec    <- sapply(parts, `[`, 3)      # third field = spatial domain



  # Desired spatial domain order
  domain_levels <- c("L1/M", "L2/3", "L3/4", "L5", "L6", "WMtz", "WM")
  domain_vec <- factor(domain_vec, levels = domain_levels)

  # Diagnosis order
  diagnosis_levels <- c("NTC", "SCZ")
  diagnosis_vec <- factor(diagnosis_vec, levels = diagnosis_levels)


  col_order <- order(domain_vec, diagnosis_vec)
  heatmap_matrix_ordered <- heatmap_matrix[, col_order]

  domain_vec    <- domain_vec[col_order]
  diagnosis_vec <- diagnosis_vec[col_order]

  # # create colour palette to annotate spds 
  # spatial_domain_vec <- sapply(strsplit(colnames(heatmap_matrix), "_"), `[`, 2)
  #domain_colors <- set_names(Polychrome::palette36.colors(7)[seq.int(7)],
  #domain_levels)

  domain_colors <- c(
    "L1/M" = "#FEAF16",
    "L2/3" = "#3283FE",
    "L3/4" = "#E4E1E3",
    "L5"   = "#16FF32",
    "L6"   = "#F6222E",
    "WMtz" = "#5A5156",
    "WM"   = "#FE00FA"
)

  diagnosis_colors <- c("NTC" = "steelblue", "SCZ" = "firebrick")


  top_ha <- HeatmapAnnotation(
    Domain = domain_vec,
    Diagnosis = diagnosis_vec,
    col = list(
      Domain = domain_colors,
      Diagnosis = diagnosis_colors
    ),
    annotation_height = unit.c(unit(5, "mm"), unit(5, "mm")),
    annotation_legend_param=list(Domain=list(nrow=1, direction="horizontal",
        labels_gp=gpar(fontsize=16),
        title_gp=gpar(fontsize=16)), 
        Diagnosis=list(nrow=1, direction="horizontal",
        labels_gp=gpar(fontsize=16),
        title_gp=gpar(fontsize=16))),
    annotation_name_gp = gpar(fontsize = 16)
  )
  # replace colnames with just the brnum
  colnames(heatmap_matrix_ordered) <- sapply(strsplit(colnames(heatmap_matrix_ordered), "_"), `[`, 1)
  print(colnames(heatmap_matrix_ordered))
  ht <- Heatmap(heatmap_matrix_ordered,
    name = paste("Expression of", gene),
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    bottom_annotation = top_ha,
    heatmap_legend_param = list(
               direction = "horizontal", 
               labels_gp=gpar(fontsize=12),
               title_gp=gpar(fontsize=16)),
    row_names_gp = gpar(fontsize = 16),
    column_title=sprintf("eGene %s found in %s", gene, domain)
  )
  return(ht)
}

# Make heatmaps for all eGenes in spe
pdf(here("plots", "10_eQTL", "eGenes_in_Xenium_heatmaps_by_celltype_domain.pdf"), width=15, height=10)
for(i in 1:nrow(eGenes_in_spe)){
  egene <- eGenes_in_spe[i, ]
  egene_name <- egene$gene_name
  egene_domain <- egene$spd

  ht <- plot_heatmap(spe, egene_name, egene_domain)
  draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend=TRUE)
}

ht <- plot_heatmap(spe, "MAPK3", "Neuropil")
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend=TRUE)

dev.off()




