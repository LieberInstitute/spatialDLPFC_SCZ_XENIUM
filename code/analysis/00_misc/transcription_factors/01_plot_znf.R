suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
  library(escheR)
})

####################################################################################################
# Pseudobulk the SPE to the donor-domain level using the smoothed labels from label transfer.
# Reference: https://github.com/LieberInstitute/spatialDLPFC_SCZ/blob/main/code/analysis/pseudobulk_dx/create_pb_data.R
####################################################################################################

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
                          V1 %in% c(12) ~ "In: VIP, LAMP5",
                          V1 %in% c(9) ~ "In: SST, PVALB",
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
  mutate(domain_annotations = case_when(predictions_smooth == "spd07" ~ "L1",
                                        predictions_smooth == "spd06" ~ "L2/3",
                                        predictions_smooth == "spd02" ~ "L3/4",
                                        predictions_smooth == "spd05" ~ "L5",
                                        predictions_smooth == "spd03" ~ "L6",
                                        predictions_smooth == "spd01" ~ "WMtz",
                                        predictions_smooth == "spd04" ~ "WM",
                                          TRUE ~ NA))

spe$domain_annotations <- domain_annotations$domain_annotations
# Add some slide and date related metadata
slide_id <- unlist(lapply(strsplit(spe$Sample, split="/"), "[", 9))
slide_id <- unlist(lapply(strsplit(slide_id, split="__Br"), "[", 1))
spe$slide_id <- slide_id

run_date <- unlist(lapply(strsplit(spe$Sample, split="/"), "[", 7))
spe$run_date <- run_date

gene_plot <- "ZNF804A"
brnums <- unique(spe$BrNum)

spe$is_neuron <- ifelse(spe$annots %in% c("L2/3 Ex", "L4/5 Ex", "L5 Ex", "L6 Ex", "In: SST, PVALB", "In: VIP, LAMP5"), TRUE, FALSE)

pdf(here("plots", "00_misc", "transcription_factors", "ZNF804A_counts_logcounts.pdf"), 
            width=24, height=18)
for(i in 1:length(brnums)){ # change later
    sample <- brnums[i]
    spe_sub <- spe[, spe$BrNum == sample]
    spe_sub <- scuttle::logNormCounts(spe_sub)
    
    colData(spe_sub)[[paste0(gene_plot, "_counts")]] <- counts(spe_sub)[gene_plot, ]
    colData(spe_sub)[[paste0(gene_plot, "_logcounts")]] <- logcounts(spe_sub)[gene_plot, ]

    p_counts <- make_escheR(spe_sub, y_reverse=FALSE) %>%
        add_fill(paste0(gene_plot, "_counts"), point_size=3) +
        #add_ground("is_neuron", stroke=0.75)+
        scale_fill_continuous(low="white", high="black")+
        #scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "transparent"))+
        ggtitle(paste0("Counts of ", gene_plot, " in sample ", sample))
     
    p_logcounts <- make_escheR(spe_sub, y_reverse=FALSE) %>%
        add_fill(paste0(gene_plot, "_logcounts"), point_size=3) +
        #add_ground("is_neuron", stroke=0.75)+
        scale_fill_continuous(low="white", high="black")+
        #scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "transparent"))+
        ggtitle(paste0("Counts of ", gene_plot, " in sample ", sample))

    print(patchwork::wrap_plots(p_counts, p_logcounts, ncol=2))
}
dev.off()


############# Make bar chart for gene_plot positive proportion ###########################
colData(spe)[[paste0(gene_plot, "_counts")]] <- counts(spe)[gene_plot, ]
spe$is_gene_plot_pos <- ifelse(spe[[paste0(gene_plot, "_counts")]] > 0, TRUE, FALSE)

df <- colData(spe) %>%
    as.data.frame() %>%
    mutate(gene_plot_info = case_when(is_gene_plot_pos == TRUE ~ annots,
                                    is_gene_plot_pos == FALSE ~ "Not expressed")) %>%
      group_by(domain_annotations) %>%
      mutate(domain_total = n()) %>%
      group_by(domain_annotations, annots) %>%
      summarise(
            n_total = n(),
            n_gene_pos = sum(is_gene_plot_pos),
            prop_gene_pos = n_gene_pos / first(domain_total),
            prop_gene_neg = sum(!is_gene_plot_pos)/first(domain_total),
            .groups = "drop"
        )

df <- df %>% pivot_longer(cols=c(prop_gene_pos, prop_gene_neg), 
                          names_to="has_expression", 
                          values_to="prop_gene_pos") %>%
    mutate(has_expression = ifelse(has_expression == "prop_gene_pos", "ZNF804A+", "ZNF804A-"))

pdf(here("plots", "00_misc", "transcription_factors", "ZNF804A_proportion.pdf"), 
            width=10, height=8)

p <- ggplot(df, aes(x=domain_annotations, y=prop_gene_pos, fill=annots)) +
    geom_bar(stat="identity", position="fill") +
    #scale_fill_manual(values=c("TRUE"="red", "FALSE"="grey"), name=paste0(gene_plot, " positive")) +
    labs(x="Spatial Domain", y="Proportion of cells") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    facet_wrap(~has_expression)

print(p)
dev.off()



