suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
  library(escheR)
  library(readxl)
  library(edgeR)
  library(cacoa)
  library(cowplot)
  library(scran)
  library(bluster)
  library(dendextend)
})

####################################################################################################
# Based on 03_compositional_analysis.R, we saw that 
# L6Ex cells are enriched in L2/3 in SCZ compared to NTC
# Explore the L6Ex cells in L2/3.
########################################################################
spe <- readRDS(here("processed-data", "07_cell_type_de", 
            "cleaned_spe_N24_with_cell_type_and_spds.RDS"))

spe$cell_types <- colData(spe)[["Banksy-clust_M0_lam0.1_k50_res0.7-cell-types"]]
brnums <- unique(spe$BrNum)
layers <- unique(spe$domain_annotations)

spe <- scuttle::logNormCounts(spe)
spe$is_l6ex <- spe$cell_types == "L6 Ex"
spe$isLayer <- spe$domain_annotations %in% c("L2/3")

#---------------------------------------
# Try a k=2 clustering on all L6 Ex cells in all layers/tissues to see if they 
# naturally separate out into those in L2/3 and those in L6
#---------------------------------------


# pdf(here("plots", "08_cell_type_density", 
#          "l6ex_on_tissue.pdf"))
# for (i in 1:length(brnums)){
#   brnum <- brnums[i]
#   spe_sub <- spe[, spe$BrNum == brnum]
  
#   p <- make_escheR(spe_sub) %>%
#       add_ground("isLayer", stroke=0.75) %>%
#       add_fill("is_l6ex", point_size=0.8) +
#       ggtitle(paste0(brnum, " ", spe_sub$Dx[1]))
#   print(p)
# }
# dev.off()

# Make a boxplot of the proportion of L6 Ex cells in L2/3  

spe_l23 <- spe[,spe$domain_annotations == "L2/3"]
spe_l23$Br_Dx <- paste(spe_l23$BrNum, spe_l23$Dx, sep="_")

pdf(here("plots", "08_cell_type_density", 
         "l6ex_boxplot.pdf"))
props <- table(spe_l23$is_l6ex, spe_l23$Br_Dx)
props <- props/colSums(props)
props <- as.data.frame(props)
colnames(props) <- c("is_l6ex", "BrNum", "proportion")

props <- props %>%
  filter(is_l6ex == TRUE) %>%
  mutate(Dx= ifelse(grepl("SCZ", BrNum), "SCZ", "NTC"))%>%
  mutate(Dx=factor(Dx, levels=c("NTC", "SCZ"))) 

p <- ggplot(props, aes(x=BrNum, y=proportion, fill=Dx))+
  geom_bar(stat="identity")+
  ggtitle("Proportion of cells in L2/3 that are L6 Ex")+
  scale_fill_manual(values=c("NTC"="steelblue", "SCZ"="firebrick"))+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45, hjust=1))
  
print(p)

dev.off()
print(props)
# print minimum and max proportion of L6 Ex cells in L2/3 for SCZ donors
props_scz <- props %>% filter(Dx=="SCZ")
print(paste0("Min proportion of L6 Ex cells in L2/3 for SCZ donors: ", min(props_scz$proportion)))
print(paste0("Max proportion of L6 Ex cells in L2/3 for SCZ donors: ", max(props_scz$proportion)))

# spe_l6ex <- spe[,spe$cell_types == "L6 Ex"]
# spe_l6ex <- scater::runPCA(spe_l6ex, ncomponents=50)
# spe_l6ex$isLayer <- spe_l6ex$domain_annotations %in% c("L2/3")


# pdf(here('plots', '08_cell_type_density', 
#          'l6ex_pca.pdf'), width=10, height=10)
# scater::plotPCA(spe_l6ex, colour_by="isLayer", ncomponents=5)
# dev.off()

# hclust_path <- here("processed-data", "08_cell_type_density",
#                 "l6ex_hclust.RDS")
# if (!file.exists(hclust_path)){
#   set.seed(1345)
#   hclust_l6ex <- clusterCells(spe_l6ex, use.dimred="PCA",
#       BLUSPARAM=HclustParam(method="ward.D2", cut.dynamic=TRUE,
#           cut.params=list(minClusterSize=10, deepSplit=1)))
#   print(hclust_l6ex)

#   saveRDS(hclust_l6ex, hclust_path)
# }else{
#   hclust_l6ex <- readRDS(hclust_path)
#   spe_l6ex$hclust <- factor(hclust_l6ex)

#   # Merge spe_l6ex colData back with the big spe object
#   clusts <- colData(spe_l6ex)[, c("hclust", "cell_id")]
#   clusts <- merge(clusts, colData(spe), by="cell_id", all.y=TRUE)
#   clusts$hclust <- ifelse(is.na(clusts$hclust), "Other", clusts$hclust)
  
#   spe$hclust <- clusts$hclust
#   spe$isLayer <- spe$domain_annotations %in% c("L2/3")

#   # Create colour palette for hclust
#   hclust_colors <- c(
#     "1" = "#A6CEE3",
#     "2" = "#1F78B4",
#     "3" = "#B2DF8A",
#     "4" = "#33A02C",
#     "5" = "#FB9A99",
#     "6" = "#E31A1C",
#     "7" = "#FDBF6F",
#     "8" = "#FF7F00",
#     "9" = "#CAB2D6",
#     "10" = "#6A3D9A",
#     "11" = "#B15928",
#     "12" = "#FF69B4", 
#     "Other" = "#D3D3D3", # Grey for other clusters
#     "13" = "#FFB6C1", # Light pink for additional cluster,
#     "14" = "#8B4513", # Saddle brown for another cluster,
#     "15" = "#FFD700", # Gold for another cluster
#     "16" = "#FF4500", # Orange red for another cluster
#     "17" = "#2E8B57", # Sea green for another cluster
#     "18" = "#8A2BE2"  # Blue violet for another cluster
# )
#   pdf(here("plots", "08_cell_type_density", 
#            "l6ex_hclust.pdf"), width=15, height=15)
#   for (i in 1:length(brnums)){
#     brnum <- brnums[i]
#     spe_sub <- spe[, spe$BrNum == brnum]
    
#     # Create colour palette for hclust
#     p <- make_escheR(spe_sub, point_size=1.1, stroke=0.75) %>%
#         add_ground("isLayer")%>%
#         add_fill("hclust")+
#         scale_fill_manual(values=hclust_colors) +
#         ggtitle(paste0(brnum, " ", spe_sub$Dx[1])) 
#   print(p)

#   }
#   dev.off()
# }

# tree <- hclust_l6ex$objects$hclust
# tree$labels <- seq_along(tree$labels)
# dend <- as.dendrogram(tree, hang=0.1)

# combined.fac <- paste0(spe_l6ex$domain_annotations == "L2/3")

# labels_colors(dend) <- c(
#     "TRUE"="blue",
#     "FALSE"="red"
   
# )[combined.fac][order.dendrogram(dend)]

# pdf("test-cluster_l6ex.pdf", width=10, height=10)
# plot(dend)

# dev.off()

#spe_l6ex$hclust <- factor(hclust_l6ex)
