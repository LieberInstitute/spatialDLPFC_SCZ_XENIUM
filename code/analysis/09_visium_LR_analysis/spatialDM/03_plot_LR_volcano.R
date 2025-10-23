library(here)
library(ggplot2)
library(patchwork)


me <- "vasc"

# Load the spatialDM results
diff <- read.csv(here("processed-data", 
                              "09_visium_LR_analysis", 
                              paste0("diff_", me, ".csv")), sep="\t")
colnames(diff) <- c("LR_pair", "zscore_diff")

fdr <- read.csv(here("processed-data", 
                              "09_visium_LR_analysis", 
                              paste0("diff_fdr_", me, ".csv")), sep="\t")
diff$fdr <- fdr$X0



pdf(here("plots", "09_visium_LR_analysis", 
            paste0("LR_volcano_", me, ".pdf")), width = 6, height = 5)
p <- ggplot(diff, aes(x=zscore_diff, y=-log10(fdr))) +
  geom_point(alpha=0.5) +
  theme_minimal() +
  xlab("Z-score difference") +
  ylab("-log10(FDR)") +
  ggtitle(paste("Ligand-Receptor interactions -", me)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(0.1), linetype="dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype="dashed", color = "blue")
print(p)
dev.off()
