library(here)
library(tidyverse)
library(SpatialExperiment)
library(ggplot2)
library(cowplot)

#####################################################################################################
# Reference: https://github.com/LieberInstitute/spatialDLPFC_SCZ/blob/abb428b8ff2765287b3a8f811ef77890b574fe85/code/analysis/06_SpD_prop/plot_spd_proportion_ordered.R#L246
#####################################################################################################

spe <- readRDS(here("processed-data", "01_build_spe", "raw_spe_N24.RDS"))
spds <- read.csv(here("processed-data", "04_label_transfer", "label_transfer_N24_k50_smoothed_labels.csv"))

rownames(spds) <- spds$X 
spds$X <- NULL
colnames(spds) <- c("predictions_smooth")
spe <- spe[,colnames(spe) %in% rownames(spds)]
colData(spe) <- merge(colData(spe), spds, by="row.names", all.x=TRUE)

col_data_df <- as.data.frame(colData(spe))

total_spots_per_sample <- col_data_df |>
  group_by(BrNum) |>
  summarize(total_spots = n())


n_spots_per_spd <- col_data_df |>
  group_by(
    BrNum,
    predictions_smooth
  ) |>
  summarize(n_spots = n()) |>
  ungroup() |>
  left_join(total_spots_per_sample, by = "BrNum") |>
  mutate(proportion = n_spots / total_spots) |>
  select(-total_spots)

demo_df <- colData(spe) |>
    as.data.frame() |>
    select(BrNum, Dx, Sex) |>
    mutate(sample_label = paste0(BrNum, "_", Dx))

## Make plots ----

pdf(here("plots", "04_label_transfer", "spd_props.pdf"))
ret_p <- n_spots_per_spd |>
  mutate(PRECAST_07 = factor(predictions_smooth, levels = sprintf("spd%02d", c(7, 6, 2, 5, 3, 1, 4)))) |>
  ggplot() +
  geom_bar(
    aes(x = BrNum, y = proportion, fill = PRECAST_07),
    stat = "identity",
    position = "fill"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  labs(
    x = "BrNum",
    y = "Proportion",
    fill = "Spds"
  ) +
  scale_fill_manual(
    name = "Spatial Domain",
    # change fill palette
    values = set_names(
      Polychrome::palette36.colors(7)[seq.int(7)],
      unique(spe$predictions_smooth) |> sort()
    ))+
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
    )

dx_p <- demo_df |>
  ggplot() +
  geom_point(
    aes(
      x = sample_label,
      y = 1,
      color = Dx,
      shape = Sex
    ),
    size = 3
  ) +
  scale_color_manual(
    values = c("NTC" = "blue", "SCZ" = "red"),
    guide = "none"
  )+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )+
    # theme_void() +
    theme(
      # legend.position = "none",
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.background = element_rect(fill = "transparent", color = NA),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA)
    )
plot_grid(ret_p, dx_p,
  align = "v",
  ncol = 1,
  rel_heights = c(0.8, 0.2))
dev.off()




