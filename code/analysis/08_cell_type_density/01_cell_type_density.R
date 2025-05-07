suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(scater)
  library(tidyverse)
  library(ggrepel)
  library(sessioninfo)
  library(alphahull)
  library(escheR)
  library(sf)
  library(readxl)
  library(SpatialFeatureExperiment)
  library(Voyager)
  library(rjson)
  library(concaveman)
})

####################################################################################################
# Compute and plot the proportion of each cell type in the spatial domains
####################################################################################################


spe <- readRDS(here("processed-data", "07_cell_type_de", 
            "cleaned_spe_N24_with_cell_type_and_spds.RDS"))

spe$cell_types <- colData(spe)[["Banksy-clust_M0_lam0.1_k50_res0.7-cell-types"]]

brnums <- unique(spe$BrNum)
sample_paths <- read.csv(here("raw-data", "experiment_info", "sample_paths.csv"))
sample_info <- read_excel(here("raw-data", "experiment_info", "PNN_64_collection_Summary_final.xlsx"))


sample_info <- sample_info %>%
    filter(Brnumbr %in% sample_paths$BrNum) %>%
    rename(BrNum=Brnumbr)

sample_paths <- merge(sample_paths, sample_info, by = "BrNum")


densities <- list()
pdf(here("plots", "08_cell_type_density", "test_density.pdf"),
        height=30, width=30)
for(i in 1:length(brnums)){
    sample_path <- sample_paths$path[i]
    br_use <- sample_paths$BrNum[i]

    spe_sub <- spe[,spe$BrNum == br_use]
    cell_poly_path <- here::here(sample_path, "cell_boundaries.csv.gz")
    cell_poly <- vroom::vroom(cell_poly_path)

    names(cell_poly)[1] <- "ID"
    cells_sf <- df2sf(cell_poly, c("vertex_x", "vertex_y"), geometryType = "POLYGON")

    if(!all(sf::st_is_valid(cells_sf))){
      ind_invalid <- !sf::st_is_valid(cells_sf)
      cells_sf[ind_invalid,] <- nngeo::st_remove_holes(sf::st_buffer(cells_sf[ind_invalid,], 0))
    }
    
    cells_sf_union <- st_union(cells_sf)
    polygon_lines <- st_cast(cells_sf_union, "MULTILINESTRING")  # convert boundaries to lines
    points <- st_cast(polygon_lines, "POINT")                    # convert lines to points
    points_sf <- st_as_sf(data.frame(geometry = points), crs = st_crs(cells_sf_union))
    
    concave <- concaveman(points_sf)

    plot(st_geometry(cells_sf), col = 'lightblue')
    plot(concave, add = TRUE, border = 'darkgreen', lwd = 2)

    area <- st_area(concave)
    cell_type_densities <- as.data.frame(table(spe_sub$cell_types)/area)
    cell_type_densities$BrNum <- br_use
    colnames(cell_type_densities) <- c("cell_type", "density", "BrNum")
    cell_type_densities$area <- area
    cell_type_densities$Dx <- unique(spe_sub$Dx)

   
    densities[[i]] <- cell_type_densities


}
dev.off()


all_densities <- do.call(rbind, densities)
all_densities <- as.data.frame(all_densities)

pdf(here("plots", "08_cell_type_density", "whole_tissue_density_boxplots.pdf"))
ggplot(all_densities, aes(x=Dx, y=density))+
  geom_boxplot()+
  facet_wrap(~cell_type, scales="free")
dev.off()
