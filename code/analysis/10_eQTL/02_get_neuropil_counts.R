#install.packages("pliman")
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
  library(pliman) # for otsu thresholding
 })
#########################################################
# Read in the transcript file for the given sample index
# and create 75x75 square bins. Calculate the 
# proportion of MOBP transcript in each of the bins
# and plot them on the tissue to see where the WM region is
#########################################################

args <- commandArgs(trailingOnly=TRUE)
sample_idx <- as.numeric(args[[1]])
print(sample_idx)

files_dir <- "/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/raw-data/xenium"
xen_runs <- list.files(files_dir, full.names=TRUE)
tissue_dirs <- lapply(xen_runs, list.files, full.names=TRUE)
tissue_dirs <- unlist(tissue_dirs)
transcript_files <- here(tissue_dirs, "transcripts.parquet")


# Extract Br number from file names
brnums <- gsub(".*Br(\\d+).*", "\\1", tissue_dirs)
brnums <- paste0("Br", brnums)



pdf(here("plots", "10_eQTL",sprintf("02_%s_neuropil_MAPK3.pdf.pdf", brnums[[sample_idx]])),
    width = 21, height = 6)

    print(brnums[[sample_idx]])
    transcripts <- arrow::read_parquet(transcript_files[sample_idx])
    transcripts <- st_as_sf(transcripts, 
    coords = c("x_location", "y_location"), crs=NA)
    
    # Create hexagonal bins
    sfc <- st_make_grid(st_bbox(transcripts), square=TRUE, cellsize = 75)
    grid_sf <- st_sf(geometry = sfc)
    grid_sf$grid_id <- seq_len(nrow(grid_sf))
    

    # Calculate mean MOBP counts in each bin

    mobp_transcripts <- transcripts %>%
        filter(feature_name == "MOBP" & qv > 20)


    # --- 1. Count MOBP transcripts per bin ---
    mobp_idx <- st_intersects(grid_sf, mobp_transcripts)
    #print(head(mobp_idx))
    MOBP_counts <- lengths(mobp_idx)

    # --- 1a. Count MAPK3 transcripts per bin, OUTSIDE the cell ---
    MAPK3_counts <- transcripts %>%
        filter(feature_name == "MAPK3" & qv > 20) %>%
        filter(cell_id == "UNASSIGNED")
    mapk3_idx <- st_intersects(grid_sf, MAPK3_counts)
    MAPK3_counts <- lengths(mapk3_idx)

    # --- 2. Count ALL transcripts per bin (library size) ---
    all_idx <- st_intersects(grid_sf, transcripts)
    total_counts <- lengths(all_idx)

    # --- 3. Compute normalized counts ---
    grid_sf <- grid_sf %>%
    mutate(
        MOBP_counts = MOBP_counts,
        total_counts = total_counts,
        MOBP_norm = ifelse(total_counts > 0, MOBP_counts / total_counts, 0),
        MAPK3_counts = MAPK3_counts
    )

    print(summary(grid_sf$MOBP_norm))

    # --- 4. Otsu thresholding to identify high MOBP bins ---
    otsu_thresh <- pliman::otsu(grid_sf$MOBP_norm)
    print(otsu_thresh)

    # Also compute Otsu threshold for in-tissue bins based on total counts
    in_tissue_otsu_thresh<- pliman::otsu(grid_sf$total_counts)
    print(in_tissue_otsu_thresh)
 


    p <- ggplot() +
        geom_sf(data = grid_sf, aes(fill = MOBP_norm), color = "grey70", size = 0.1) +
        scale_fill_viridis_c(option = "plasma") +
        coord_sf() +
        theme_minimal() +
        labs(title = sprintf("MOBP transcript counts per 75x75 bin for %s", brnums[[sample_idx]]),
            fill = "MOBP normalized counts")
    
    grid_sf <- grid_sf %>%
    mutate(
        passes_otsu = MOBP_norm > otsu_thresh,
        in_tissue = total_counts > in_tissue_otsu_thresh | passes_otsu
    )

    

    # saveRDS(grid_sf, here("processed-data", "QC", "07_segmentation_free_analysis", "SCZ", "03_segment_WM",
    #     sprintf("%s-MOBP-bins.rds", brnums[[sample_idx]])))

    p2 <- ggplot(grid_sf) +
    geom_sf(aes(fill = passes_otsu), color = "grey70", size = 0.1) +
    scale_fill_manual(
        values = c("FALSE" = "grey90", "TRUE" = "orange"),
        labels = c("Below", "Above"),
        name = sprintf("MOBP > %.3f", otsu_thresh)
    ) +
    coord_sf() +
    theme_minimal() +
    labs(
        title = sprintf("MOBP-high bins (> Otsu threshold %.3f)", otsu_thresh),
        subtitle = sprintf("Total %d bins above threshold", sum(grid_sf$passes_otsu))
    )

    p3 <- ggplot(grid_sf) +
        geom_sf(aes(fill = in_tissue), color = "grey70", size = 0.1) +
        scale_fill_manual(
            values = c("FALSE" = "transparent", "TRUE" = "grey70"),
            labels = c("Below", "Above"),
            name = sprintf("total_counts > %.3f", in_tissue_otsu_thresh)
        ) +
        coord_sf() +
        theme_minimal() +
        labs(
            title = sprintf("In tissue bins (> Otsu threshold %.3f)", in_tissue_otsu_thresh),
            subtitle = sprintf("Total %d bins above threshold", sum(grid_sf$in_tissue))
        )
    
    p4 <- ggplot(grid_sf) +
        geom_sf(aes(fill = MAPK3_counts), color = "grey70", size = 0.1) +
        scale_fill_viridis_c(option="plasma") +
        coord_sf() +
        theme_minimal() +
        labs(
            title = sprintf("Raw MAPK3 counts", in_tissue_otsu_thresh))
        
    print(p + p2 + p3 + p4)
    rm(transcripts)
    rm(grid_sf)
    gc()
dev.off()
