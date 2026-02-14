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
  library(data.table)
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



pdf(here("plots", "10_eQTL",sprintf("02_%s_neuropil_MAPK3.pdf", brnums[[sample_idx]])),
    width = 21, height = 6)

    print(brnums[[sample_idx]])
    transcripts <- arrow::read_parquet(transcript_files[sample_idx])
    transcripts <- st_as_sf(transcripts, 
        coords = c("x_location", "y_location"), crs=NA, remove=FALSE)
    
    # Create square bins
    sfc <- st_make_grid(st_bbox(transcripts), square=TRUE, cellsize = 75)
    grid_sf <- st_sf(geometry = sfc)
    grid_sf$grid_id <- seq_len(nrow(grid_sf))
    

    # Calculate mean MOBP counts in each bin

    transcripts <- transcripts %>%
        filter(is_gene & qv > 20)
    
    tx_dt <- as.data.table(transcripts)
    cellsize <- 75
    tx_dt[, `:=`(
        xbin = floor((x_location - min(x_location)) / cellsize) + 1,
        ybin = floor((y_location - min(y_location)) / cellsize) + 1)]


    print("tx_dt made")
    Sys.time()

    grid_dt <- as.data.table(grid_sf)

    grid_centroids <- st_centroid(grid_sf$geometry)
    grid_xy <- st_coordinates(grid_centroids)

    grid_dt[, `:=`(
        xbin = floor((grid_xy[,1] - min(tx_dt$x_location)) / cellsize) + 1,
        ybin = floor((grid_xy[,2] - min(tx_dt$y_location)) / cellsize) + 1
    )]
    
    print("grid_dt made")
    Sys.time()
    grid_dt[, grid_id := seq_len(.N)]
    tx_dt <- merge(tx_dt, grid_dt[, .(xbin, ybin, grid_id)], by = c("xbin", "ybin"), all.x = TRUE)


    print("merge completed")
    Sys.time()
    total_counts <- tx_dt[!is.na(grid_id), 
                      .(total_counts = .N), 
                      by = grid_id]

    mobp_counts <- tx_dt[feature_name == "MOBP" & qv > 20 & !is.na(grid_id), 
                     .(MOBP_counts = .N), 
                     by = grid_id]

    extrasomatic_counts <- tx_dt[cell_id == "UNASSIGNED" & qv > 20 & 
                             !grepl("NegControl|Unassigned", feature_name) & 
                             !is.na(grid_id),
                             .(n_transcripts = .N),
                             by = .(grid_id, feature_name)]
    gene_counts_wide <- dcast(extrasomatic_counts, 
                          grid_id ~ feature_name, 
                          value.var = "n_transcripts", 
                          fill = 0)

    grid_dt <- merge(grid_dt, total_counts, by = "grid_id", all.x = TRUE)
    grid_dt <- merge(grid_dt, mobp_counts,  by = "grid_id", all.x = TRUE)
    grid_dt <- merge(grid_dt, gene_counts_wide, by = "grid_id", all.x = TRUE)

    # Fill NAs with 0
    cols_to_fill <- c("total_counts","MOBP_counts", setdiff(names(grid_dt), names(grid_sf)))
    for (col in cols_to_fill) grid_dt[is.na(get(col)), (col) := 0]

    grid_dt <- grid_dt %>%
        mutate(MOBP_norm = ifelse(total_counts > 0, MOBP_counts / total_counts, 0))

    print(summary(grid_sf$MOBP_norm))

    # --- 4. Otsu thresholding to identify high MOBP bins ---
    otsu_thresh <- pliman::otsu(grid_dt$MOBP_norm)
    print(otsu_thresh)
    grid_dt[, passes_otsu := MOBP_norm > otsu_thresh]

    # Also compute Otsu threshold for in-tissue bins based on total counts
    in_tissue_otsu_thresh <- pliman::otsu(grid_dt$total_counts)
    print(in_tissue_otsu_thresh)
    grid_dt[, in_tissue := total_counts > in_tissue_otsu_thresh]

    
    # Compute bin coordinates from xbin / ybin
    cellsize <- 75
    grid_dt[, `:=`(
        xmin = (xbin - 1) * cellsize + min(tx_dt$x_location),
        xmax = xbin * cellsize + min(tx_dt$x_location),
        ymin = (ybin - 1) * cellsize + min(tx_dt$y_location),
        ymax = ybin * cellsize + min(tx_dt$y_location)
        )]

    print(head(grid_dt))

    p <- ggplot(grid_dt) +
        geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = MOBP_norm),
                    color = "grey70", size = 0.1) +
        scale_fill_viridis_c(option = "plasma") +
        coord_fixed() +
        theme_minimal() +
        labs(title = "MOBP normalized counts per 75x75 bin",
            fill = "MOBP_norm")
        

    # keep only grey matter bins
    grid_dt_grey <- grid_dt %>%
        filter(!passes_otsu) 

    
    print(head(grid_dt_grey))


    saveRDS(grid_dt_grey, here("processed-data", "10_eQTL", "neuropil",
        sprintf("%s-neuropil-mapk3-counts_grey_matter.RDS", brnums[[sample_idx]])))

    p2 <- ggplot(grid_dt) +
         geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = passes_otsu),
                    color = "grey70", size = 0.1) +
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

    p3 <- ggplot(grid_dt) +
         geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = in_tissue),
                    color = "grey70", size = 0.1) + 
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
    
    p4 <- ggplot(grid_dt) +
         geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = MAPK3),
                    color = "grey70", size = 0.1) +
        scale_fill_viridis_c(option="plasma") +
        coord_sf() +
        theme_minimal() +
        labs(title = "Raw MAPK3 counts")
        
    print(p + p2 + p3 + p4)
    rm(transcripts)
    rm(grid_sf)
    gc()
dev.off()
