library(here)

#####################################################################################
# From the paths of the raw data files, get the brnums 
# for each of the Xenium samples and save them into a file.
# Locations:
#/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2024-11-12_Xenium-pnn-SKCCC/
#/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2024-11-18_xenium-pnn-SKCCC/
#/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2024-12-06_xenium-pnn-SKCCC/
#/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2024-12-12_xenium-pnn-SKCCC/
#####################################################################################

raw_path_template <- "/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/%s_xenium-pnn-SKCCC/"
run_dates <- c("2024-12-12", "2024-12-06", "2024-11-18", "2024-11-12")
raw_paths <- lapply(raw_path_template, sprintf, run_dates)

full_paths <- lapply(raw_paths, function(x) {
  list.files(x, full.names = TRUE)
})

paths <- lapply(strsplit(split="//", unlist(full_paths)), "[", -1)
paths <- unlist(paths)

raw_path_with_typo <- "/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2024-11-12_Xenium-pnn-SKCCC/"
full_paths_with_typo <- list.files(raw_path_with_typo, full.names = TRUE)
paths_with_typo <- lapply(strsplit(split="//", unlist(full_paths_with_typo)), "[", -1)
paths_with_typo <- unlist(paths_with_typo)

paths <- c(paths, paths_with_typo)

full_paths <- c(unlist(full_paths), unlist(full_paths_with_typo))

# Extract the Br**** part from each brnum
brnums <- unlist(lapply(strsplit(paths, split="__"), "[", 3))

sample_info_df <- data.frame(cbind(BrNum=brnums, path=full_paths))
write.csv(sample_info_df, here("raw-data", "experiment_info", "sample_paths.csv"))
