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

paths <- lapply(raw_paths, function(x) {
  list.files(x, full.names = TRUE)
})

paths <- lapply(strsplit(split="//", unlist(paths)), "[", -1)
paths <- unlist(paths)

raw_path_with_typo <- "/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2024-11-12_Xenium-pnn-SKCCC/"
paths_with_typo <- list.files(raw_path_with_typo, full.names = TRUE)
paths_with_typo <- lapply(strsplit(split="//", unlist(paths_with_typo)), "[", -1)
paths_with_typo <- unlist(paths_with_typo)

paths <- c(brnums, brnums_with_typo)

# Extract the Br**** part from each brnum
brnums <- unlist(lapply(strsplit(brnums, split="__"), "[", 3))
