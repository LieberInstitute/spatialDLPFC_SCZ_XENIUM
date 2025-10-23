library(here)
library(tidyverse)

# read in the data
me <- "neuropil"
files <- list.files(here("processed-data", "09_visium_LR_analysis", me), full.names=TRUE)

dfs <- lapply(files, read.csv, sep="\t")

brnums <- strsplit(files, "global_res_")
brnums <- unlist(lapply(brnums, "[", 2))

for(i in 1:length(dfs)){
    df <- dfs[[i]]
    df$BrNum <- strsplit(brnums[[i]], "_",)[[1]][[1]]
    df$Dx <- strsplit(brnums[[i]], "_",)[[1]][[2]]
    dfs[[i]] <- df
}
res <- do.call(rbind, dfs)

res %>% group_by(BrNum) %>%
    summarise(n_selected = sum(selected=="True"))

scz_lrs <- res %>% filter(Dx=="scz") %>%
    group_by(X) %>%
    summarise(n_selected = sum(selected=="True"))%>%
    arrange(desc(n_selected))


ntc_lrs <- res %>% filter(Dx=="ntc") %>%
    group_by(X) %>%
    summarise(n_selected = sum(selected=="True"))%>%
    arrange(desc(n_selected))

all_lrs <- res %>% group_by(X) %>%
    summarise(n_selected = sum(selected=="True"))%>%
    arrange(desc(n_selected)) %>%
    filter(n_selected >0) %>%
    filter(n_selected == 63)

all_scz_lrs <- res %>% group_by(X) %>%
    filter(Dx=="SCZ") %>%
    summarise(n_selected = sum(selected=="True"))%>%
    arrange(desc(n_selected)) 

all_lrs <- all_lrs %>%
    separate(X, into=c("Ligand", "Receptor1", "Receptor2"), sep="_")

print(all_lrs)

# read in spg-degs from the visium analysis:
degs <- read.csv("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/code/analysis/dx_deg_spg_neun/neun-dx_DEG-GM.csv")

degs <- degs %>%
    filter(p_value_scz < 0.05)

all_lrs$Ligand %in% degs
all_lrs$Receptor1 %in% degs
all_lrs$Receptor1 %in% degs