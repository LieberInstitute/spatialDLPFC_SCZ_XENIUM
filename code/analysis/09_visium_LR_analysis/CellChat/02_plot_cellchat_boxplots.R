library(dplyr)
library(here)
library(ggplot2)
library(patchwork)
library(tidyverse)


#####################################################################
# Read in the CellChat results csvs and make boxplots of 
# number of LR pairs for NTC and SCZ within each microenvironment
#####################################################################


files <- list.files(here(
    "processed-data", "09_visium_LR_analysis", "CellChat"
), full.names = TRUE)

microenvs <- c("neuropil", "neun", "pnn", "vasc")
p_list <- list()
for (j in 1:length(microenvs)){
    me_use <- microenvs[[j]]

    files_use <- files[grepl(me_use, files)]
    df_list <- lapply(files_use, function(x) {
        df <- read.csv(x)
        df$microenv <- me_use

        print(df)

        df_within <- df %>%
            group_by(Dx, microenv, BrNum) %>%
            summarize(
                n_LR_pairs = sum(pval < 0.05 & source & target, na.rm = TRUE),
                .groups = "drop"
        )
        df_out <- df %>%
            group_by(Dx, microenv, BrNum) %>%
            summarize(
                n_LR_pairs = sum(pval < 0.05 & !(source & target) & !(!source & !target), na.rm = TRUE),
                .groups = "drop"
        )
         
        df_within$interaction_type <- "within_microenv"
        df_out$interaction_type <- "outside_microenv"
        df <- rbind(df_within, df_out)  
         
        return(df)


    })

    df <- do.call(rbind, df_list)
    p_list[[j]] <- ggplot(df, aes(x=Dx, y=n_LR_pairs))+
        geom_boxplot()+
        ggtitle(me_use)+
        xlab("Diagnosis")+
        ylab("Number of significant LR pairs")+
        facet_wrap(~interaction_type)

}


pdf(here("plots", "09_visium_LR_analysis", "CellChat", "cellchat_boxplots_microenvs.pdf"),
    width=8, height=6)

patchwork::wrap_plots(p_list, ncol=2)
dev.off()