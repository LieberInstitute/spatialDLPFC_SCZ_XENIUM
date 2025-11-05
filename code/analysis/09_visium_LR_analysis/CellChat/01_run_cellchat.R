library(scater)
library(CellChat)
library(Seurat)
library(dplyr)
library(here)


#####################################################################
# Run CellChat across each of the samples, using
# spots from all microenvironments.
#####################################################################


source(here("code", "analysis", "09_visium_LR_analysis", "CellChat", "modeling.R"))
args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[[1]])


spe <- readRDS(
  here::here(
    "/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/",
    "processed-data/rds/01_build_spe",
    "fnl_spe_kept_spots_only.rds"
  )
)

# Check if any spots belong to multiple microenvironments
table(spe$microenv)
table(rowSums(as.data.frame(colData(spe))[,c("neuropil_pos", "neun_pos", "pnn_pos", "vasc_pos")]))
#      0      1      2      3      4 
#  74195 175675  28283   1632     21 


brnums <- unique(spe$BrNumbr)

br_use <- brnums[[i]]
spe_use <- spe[,spe$BrNumbr == br_use]
print(sprintf("Running CellChat on sample %s", br_use))

data.input <- logcounts(spe_use)
rownames(data.input) <- rowData(spe_use)$gene_name # Convert emsembl id to gene sym


# loop over samples and run cellchat for each sample
microenvs <- c("neuropil", "neun", "pnn", "vasc")
for (j in 1:length(microenvs)){

    me_use <- microenvs[[j]]
    print(sprintf("  Microenvironment: %s", me_use))

    # create object for cellchat
    colData(spe_use)$is_microenv <- colData(spe_use)[,grepl(me_use, colnames(colData(spe)))]

    meta <- as.data.frame(SingleCellExperiment::colData(spe_use))
    meta$samples <- meta$sample_id |> factor()
    cellchat <- createCellChat(
      object = data.input, meta = meta,
      group.by = "is_microenv")

    unique(cellchat@idents)

    # add metadata to cellchat object
    cellchat <- addMeta(cellchat, meta = meta)
    cellchat <- setIdent(cellchat, ident.use = "is_microenv") 
    groupSize <- as.numeric(table(cellchat@idents))

    # Attach database
    cellchat@DB <- CellChatDB.human

    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat, do.DE = FALSE)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = 5)
    
    comm_df <- subsetCommunication(cellchat)

    print(head(comm_df))


    write.csv(comm_df,
      file = here::here(
        "results",
        "09_visium_LR_analysis",
        "CellChat",
        sprintf("cellchat_communications_BrNumbr_%s_microenv_%s.csv",
          br_use,
          me_use)
      ),
      row.names = FALSE
    )


}

    

