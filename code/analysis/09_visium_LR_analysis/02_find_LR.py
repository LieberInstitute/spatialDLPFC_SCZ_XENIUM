# use the ~/mambaforge/envs/scSLAT/ conda environment
import os
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc

import spatialdm as sdm
from spatialdm.datasets import dataset
import spatialdm.plottings as pl

from spatialdm.diff_utils import *

# tutorial: https://spatialdm.readthedocs.io/en/latest/api.html
# https://spatialdm.readthedocs.io/en/latest/differential_test_intestine.html
# load the neuronal data
adata = sc.read_h5ad("/dcs05/lieber/marmaypag/spatialDLPFC_SCZ_XENIUM_LIBD4100/spatialDLPFC_SCZ_XENIUM/processed-data/09_visium_LR_analysis/sce_neurons.h5ad")

# split each sample into its own anndata object
samples = adata.obs['BrNumbr'].unique().tolist()
adatas = []
for sample in samples:
    adatas.append(adata[adata.obs['BrNumbr'] == sample].copy())
    
# run spatialDM for each sample
for adata in adatas:
    # get spatial coords from coord_x and coord_y in coldata
    adata.obsm['spatial'] = adata.obs[['coord_x', 'coord_y']].to_numpy()
    adata.X = adata.layers['logcounts']
    
    try:
        sdm.weight_matrix(adata, l=110, cutoff=0.2, single_cell=False) # weight_matrix by rbf kernel
        sdm.extract_lr(adata, 'human', min_cell=3)      # find overlapping LRs from CellChatDB
        sdm.spatialdm_global(adata, 1000, specified_ind=None, method='z-score', nproc=1)     # global Moran selection
        sdm.sig_pairs(adata, method='z-score', fdr=True, threshold=0.1)     # select significant pairs
        sdm.spatialdm_local(adata, n_perm=1000, method='z-score', specified_ind=None, nproc=1)     # local spot selection
        sdm.sig_spots(adata, method='z-score', fdr=False, threshold=0.1)
        
        print(adata)
    except:
        print("no LR pairs found for sample" + adata.obs['BrNumbr'][0])
        continue

concat=concat_obj(adatas, samples, 'human', 'z-score', fdr=False)
print(concat)