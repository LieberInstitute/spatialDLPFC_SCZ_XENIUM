# use the ~/mambaforge/envs/scSLAT/ conda environment
import os
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import scipy
import seaborn as sns
import matplotlib.pyplot as plt


import spatialdm as sdm
from spatialdm.datasets import dataset
import spatialdm.plottings as pl
from scipy.stats.mstats import gmean
from argparse import ArgumentParser
from spatialdm.diff_utils import *

# tutorial: https://spatialdm.readthedocs.io/en/latest/api.html
# https://spatialdm.readthedocs.io/en/latest/differential_test_intestine.html
# load the neuronal data

parser = ArgumentParser()
parser.add_argument("-me", "--microenv", type=str, help="which microenvironment to run spatialDM on")
args = parser.parse_args()
me = args.microenv

print(me, flush=True)
# microenvs = ["neurons", "neuropil", "pnn", "vasc"]
#for me in microenvs:
print(me, flush=True)
adata = sc.read_h5ad("/dcs05/lieber/marmaypag/spatialDLPFC_SCZ_XENIUM_LIBD4100/spatialDLPFC_SCZ_XENIUM/processed-data/09_visium_LR_analysis/sce_" + me + ".h5ad")
adata.obs["BrNumbr"] = adata.obs["BrNumbr"].astype("str") + "_" + adata.obs["dx"].astype("str")
res_fname = "/dcs05/lieber/marmaypag/spatialDLPFC_SCZ_XENIUM_LIBD4100/spatialDLPFC_SCZ_XENIUM/processed-data/09_visium_LR_analysis/sce_" + me + "_spatialDM.h5ad"

# Order the samples by Dx
adata.obs["dx"] = pd.Categorical(
    adata.obs["dx"], 
    categories=["scz", "ntc"], 
    ordered=True
)
adata = adata[adata.obs.sort_values("dx").index, :]

# split each sample into its own anndata object
samples = adata.obs['BrNumbr'].unique().tolist()
#samples = samples[1:4] # test a few samples, remove this line later!!
adatas = []
for sample in samples:
    adatas.append(adata[adata.obs['BrNumbr'] == sample].copy())
    
    
def differential_volcano(sample, pairs=None, legend=None, xmax = 25, xmin = -20):
    """
    Volcano plot for a differential obj
    :param sample: concatenated AnnData after running spatialdm separately
    :param pairs: list
    list of pairs to be highlighted in the volcano plot, e.g. ['SPP1_CD44'] or ['SPP1_CD44','ANGPTL4_SDC2']
    :param legend: list
    list of specified names for each side of the volcano plot
    :param xmax: float
    max z-score difference
    :param xmin: float
    min z-score difference
    :return: ax: matplotlib Axes.
    """
    if pairs is not None:
        color_codes = generate_colormap(max(10, len(pairs)+8))[8:]
    q1 = sample.uns['q1']
    q2 = sample.uns['q2']
    fdr_co = sample.uns['fdr_co']
    _range = np.arange(1, sample.uns['n_sub'])
    diff_cp = sample.uns['diff'].copy()
    diff_cp = np.where((diff_cp>xmax), xmax, diff_cp)
    diff_cp = np.where((diff_cp<xmin), xmin, diff_cp)
    plt.scatter(diff_cp[sample.uns['tf_df'].sum(1).isin(_range)],
                -np.log10(sample.uns['diff_fdr'])[sample.uns['tf_df'].sum(1).isin(_range)], s=10, c='grey')
    keys = sample.uns.keys()
    conditions = []
    for key in keys:
        if key.endswith('_specific'):
            conditions.append(key.replace('_specific',''))
    label = 'difference between z-score of {0[0]} and {0[1]}'.format(conditions)
    plt.xlabel(label)
    plt.ylabel('differential fdr (log-likelihood, -log10)')
    plt.xlim([xmin-1,xmax+1])
    plt.scatter(diff_cp[(diff_cp>q1) & (sample.uns['diff_fdr']<fdr_co) & \
                           (sample.uns['tf_df'].sum(1).isin(_range))],
                -np.log10(sample.uns['diff_fdr'])[(diff_cp>q1) & (sample.uns['diff_fdr']<fdr_co) & \
                           (sample.uns['tf_df'].sum(1).isin(_range))], s=10,c='tab:orange')
    plt.scatter(diff_cp[(diff_cp<q2) & (sample.uns['diff_fdr']<fdr_co) & \
                           (sample.uns['tf_df'].sum(1).isin(_range))],
                -np.log10(sample.uns['diff_fdr'])[(diff_cp<q2) & (sample.uns['diff_fdr']<fdr_co)& \
                           (sample.uns['tf_df'].sum(1).isin(_range))], s=10,c='tab:green')
    if type(pairs)!=type(None):
        for i,pair in enumerate(pairs):
            plt.scatter(diff_cp[sample.uns['p_df'].index==pair],
                        -np.log10(sample.uns['diff_fdr'])[sample.uns['p_df'].index==pair], c=color_codes[i])
    plt.legend(np.hstack(([''], legend, pairs)))
    return(plt)
    
# run spatialDM for each sample
if os.path.exists(res_fname):
    print("file exists, reading it in", flush=True)
    concat = sc.read_h5ad(res_fname)
else:
    print("file does not exist, running spatialDM")
    for adata in adatas:
        # get spatial coords from coord_x and coord_y in coldata
        adata.obsm['spatial'] = adata.obs[['coord_x', 'coord_y']].to_numpy()
        adata.X = adata.layers['logcounts'].toarray()
        adata.var_names_make_unique()
        sdm.weight_matrix(adata, l=75, cutoff=0.2, single_cell=False) # weight_matrix by rbf kernel
        sdm.extract_lr(adata, 'human', min_cell=3, mean="algebra") 
        sdm.spatialdm_global(adata, 1000, specified_ind=None, method='z-score', nproc=1)  
        sdm.sig_pairs(adata, method='z-score', fdr=True, threshold=0.1)
        #print(adata.uns['global_res'][adata.uns['global_res']['selected'] == True], flush=True)
        # sdm.spatialdm_global(adata, 1000, specified_ind=None, method='z-score', nproc=1)     # global Moran selection
        # sdm.sig_pairs(adata, method='z-score', fdr=True, threshold=0.1)     # select significant pairs
        # sdm.spatialdm_local(adata, n_perm=1000, method='z-score', specified_ind=None, nproc=1)     # local spot selection
        # sdm.sig_spots(adata, method='z-score', fdr=False, threshold=0.1)


    concat=concat_obj(adatas, samples, 'human', 'z-score', fdr=False)
    print(concat, flush=True)

    
    print(concat, flush=True)
    print(concat.uns['p_df'].head(), flush=True)
    print(concat.uns['tf_df'].head(), flush=True)
    print(concat.uns['zscore_df'].head(), flush=True)
    
    


    fig=sns.clustermap(1-concat.uns['p_df'])
    fig.savefig("/dcs05/lieber/marmaypag/spatialDLPFC_SCZ_XENIUM_LIBD4100/spatialDLPFC_SCZ_XENIUM/plots/09_visium_LR_analysis/spatialDM_clustermap_" + me + ".png", dpi=300)

    # DE analysis

    
    dx = np.array(list(map(lambda item: 1 if "scz" in item else 0, samples)))
    print(dx, flush=True)
    print(samples, flush=True)
    print(concat.uns['zscore_df'].loc[:, np.array(samples)[dx == 1]].mean(1), flush=True)
    differential_test(cdata=concat, subset=samples, conditions=dx)
    group_differential_pairs(concat, 'scz', 'ntc')
    
    #fig=pl.differential_dendrogram(concat)
    _range = np.arange(1, concat.uns['n_sub'])
    fig = sns.clustermap(1-concat.uns['p_df'].loc[(concat.uns['p_val']<0.1) & (concat.uns['tf_df'].sum(1).isin(_range)),
                                     concat.uns['subset']], col_cluster=False)
    fig.savefig("/dcs05/lieber/marmaypag/spatialDLPFC_SCZ_XENIUM_LIBD4100/spatialDLPFC_SCZ_XENIUM/plots/09_visium_LR_analysis/spatialDM_dendrogram_" + me + ".png", dpi=300)
    
    volcano=differential_volcano(concat, legend=['scz specific', 'ntc specific'])
    volcano.savefig("/dcs05/lieber/marmaypag/spatialDLPFC_SCZ_XENIUM_LIBD4100/spatialDLPFC_SCZ_XENIUM/plots/09_visium_LR_analysis/spatialDM_volcano_" + me + ".png", dpi=300)
    
    # concat.uns["ligand"] = concat.uns["ligand"].astype(str)
    # concat.uns["receptor"] = concat.uns["receptor"].astype(str)
    # concat.uns["geneInter"] = concat.uns["geneInter"].astype(str)
    # concat.uns["p_df"] = concat.uns["p_df"].astype(str)
    # concat.uns["tf_df"] = concat.uns["tf_df"].astype(str)
    # concat.uns["zscore_df"] = concat.uns["zscore_df"].astype(str)
    # concat.uns["geneInter"] = concat.uns["geneInter"].to_numpy()
    # concat.uns["ligand"] = concat.uns["ligand"].to_numpy()
    # concat.uns["receptor"] = concat.uns["receptor"].to_numpy()
    # concat.uns["p_df"] = concat.uns["p_df"].to_numpy()
    # concat.uns["tf_df"] = concat.uns["tf_df"].to_numpy()
    # concat.uns["zscore_df"] = concat.uns["zscore_df"].to_numpy()
    
    concat.uns['geneInter'] = concat.uns['geneInter'].fillna('NA')
    concat.uns['ligand'] = concat.uns['ligand'].fillna('NA')
    concat.uns['receptor'] = concat.uns['receptor'].fillna('NA')
    concat.uns['p_df'] = concat.uns['p_df'].fillna('NA')
    concat.uns['tf_df'] = concat.uns['tf_df'].fillna('NA')
    concat.uns['zscore_df'] = concat.uns['zscore_df'].fillna('NA')
    
    fdr_df = pd.DataFrame(concat.uns['diff_fdr'])
    fdr_df.to_csv("/dcs05/lieber/marmaypag/spatialDLPFC_SCZ_XENIUM_LIBD4100/spatialDLPFC_SCZ_XENIUM/processed-data/09_visium_LR_analysis/diff_fdr_" + me + ".csv", sep="\t")
    
    tf_df = pd.DataFrame(concat.uns['tf_df'])
    tf_df.to_csv("/dcs05/lieber/marmaypag/spatialDLPFC_SCZ_XENIUM_LIBD4100/spatialDLPFC_SCZ_XENIUM/processed-data/09_visium_LR_analysis/tf_df_" + me + ".csv", sep="\t")
    
    diff_df = pd.DataFrame(concat.uns['diff'])
    diff_df.to_csv("/dcs05/lieber/marmaypag/spatialDLPFC_SCZ_XENIUM_LIBD4100/spatialDLPFC_SCZ_XENIUM/processed-data/09_visium_LR_analysis/diff_" + me + ".csv", sep="\t")
    
    
    
    print(concat, flush=True)
    concat.write(res_fname)
    #sdm.write_spatialdm_h5ad(concat, res_fname)
    
    #concat.write_h5ad("/dcs05/lieber/marmaypag/spatialDLPFC_SCZ_XENIUM_LIBD4100/spatialDLPFC_SCZ_XENIUM/processed-data/09_visium_LR_analysis/sce_" + me + "_spatialDM.h5ad")
    