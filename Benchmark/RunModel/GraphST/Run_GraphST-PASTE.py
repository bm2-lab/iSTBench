import os
import argparse
import torch
import pandas as pd
import scanpy as sc
import anndata as ad
from sklearn import metrics
import multiprocessing as mp
import matplotlib.pyplot as plt
from GraphST import GraphST
import numpy as np
from GraphST.utils import clustering

import sys
import math
import scipy
import seaborn as sns
from matplotlib import style
import matplotlib
import sklearn
import networkx as nx
import ot
import paste as pst
import plotly.express as px
import plotly.io as pio



def process_data(input_file, output_file, sample, nclust, hvgs, run_normalization, tool,device, seed):
    np.random.seed(seed)
    #os.environ['R_HOME'] = '/usr/lib/R'

    adata = sc.read_h5ad(input_file)
    adata.var_names_make_unique()
    adata.obs['data'] = adata.obs['slices'].astype(str)
    count_data = adata.X.copy()

    if run_normalization:
        total_counts_per_cell = np.sum(adata.X, axis=1) 

        sc.pp.normalize_total(adata, target_sum=total_counts_per_cell)
        sc.pp.log1p(adata)
        sc.pp.scale(adata, max_value=10)
        sc.pp.highly_variable_genes(adata, n_top_genes=hvgs)
        
        highly_variable_data = adata[:, adata.var['highly_variable']]
        highly_variable_data.obs_names = highly_variable_data.obs_names + '_' + highly_variable_data.obs['slices'].astype(str)
    else:
        highly_variable_data = adata
        highly_variable_data.obs_names = highly_variable_data.obs_names + '_' + highly_variable_data.obs['slices'].astype(str)

    
    slices_values = highly_variable_data.obs['slices'].unique()
    adata_list = {}
    for slice_value in slices_values:
      adata_slice = highly_variable_data[highly_variable_data.obs['slices'] == slice_value]
      adata_list[slice_value] =  adata_slice
      
    sample_groups = [slices_values]
    layer_groups = [[adata_list[sample] for sample in slices_values]]
    
    alpha = 0.1
    res_df = pd.DataFrame(columns=['Sample','Pair','Kind'])
    pis = [[None for i in range(len(layer_groups[j])-1)] for j in range(len(layer_groups))]

    for j in range(len(layer_groups)):
       for i in range(len(layer_groups[j])-1):
           pi0 = pst.match_spots_using_spatial_heuristic(layer_groups[j][i].obsm['spatial'],layer_groups[j][i+1].obsm['spatial'],use_ot=True)
           pis[j][i] = pst.pairwise_align(layer_groups[j][i], layer_groups[j][i+1],alpha=alpha,G_init=pi0,norm=True,verbose=False)
           res_df.loc[len(res_df)] = [j,i,'PASTE']
    
    paste_layer_groups = [pst.stack_slices_pairwise(layer_groups[j], pis[j]) for j in range(len(layer_groups))]
    combined_adata = ad.concat(paste_layer_groups[0])
    combined_adata.var = highly_variable_data.var
    combined_adata.uns = highly_variable_data.uns

    torch.set_num_threads(10)
    model = GraphST.GraphST(combined_adata, device=device, random_seed = seed)
    adata = model.train()

    if tool == 'mclust':
        clustering(adata, nclust, method=tool, refinement=True)
    elif tool in ['leiden', 'louvain']:
        clustering(adata, nclust, method=tool, start=0.1, end=2.0, increment=0.01)
    
    adata.obsm["X_embedding"] = adata.obsm["emb"]
    adata.obs['predicted_domain'] = adata.obs['domain']
    
    adata.write(f"{output_file}/{sample}.h5ad")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')

    parser.add_argument('--input_file', type=str, help='input file path')
    parser.add_argument('--output_file', type=str, help='output file path')
    parser.add_argument('--sample', type=str, help='sample name')
    parser.add_argument('--nclust', type=int, help='number of clusters')
    parser.add_argument('--hvgs',type=int, default=3000, help='number of highly variable genes')
    parser.add_argument('--runNormalization', type=str, default='True', help='whether to run normalization')
    parser.add_argument('--tool', type=str, default = "mclust", help='clustering tool')
    parser.add_argument('--device', type=str, default = "cpu", help='gpu or cpu')
    parser.add_argument('--seed', type=int, default = 41, help='random seed')

    args = parser.parse_args()

    run_normalization = args.runNormalization.lower() in ('true', 't', 'yes', 'y', '1')

    process_data(args.input_file, args.output_file, args.sample, args.nclust, args.hvgs, run_normalization, args.tool, args.device, args.seed)
