# Import necessary libraries and modules
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

# Function to process the input data, normalize if specified, perform PASTE, and run clustering
# This function reads the input data, applies normalization (if specified), selects highly variable genes, and performs the PASTE spatial alignment
# After that, it trains the GraphST model and performs clustering using the selected method.
def process_data(input_file, output_file, sample, nclust, hvgs, run_normalization, tool, device, seed):
    # Parameters:
    # input_file: Path to the multi-slice data that needs to be integrated. All slices should be contained in a single h5ad file. Refer to the example data for the specific format. This is a required parameter.
    # output_file: Path to the directory where the output results will be saved. This is a required parameter.
    # sample: Name for the output result files. It is recommended to use the model name, for example, setting it as "GraphST" will result in output files named "GraphST.h5ad". This is a required parameter.
    # nclust: Number of identified domains. This is a required parameter.
    # hvgs: The number of highly variable genes identified if data normalization is performed. This is an optional parameter, with the default set to 3000.
    # runNormalization: A boolean variable indicating whether data normalization should be performed before integration. This is an optional parameter, with the default set to True.
    # tool: Methods for spatial clustering, which can be set to "mclust", "leiden" and "louvain". This is an optional parameter, with the default set to "mclust".
    # device: The device on which the model runs can be set to "cpu" and "gpu". This is an optional parameter, with the default set to "cpu".
    # seed: Random seed for model execution. This is an optional parameter, with the default set to 0.

    np.random.seed(seed)  # Set the random seed for reproducibility
    
    # Set the environment variable for R (required for Scanpy)
    # Set R based on your environment
    # os.environ['R_HOME'] = '/opt/R/4.3.2/lib/R'

    # Read input data and prepare AnnData object
    adata = sc.read_h5ad(input_file)
    adata.var_names_make_unique()  
    adata.obs['data'] = adata.obs['slices'].astype(str)  
    count_data = adata.X.copy() 

    # Perform normalization if specified
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

    # Prepare data for PASTE: Split the data into slices and align them
    slices_values = highly_variable_data.obs['slices'].unique() 
    adata_list = {} 
    for slice_value in slices_values:
        adata_slice = highly_variable_data[highly_variable_data.obs['slices'] == slice_value]
        adata_list[slice_value] = adata_slice  

    sample_groups = [slices_values]  
    layer_groups = [[adata_list[sample] for sample in slices_values]] 

    # Initialize variables for the PASTE alignment process
    alpha = 0.1  # Regularization parameter for PASTE
    res_df = pd.DataFrame(columns=['Sample', 'Pair', 'Kind']) 
    pis = [[None for i in range(len(layer_groups[j]) - 1)] for j in range(len(layer_groups))]  # List for storing pairwise alignments

    # Perform pairwise alignment using PASTE
    for j in range(len(layer_groups)):
        for i in range(len(layer_groups[j]) - 1):
            pi0 = pst.match_spots_using_spatial_heuristic(layer_groups[j][i].obsm['spatial'], layer_groups[j][i + 1].obsm['spatial'], use_ot=True)  
            pis[j][i] = pst.pairwise_align(layer_groups[j][i], layer_groups[j][i + 1], alpha=alpha, G_init=pi0, norm=True, verbose=False)  
            res_df.loc[len(res_df)] = [j, i, 'PASTE']  

    # Stack the aligned slices and concatenate them into a single AnnData object
    paste_layer_groups = [pst.stack_slices_pairwise(layer_groups[j], pis[j]) for j in range(len(layer_groups))]
    combined_adata = ad.concat(paste_layer_groups[0]) 
    combined_adata.var = highly_variable_data.var 
    combined_adata.uns = highly_variable_data.uns 

    # Set the number of threads for PyTorch
    torch.set_num_threads(10)

    # Initialize and train the GraphST model
    model = GraphST.GraphST(combined_adata, device=device, random_seed=seed)
    adata = model.train()

    # Perform clustering using the specified tool
    if tool == 'mclust':
        clustering(adata, nclust, method=tool, refinement=True)  # Use MCLUST for clustering
    elif tool in ['leiden', 'louvain']:
        clustering(adata, nclust, method=tool, start=0.1, end=2.0, increment=0.01)  # Use Leiden or Louvain for clustering

    # Add embeddings and predicted domains to the AnnData object
    adata.obsm["X_embedding"] = adata.obsm["emb"]
    adata.obs['predicted_domain'] = adata.obs['domain']

    # Save the results to an h5ad file
    adata.write(f"{output_file}/{sample}.h5ad")

# Main function to parse arguments and run the data processing pipeline
# This function handles command-line input, sets up the necessary arguments, and calls the process_data function
if __name__ == "__main__":
    # Set up the argument parser for command-line inputs
    parser = argparse.ArgumentParser(description='Process some integers.')

    # Define arguments for input file, output file, sample name, etc.
    parser.add_argument('--input_file', type=str, help='input file path')
    parser.add_argument('--output_file', type=str, help='output file path')
    parser.add_argument('--sample', type=str, help='sample name')
    parser.add_argument('--nclust', type=int, help='number of clusters')
    parser.add_argument('--hvgs', type=int, default=3000, help='number of highly variable genes')
    parser.add_argument('--runNormalization', type=str, default='True', help='whether to run normalization')
    parser.add_argument('--tool', type=str, default="mclust", help='clustering tool')
    parser.add_argument('--device', type=str, default="cpu", help='gpu or cpu')
    parser.add_argument('--seed', type=int, default=41, help='random seed')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Convert normalization flag to boolean
    run_normalization = args.runNormalization.lower() in ('true', 't', 'yes', 'y', '1')

    # Call the process_data function with the parsed arguments
    process_data(args.input_file, args.output_file, args.sample, args.nclust, args.hvgs, run_normalization, args.tool, args.device, args.seed)
