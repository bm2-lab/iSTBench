import os
import argparse
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp
import matplotlib.pyplot as plt
from GraphST import GraphST
import numpy as np
from GraphST.utils import clustering
import time
from memory_profiler import memory_usage
import csv

def run_model(highly_variable_data, device, nclust, tool, seed):
    np.random.seed(seed)
    def model_execution():
        model = GraphST.GraphST(highly_variable_data, device=device, random_seed = seed)
        adata = model.train()

        if tool == 'mclust':
            clustering(adata, nclust, method=tool, refinement=True)
        elif tool in ['leiden', 'louvain']:
            clustering(adata, nclust, method=tool, start=0.1, end=2.0, increment=0.01)

        return adata
      
    start_time = time.time()
    mem_usage = memory_usage((model_execution,), max_usage = True, retval=True)
    end_time = time.time()

    adata = mem_usage[1]
    mem_usage_values = mem_usage[0] 
    return mem_usage_values, end_time - start_time, adata

def process_data(input_file, output_file, sample, nclust, hvgs, run_normalization, tool, device, seed):
    #os.environ['R_HOME'] = '/usr/lib/R'

    adata = sc.read_h5ad(input_file)
    adata.var_names_make_unique()
    adata.obs['data'] = adata.obs['slices'].astype(str)
    count_data = adata.X.copy()

    if run_normalization:
        total_counts_per_cell = np.sum(adata.X, axis=1)
        sc.pp.log1p(adata)
        
        sc.pp.normalize_total(adata, target_sum=total_counts_per_cell)
        sc.pp.scale(adata, max_value=10)
        sc.pp.highly_variable_genes(adata, n_top_genes=hvgs)
        
        highly_variable_data = adata[:, adata.var['highly_variable']]
        highly_variable_data.obs_names = highly_variable_data.obs_names + '_' + highly_variable_data.obs['slices'].astype(str)
    else:
        highly_variable_data = adata
        highly_variable_data.obs_names = highly_variable_data.obs_names + '_' + highly_variable_data.obs['slices'].astype(str)

    torch.set_num_threads(10)
    memory_usage_value, execution_time, adata = run_model(highly_variable_data, device, nclust, tool, seed)

    # Save the statistics to a CSV file
    stats_file = os.path.join(output_file, f"{sample}_stats.csv")
    with open(stats_file, "w", newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Memory_Usage_MiB", "Execution_Time_s"])
        csvwriter.writerow([memory_usage_value, execution_time])
    
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
