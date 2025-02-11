import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["CUDA_VISIBLE_DEVICES"] = "0"

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
import time
import csv

from memory_profiler import memory_usage
from sklearn.cluster import KMeans
from threadpoolctl import threadpool_limits

def get_windows(job, n_neighbors, tissue_group, exps, X, Y):
    start_time, idx, tissue_name, indices = job
    job_start = time.time()
    print("Starting:", str(idx+1) + '/' + str(len(exps)), ': ' + exps[idx])
    tissue = tissue_group.get_group(tissue_name)
    to_fit = tissue.loc[indices][[X,Y]].values
    fit = NearestNeighbors(n_neighbors=n_neighbors).fit(tissue[[X,Y]].values)
    m = fit.kneighbors(to_fit)
    m = m[0], m[1]
    args = m[0].argsort(axis=1)
    add = np.arange(m[1].shape[0]) * m[1].shape[1]
    sorted_indices = m[1].flatten()[args + add[:, None]]
    neighbors = tissue.index.values[sorted_indices]
    end_time = time.time()
    print("Finishing:", str(idx+1) + "/" + str(len(exps)), ": " + exps[idx], end_time - job_start, end_time - start_time)
    return neighbors.astype(np.int32)

def run_model(cells, ks, n_neighborhoods,X, Y, reg,sum_cols,values,keep_cols, random_seed):
    np.random.seed(random_seed)
    def model_execution():
        n_neighbors = max(ks)
        tissue_group = cells[[X, Y, reg]].groupby(reg)
        exps = list(cells[reg].unique())
        tissue_chunks = [(time.time(), exps.index(t), t, a) for t, indices in tissue_group.groups.items() for a in np.array_split(indices, 1)]
        tissues = [get_windows(job, n_neighbors, tissue_group, exps, X, Y) for job in tissue_chunks]
        out_dict = {}
        for k in ks:
          for neighbors, job in zip(tissues, tissue_chunks):
            chunk = np.arange(len(neighbors)) # indices
            tissue_name = job[2]
            indices = job[3]
            window = values[neighbors[chunk, :k].flatten()].reshape(len(chunk), k, len(sum_cols)).sum(axis=1)
            out_dict[(tissue_name, k)] = (window.astype(np.float16), indices)
        windows = {}
        for k in ks:
          window = pd.concat([pd.DataFrame(out_dict[(exp, k)][0], index=out_dict[(exp, k)][1].astype(int), columns=sum_cols) for exp in exps], axis=0)
          window = window.loc[cells.index.values]
          window = pd.concat([cells[keep_cols], window], axis=1)
          windows[k] = window
        k_centroids = {}
        windows2 = windows[10]
        # km = MiniBatchKMeans(n_clusters=n_neighborhoods, random_state=0)
        # labelskm = km.fit_predict(windows2[sum_cols].values)
        km = KMeans(n_clusters=n_neighborhoods, random_state=random_seed)
        labelskm = km.fit_predict(windows2[sum_cols].values)
        k_centroids[k] = km.cluster_centers_
        cells['predicted_domain'] = labelskm
        cells["predicted_domain"] = cells["predicted_domain"].astype('category')
        return cells
    start_time = time.time()
    mem_usage = memory_usage((model_execution,), max_usage = True, retval=True)
    end_time = time.time()
    adata = mem_usage[1] 
    mem_usage_values = mem_usage[0] 
    return mem_usage_values, end_time - start_time, adata

def process_data(input_file, output_file, sample, nclust, runNormalization=True,hvgs = 4000, n_pcs=50, seed = 0):
    ks = [10]
    X = 'X'
    Y = 'Y'
    reg = 'slices'
    cluster_col = 'ct'
    n_neighborhoods = nclust
    keep_cols = [X, Y, reg, cluster_col]

    # Read data
    adata = sc.read_h5ad(input_file)

    # Preprocessing
    if 'celltype' in adata.obs.columns:
        celltype_codes, unique_types = pd.factorize(adata.obs['celltype'].astype(str))
        adata.obs['ct'] = celltype_codes
        adata.obs['ct'] = adata.obs['ct'].astype('category')
        if runNormalization:
          sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=hvgs)
          sc.pp.normalize_total(adata, inplace=True)
          sc.pp.log1p(adata)
          sc.pp.pca(adata, n_comps=n_pcs)
        else:
          sc.pp.pca(adata, n_comps=n_pcs)
    else:
        if runNormalization:
            sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=hvgs)
            sc.pp.normalize_total(adata, inplace=True)
            sc.pp.log1p(adata)
            sc.pp.pca(adata, n_comps=n_pcs)
            sc.pp.neighbors(adata)
            sc.tl.leiden(adata, resolution=2, key_added='ct', random_state=666)
            adata.obs['ct'] = adata.obs['ct'].astype('category')
        else:
            sc.pp.pca(adata, n_comps=n_pcs)
            sc.pp.neighbors(adata)
            sc.tl.leiden(adata, resolution=2, key_added='ct', random_state=666)
            adata.obs['ct'] = adata.obs['ct'].astype('category')

    cells = pd.DataFrame(adata.obs)
    cells.slices = cells.slices.astype(str)
    cells.reset_index(drop=True, inplace=True)
    cells = pd.concat([cells, pd.get_dummies(cells[cluster_col])], axis=1)
    sum_cols = cells[cluster_col].unique()
    values = cells[sum_cols].values
    
    memory_stats, time_stats,cells = run_model(cells, ks, n_neighborhoods, X, Y, reg, sum_cols,values,keep_cols,seed)
    
    # Write memory usage and execution time to CSV file
    stats_file = os.path.join(output_file, f"{sample}_stats.csv")
    with open(stats_file, "w", newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Memory_Usage_MiB", "Execution_Time_s"])
        csvwriter.writerow([memory_stats, time_stats])
        
    cells["predicted_domain"] = cells["predicted_domain"].astype('category')
    adata.obs["predicted_domain"] = cells["predicted_domain"].values
    adata.obsm["X_embedding"] = adata.obsm['X_pca']
    adata.write(f"{output_file}/{sample}.h5ad")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process data")
    parser.add_argument("--input_file", type=str, help="Input file path")
    parser.add_argument("--output_file", type=str, help="Output file directory")
    parser.add_argument("--sample", type=str, help="Sample name")
    parser.add_argument("--nclust", type=int, help="number of neighborhoods to find")
    parser.add_argument('--runNormalization', type=bool, default=True, help="Run normalization")
    parser.add_argument("--hvgs", type=int, default=4000, help="Number of principal components")
    parser.add_argument("--n_pcs", type=int, default=50, help="Number of principal components")
    parser.add_argument("--seed", type=int, default=0, help="Number of principal components")
    
    
    args = parser.parse_args()

    process_data(args.input_file, args.output_file, args.sample, args.nclust, args.runNormalization, args.hvgs, args.n_pcs, args.seed)
