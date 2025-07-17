import warnings
warnings.filterwarnings("ignore")
import STAligner
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import scipy.linalg
import torch
import argparse
import random
import time
from memory_profiler import memory_usage
import csv


def run_staligner_and_cluster(adata_concat, used_device, nclust, sample, output_file):
    def _run():
        adata_concat_ = STAligner.train_STAligner(adata_concat, verbose=True, knn_neigh=50, device=used_device)
        STAligner.mclust_R(adata_concat_, num_cluster=nclust, used_obsm='STAligner')
        return adata_concat_
      
    start_time = time.time()
    mem_usage = memory_usage((_run,), max_usage=True, retval=True)
    end_time = time.time()
    
    mem_used = mem_usage[0]
    adata_concat = mem_usage[1]
    run_time = end_time - start_time
    
    print(f"Memory usage (MiB): {mem_used:.2f}")
    print(f"Execution time (s): {run_time:.2f}")
    
    stats_file = os.path.join(output_file, f"{sample}_stats.csv")
    with open(stats_file, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Memory_Usage_MiB", "Execution_Time_s"])
        writer.writerow([mem_used, run_time])
        
    return adata_concat

def main(input_file, output_file, sample, batches, runNormalization, nclust, hvgs, r,seed = 1234, device = "gpu"):
  # input_file: Each slice needs to be saved separately to an h5ad file, and the path represents the location where these slices are stored. Refer to the example data for the specific format. This is a required parameter.
  # 
  # output_file: Path to the directory where the output results will be saved. This is a required parameter.
  # 
  # sample: Name for the output result files. It is recommended to use the model name, for example, setting it as "STAligner" will result in output files named "STAligner.h5ad". This is a required parameter.
  # 
  # batches: Each slice name should be specified according to its sequential order. This is a required parameter.
  # 
  # hvgs: The number of highly variable genes identified if data normalization is performed. This is an optional parameter, with the default set to 5000.
  # 
  # runNormalization: A boolean variable indicating whether data normalization should be performed before integration. This is an optional parameter, with the default set to TRUE.
  # 
  # nclust: Number of identified domains. This is a required parameter.
  # 
  # r: Radius for identify neighbors, with the default set to 50.
  # 
  # seed: Random seed for model execution. This is an optional parameter, with the default set to 1234.
  # 
  # device: The device on which the model runs can be set to "cpu" and "gpu". This is an optional parameter, with the default set to "gpu".
  
    np.random.seed(seed)
    random.seed(seed)
    batches = batches.split(",")
    
    if device == "gpu":
      used_device = torch.device('cuda:2' if torch.cuda.is_available() else 'cpu')
    else:
      used_device = device

    ###1.Load Data
    Batch_list = []
    adj_list = []
    for batche in batches:
      print(batche)
      adata = sc.read_h5ad(f"{input_file}/{batche}.h5ad")
      adata.var_names_make_unique(join="++")
      adata.obs_names = [x+'_'+batche for x in adata.obs_names]
      
      # Constructing the spatial network
      adata.X = sp.csr_matrix(adata.X)
      STAligner.Cal_Spatial_Net(adata, rad_cutoff = r) # the spatial network are saved in adata.uns[‘adj’]
      
      # Normalization
      if runNormalization:
        sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=hvgs)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata = adata[:, adata.var['highly_variable']]
      adj_list.append(adata.uns['adj'])
      Batch_list.append(adata)
    
    ###2. Concat multiple slices
    adata_concat = ad.concat(Batch_list, label="slice_name", keys=batches)
    adata_concat.obs["batch_name"] = adata_concat.obs["slice_name"].astype('category')
    print('adata_concat.shape: ', adata_concat.shape)
    
    adj_concat = np.asarray(adj_list[0].todense())
    for batch_id in range(1,len(batches)):
      adj_concat = scipy.linalg.block_diag(adj_concat, np.asarray(adj_list[batch_id].todense()))
    adata_concat.uns['edgeList'] = np.nonzero(adj_concat)
    
    ###3.Running STAligner
    adata_concat = run_staligner_and_cluster(adata_concat, used_device, nclust, sample, output_file)
    # adata_concat = STAligner.train_STAligner(adata_concat, verbose=True, knn_neigh = 50, device=used_device)
    # STAligner.mclust_R(adata_concat, num_cluster=nclust, used_obsm='STAligner')
    
    adata_concat.obs["predicted_domain"] = adata_concat.obs["mclust"].astype('category')
    adata_concat.obsm["X_embedding"] = adata_concat.obsm["STAligner"]
    
    adata_concat.uns.pop("edgeList", None)
    adata_concat.write(os.path.join(output_file, f"{sample}.h5ad"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run STAligner with given parameters.')
    parser.add_argument('--input_file', type=str, required=True, help='Input file directory.')
    parser.add_argument('--output_file', type=str, required=True, help='Output file directory.')
    parser.add_argument('--sample', type=str, help='Model name')
    parser.add_argument('--batches', type=str, required=True, help='Slices name in order')
    parser.add_argument('--runNormalization', type=bool, required=True, help='Run normalization.')
    parser.add_argument('--nclust', type=int, help='number of clusters')
    parser.add_argument('--hvgs', type=int, default=5000, help='Number of highly variable genes.')
    parser.add_argument('--r', type=int, default=50, help='Radius for identify neighbors.')
    parser.add_argument('--seed', type=int, default=1234, help='seed')
    parser.add_argument('--device', type=str, default="gpu", help='device')

    args = parser.parse_args()
    main(args.input_file, args.output_file, args.sample, args.batches, args.runNormalization, args.nclust, args.hvgs, args.r, args.seed, args.device)





