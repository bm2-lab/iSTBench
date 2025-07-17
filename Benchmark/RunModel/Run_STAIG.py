import argparse
import os
os.chdir('/Benchmark/RunModel/STAIG')
import warnings
warnings.filterwarnings('ignore')
import random
import yaml
from yaml import SafeLoader
import torch
from staig.adata_processing import LoadBatchAdata
import numpy as np
from staig.staig import STAIG
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
import ot
from tqdm import tqdm
from anndata import AnnData
from scipy.sparse import csc_matrix, csr_matrix
from scipy.spatial.distance import euclidean, cosine
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.special import softmax
from scipy.linalg import block_diag
import time
from memory_profiler import memory_usage
import csv

###load data function
def process_batch_adata(file_list, n_top_genes=1000, n_neighbors=5, image_emb=False, runNormalization=True, nPCA=64):
    def construct_interaction(input_adata, n_neighbors):
        position = input_adata.obsm['spatial']
        distance_matrix = ot.dist(position, position, metric='euclidean')
        n_spot = distance_matrix.shape[0]
        interaction = np.zeros([n_spot, n_spot])
        for i in range(n_spot):
            vec = distance_matrix[i, :]
            distance = vec.argsort()
            for t in range(1, n_neighbors + 1):
                y = distance[t]
                interaction[i, y] = 1
                
        adj = interaction + interaction.T
        adj = np.where(adj > 1, 1, adj)
        input_adata.obsm['local_graph'] = adj
        return input_adata
      
    def load_data(file_list, n_neighbors, runNormalization, n_top_genes):
        adata_list = []
        for i in file_list:
            print(i)
            if runNormalization:
              sc.pp.highly_variable_genes(i, flavor="seurat_v3", n_top_genes=n_top_genes)
              sc.pp.normalize_total(i, target_sum=1e4)
              sc.pp.log1p(i)
            else:
              i.var['highly_variable'] = True
            adata = construct_interaction(input_adata=i, n_neighbors = n_neighbors)
            adata_list.append(adata)           
        print('load all slices done')
        return adata_list
      
    def concatenate_slices(adata_list):
        highly_variable_genes_set = set(adata_list[0].var['highly_variable'][adata_list[0].var['highly_variable']].index)
        
        for adata in adata_list[1:]:
            current_set = set(adata.var['highly_variable'][adata.var['highly_variable']].index)
            highly_variable_genes_set = highly_variable_genes_set.intersection(current_set)
        adata = AnnData.concatenate(*adata_list, join='outer')
        
        adata_Vars = adata[:, adata.var.index.isin(highly_variable_genes_set)]
        if isinstance(adata_Vars.X, csc_matrix) or isinstance(adata_Vars.X, csr_matrix):
            feat = adata_Vars.X.toarray()[:, ]
        else:
            feat = adata_Vars.X[:, ]
            
        adata.obsm['feat'] = feat
        merged_adata = adata
        print('merge done')
        return merged_adata
      
    def construct_whole_graph(adata_list, merged_adata):
        matrix_list = [i.obsm['local_graph'] for i in adata_list]
        adjacency = block_diag(*matrix_list)
        merged_adata.obsm['graph_neigh'] = adjacency
        
        mask_list = [np.ones_like(i.obsm['local_graph'], dtype=int) for i in adata_list]
        mask = block_diag(*mask_list)
        merged_adata.obsm['mask_neigh'] = mask
        
    def calculate_edge_weights(merged_adata):
        graph_neigh = merged_adata.obsm['graph_neigh']
        node_emb = merged_adata.obsm['img_emb']
        num_nodes = node_emb.shape[0]
        edge_weights = np.zeros_like(graph_neigh)  
        
        for i in tqdm(range(num_nodes), desc="Calculating distances"):  
            for j in range(num_nodes):
                if graph_neigh[i, j] == 1:  
                    edge_weights[i, j] = euclidean(node_emb[i], node_emb[j])
                    
        edge_probabilities = np.zeros_like(edge_weights)
        for i in tqdm(range(num_nodes), desc="Calculating edge_probabilities"):
            non_zero_indices = edge_weights[i] != 0
            if non_zero_indices.any():  
                non_zero_weights = np.log(edge_weights[i][non_zero_indices]) 
                softmax_weights = softmax(non_zero_weights)
                edge_probabilities[i][non_zero_indices] = softmax_weights
                
        merged_adata.obsm['edge_probabilities'] = edge_probabilities
        
    def calculate_edge_weights_gene(merged_adata, nPCA):
      
        graph_neigh = merged_adata.obsm['graph_neigh']
        node_emb = merged_adata.obsm['feat']
        scaler = StandardScaler()
        embedding = scaler.fit_transform(node_emb)
        pca = PCA(n_components=nPCA, random_state=42)
        embedding = pca.fit_transform(embedding)
        node_emb = embedding
        
        num_nodes = node_emb.shape[0]
        edge_weights = np.zeros((num_nodes, num_nodes))
        
        for i in tqdm(range(num_nodes), desc="Calculating distances"):
            for j in range(num_nodes):
                if graph_neigh[i, j] == 1:  
                    edge_weights[i, j] = cosine(node_emb[i], node_emb[j])
                    
        edge_probabilities = np.zeros_like(edge_weights)
        for i in range(num_nodes):
            non_zero_indices = edge_weights[i] != 0
            if non_zero_indices.any():
                non_zero_weights = edge_weights[i][non_zero_indices]
                softmax_weights = softmax(non_zero_weights)
                edge_probabilities[i][non_zero_indices] = softmax_weights
                
        merged_adata.obsm['edge_probabilities'] = edge_probabilities
        
    adata_list = load_data(file_list, n_neighbors, runNormalization, n_top_genes)
    merged_adata = concatenate_slices(adata_list)
    construct_whole_graph(adata_list, merged_adata)
    if image_emb:
      calculate_edge_weights(merged_adata)
    else:
      calculate_edge_weights_gene(merged_adata, nPCA)
    
    return merged_adata

def run_staig_and_cluster(data, args, config, output_file, sample):
    def _run():
        staig = STAIG(args=args, config=config, single=True, refine=False)
        staig.adata = data
        staig.train()
        staig.eva()
        staig.cluster(args.label)
        return staig.adata

    start_time = time.time()
    mem_usage = memory_usage((_run,), max_usage=True, retval=True)
    end_time = time.time()
    
    mem_used = mem_usage[0]  # in MiB
    adata_out = mem_usage[1]
    run_time = end_time - start_time
    
    print(f"Memory usage (MiB): {mem_used:.2f}")
    print(f"Execution time (s): {run_time:.2f}")
    
    stats_file = os.path.join(output_file, f"{sample}_stats.csv")
    with open(stats_file, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Memory_Usage_MiB", "Execution_Time_s"])
        writer.writerow([mem_used, run_time])
        
    return adata_out

def main(input_file, output_file, sample, runNormalization, nclust, hvgs, nPCA, num_layers, k_neighbor, tau, num_epochs,seed = 1234, device = "gpu"):
  # input_file: Each slice needs to be saved separately to an h5ad file, and the path represents the location where these slices are stored. Refer to the example data for the specific format. This is a required parameter.
  # 
  # output_file: Path to the directory where the output results will be saved. This is a required parameter.
  # 
  # model_path: Path to store the intermediate results. This is a required parameter.
  # 
  # sample: Name for the output result files. It is recommended to use the model name, for example, setting it as "STAIG" will result in output files named "STAIG.h5ad". This is a required parameter.
  # 
  # hvgs: The number of highly variable genes identified if data normalization is performed. This is an optional parameter, with the default set to 3000.
  # 
  # runNormalization: A boolean variable indicating whether data normalization should be performed before integration. This is an optional parameter, with the default set to TRUE.
  # 
  # nclust: Number of identified domains. This is a required parameter.
  # 
  # nPCA: Number of PCA, with the default set to 64.
  # 
  # num_layers: Number of GCN layers.
  # 
  # k_neighbor: Number of neighboring nodes for each spot.
  # 
  # tau: Temperature parameter.
  # 
  # num_epochs: Number of epochs.
  # 
  # The original paper provides reference values for parameters such as num_layers, k_neighbor, tau, and num_epochs across different data types.

seed: Random seed for model execution. This is an optional parameter, with the default set to 1234.

device: The device on which the model runs can be set to "cpu" and "gpu". This is an optional parameter, with the default set to "gpu".
  
  if device == "gpu":
    used_device = torch.device('cuda:2' if torch.cuda.is_available() else 'cpu')
  else:
    used_device = device
  
  torch.manual_seed(seed)
  np.random.seed(seed)
  if torch.cuda.is_available():
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
  random.seed(12345)
  os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8' 
  torch.use_deterministic_algorithms(True)
  
  ###1. Load data
  batches = os.listdir(input_file)
  batches = [batch.split('.')[0] for batch in batches]
  filelist = []
  for batch in batches:
    adata_batch = sc.read_h5ad(f"{input_file}/{batch}.h5ad")
    filelist.append(adata_batch)
    
  data = process_batch_adata(file_list=filelist, n_top_genes=hvgs, n_neighbors=k_neighbor, runNormalization=runNormalization, nPCA=nPCA)
  columns_to_drop = [col for col in data.var.columns if col not in ['means', 'variances']]
  data.var = data.var.drop(columns=columns_to_drop)
  print(data)
  
  ###2. Run STAIG
  args = argparse.Namespace(
    dataset='na',
    slide='na',
    config='na',
    label=False
    )
  config = {
    'learning_rate': 0.0005,
    'num_hidden': 64,
    'num_proj_hidden': 64,
    'activation': 'prelu',
    'base_model': 'GCNConv',
    'num_layers': num_layers,
    'drop_edge_rate_1': 0.1,
    'drop_edge_rate_2': 0.3,
    'drop_feature_rate_1': 0.01,
    'drop_feature_rate_2': 0.02,
    'tau': tau,
    'tau_decay': 0.00,
    'num_epochs': num_epochs,
    'weight_decay': 1e-5,
    'num_clusters': nclust,
    'num_gene': data.obsm['feat'].shape[1] ,
    'device': used_device
    }
  
  data = run_staig_and_cluster(data, args, config, output_file, sample)

  ###3. writing result
  data.obs["predicted_domain"] = data.obs["mclust"].astype('category')
  data.obsm["X_embedding"] = data.obsm["emb"]
  data.write(os.path.join(output_file, f"{sample}.h5ad"))
  
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Run STAIG with given parameters.')
  parser.add_argument('--input_file', type=str, required=True, help='Input file directory.')
  parser.add_argument('--output_file', type=str, required=True, help='Output file directory.')
  parser.add_argument('--sample', type=str, help='Model name')
  parser.add_argument('--runNormalization', type=bool, required=True, help='Run normalization.')
  parser.add_argument('--nclust', type=int, help='number of clusters')
  parser.add_argument('--hvgs', type=int, default=3000, help='Number of highly variable genes.')
  parser.add_argument('--nPCA', type=int, default=64, help='Number of PCA, default is 64.')
  parser.add_argument('--num_layers', type=int, default=1, help='Number of GCN layers.')
  parser.add_argument('--k_neighbor', type=int, default=5, help='Number of neighboring nodes 𝑘 for each spot.')
  parser.add_argument('--tau', type=int, default=10, help='Temperature parameter.')
  parser.add_argument('--num_epochs', type=int, default=400, help='Number of epochs.')
  parser.add_argument('--seed', type=int, default=1234, help='seed')
  parser.add_argument('--device', type=str, default="gpu", help='device')
  args = parser.parse_args()
  main(args.input_file, args.output_file, args.sample, args.runNormalization, args.nclust, args.hvgs, args.nPCA, args.num_layers, args.k_neighbor, args.tau, args.num_epochs,args.seed, args.device)


