import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import torch
import scipy as sp
import sklearn.neighbors
from sklearn.metrics.pairwise import euclidean_distances
from spiral.main import SPIRAL_integration
from spiral.layers import MeanAggregator
from spiral.utils import layer_map
import time
import csv
from memory_profiler import memory_usage
import argparse

# ------------------------
# Data Preparation
# ------------------------
def prepare_input_data(input_file, output_dir, hvgs=1000, runNormalization=True):
    sample_files = os.listdir(input_file)
    sample_names = [f.split('.')[0] for f in sample_files]

    VF, MAT, meta_combine, coord_combine, adata_batch_list = [], [], [], [], []

    for name in sample_names:
        adata = sc.read_h5ad(f"{input_file}/{name}.h5ad")
        adata.var_names_make_unique()
        adata.obs['Ground Truth'] = adata.obs["original_domain"]
        adata.obs_names = [x + '_' + name for x in adata.obs_names]

        sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=hvgs)
        if runNormalization:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)

        adata.obs['batch'] = adata.obs['slices']
        mat = pd.DataFrame(adata.X, columns=adata.var_names, index=adata.obs_names)
        coord = pd.DataFrame(adata.obsm['spatial'], columns=['x', 'y'], index=adata.obs_names)
        meta = adata.obs[['Ground Truth', 'batch']].rename(columns={"Ground Truth": "celltype"})

        adata_batch_list.append(adata)
        MAT.append(mat)
        meta_combine.append(meta)
        coord_combine.append(coord)
        VF = np.union1d(VF, adata.var_names[adata.var['highly_variable']])

        meta.to_csv(f"{output_dir}/{name}_label-1.txt")
        coord.to_csv(f"{output_dir}/{name}_positions-1.txt")

    return sample_names, VF, MAT, meta_combine, coord_combine, adata_batch_list

def export_features(sample_names, MAT, VF, output_dir):
    feat_mat = []
    for i, name in enumerate(sample_names):
        mat = MAT[i].loc[:, VF]
        mat.to_csv(f"{output_dir}/{name}_features-1.txt")
        feat_mat.append(mat)
    return feat_mat

def build_knn_graphs(sample_names, feat_mat, coord_combine, meta_combine, output_dir, KNN=6):
    def Cal_Spatial_Net(adata, k_cutoff=6):
        coor = pd.DataFrame(adata.obsm['spatial'], index=adata.obs.index)
        nbrs = sklearn.neighbors.NearestNeighbors(n_neighbors=k_cutoff+1).fit(coor)
        distances, indices = nbrs.kneighbors(coor)
        edges = pd.DataFrame([
            (i, indices[i, j], distances[i, j]) 
            for i in range(indices.shape[0]) 
            for j in range(indices.shape[1])
        ], columns=['Cell1', 'Cell2', 'Distance'])
        edges = edges[edges['Distance'] > 0]
        trans = dict(zip(range(len(coor)), coor.index))
        edges['Cell1'] = edges['Cell1'].map(trans)
        edges['Cell2'] = edges['Cell2'].map(trans)
        return edges

    for i, name in enumerate(sample_names):
        adata = sc.AnnData(feat_mat[i])
        adata.X = sp.sparse.csr_matrix(adata.X)
        adata.var_names_make_unique()
        adata.obsm['spatial'] = coord_combine[i][['x', 'y']].values
        edges = Cal_Spatial_Net(adata, k_cutoff=KNN)
        np.savetxt(f"{output_dir}/{name}_edge_KNN_{KNN}.csv", edges.values[:, :2], fmt='%s')

def configure_parser(N, M, KNN, lamda, seed):
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, default=seed)
    parser.add_argument('--AEdims', type=list, default=[N, [512], 32])
    parser.add_argument('--AEdimsR', type=list, default=[32, [512], N])
    parser.add_argument('--GSdims', type=list, default=[512, 32])
    parser.add_argument('--zdim', type=int, default=32)
    parser.add_argument('--znoise_dim', type=int, default=4)
    parser.add_argument('--CLdims', type=list, default=[4, [], M])
    parser.add_argument('--DIdims', type=list, default=[28, [32, 16], M])
    parser.add_argument('--beta', type=float, default=1.0)
    parser.add_argument('--agg_class', type=str, default=MeanAggregator)
    parser.add_argument('--num_samples', type=int, default=KNN)
    parser.add_argument('--N_WALKS', type=int, default=KNN)
    parser.add_argument('--WALK_LEN', type=int, default=1)
    parser.add_argument('--N_WALK_LEN', type=int, default=KNN)
    parser.add_argument('--NUM_NEG', type=int, default=KNN)
    parser.add_argument('--epochs', type=int, default=100) 
    parser.add_argument('--batch_size', type=int, default=1024)
    parser.add_argument('--lr', type=float, default=1e-3)
    parser.add_argument('--weight_decay', type=float, default=5e-4)
    parser.add_argument('--alpha1', type=float, default=N)
    parser.add_argument('--alpha2', type=float, default=1)
    parser.add_argument('--alpha3', type=float, default=1)
    parser.add_argument('--alpha4', type=float, default=1)
    parser.add_argument('--lamda', type=float, default=lamda)
    parser.add_argument('--Q', type=float, default=10)
    return parser.parse_args([])

def find_optimal_resolution(adata, target_clusters, max_iter=100, tol=1):
    low, high = 0, 5.0
    for _ in range(max_iter):
        mid = (low + high) / 2
        sc.tl.louvain(adata, resolution=mid, key_added="predicted_domain")
        n = adata.obs['predicted_domain'].nunique()
        if abs(n - target_clusters) < tol:
            return mid
        if n < target_clusters:
            low = mid
        else:
            high = mid
    return mid

def run_spiral(params, feat_file, edge_file, meta_file):
    SPII = SPIRAL_integration(params, feat_file, edge_file, meta_file)
    SPII.train()
    return SPII

def main(input_file, output_file, model_path, sample="SPIRAL", hvgs=1000, runNormalization=True, KNN=6, data="DLPFC_sample3", seed=1234, lamda=1, nclust=7):
  # input_file: Each slice needs to be saved separately to an h5ad file, and the path represents the location where these slices are stored. Refer to the example data for the specific format. This is a required parameter.
  # 
  # output_file: Path to the directory where the output results will be saved. This is a required parameter.
  # 
  # model_path: Path to store the intermediate results. This is a required parameter.
  # 
  # sample: Name for the output result files. It is recommended to use the model name, for example, setting it as "SPIRAL" will result in output files named "SPIRAL.h5ad". This is a required parameter.
  # 
  # hvgs: The number of highly variable genes identified if data normalization is performed. This is an optional parameter, with the default set to 1000.
  # 
  # runNormalization: A boolean variable indicating whether data normalization should be performed before integration. This is an optional parameter, with the default set to TRUE.
  # 
  # KNN: In constructing the k-nearest neighbor graph, the choice of k determines the number of nearest neighbors considered for each cell or spot, with the default set to 6.
  # 
  # data: The name of the data. This is a required parameter.
  # 
  # seed: Random seed for model execution. This is an optional parameter, with the default set to 1234.
  # 
  # nclust: Number of identified domains. This is a required parameter.
  
    model_path = f"{model_path}/{data}"
    os.makedirs(model_path, exist_ok=True)
    os.environ['CUDA_VISIBLE_DEVICES'] = '1'

    sample_names, VF, MAT, meta_combine, coord_combine, adata_list = prepare_input_data(input_file, model_path, hvgs, runNormalization)
    adata_combine = anndata.concat(adata_list, join="inner")
    feat_mat = export_features(sample_names, MAT, VF, model_path)
    build_knn_graphs(sample_names, feat_mat, coord_combine, meta_combine, model_path, KNN)

    feat_file = [f"{model_path}/{n}_features-1.txt" for n in sample_names]
    edge_file = [f"{model_path}/{n}_edge_KNN_{KNN}.csv" for n in sample_names]
    meta_file = [f"{model_path}/{n}_label-1.txt" for n in sample_names]
    coord_file = [f"{model_path}/{n}_positions-1.txt" for n in sample_names]

    N = len(VF)
    M = 1 if len(sample_names) == 2 else len(sample_names)
    params = configure_parser(N, M, KNN, lamda, seed)

    start_time = time.time()
    mem_usage, SPII = memory_usage((run_spiral, (params, feat_file, edge_file, meta_file)), retval=True, max_usage=True)
    end_time = time.time()

    run_time = end_time - start_time
    mem_used = mem_usage

    print(f"Memory usage (MiB): {mem_used:.2f}")
    print(f"Execution time (s): {run_time:.2f}")

    stats_file = os.path.join(output_file, f"{sample}_stats.csv")
    with open(stats_file, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Memory_Usage_MiB", "Execution_Time_s"])
        writer.writerow([f"{mem_used:.2f}", f"{run_time:.2f}"])

    SPII.model.eval()
    all_idx = np.arange(SPII.feat.shape[0])
    all_layer, all_mapping = layer_map(all_idx.tolist(), SPII.adj, len(SPII.params.GSdims))
    all_rows = SPII.adj.tolil().rows[all_layer[0]]
    all_feature = torch.Tensor(SPII.feat.iloc[all_layer[0], :].values).float().cuda()
    all_embed, _, _, _ = SPII.model(all_feature, all_layer, all_mapping, all_rows, SPII.params.lamda, SPII.de_act, SPII.cl_act)
    embed = all_embed[2].cpu().detach().numpy()
    embed_df = pd.DataFrame(embed, index=SPII.feat.index, columns=[f'GTT_{i}' for i in range(embed.shape[1])])

    ann = anndata.AnnData(SPII.feat)
    ann.obsm['spiral'] = embed_df.iloc[:, SPII.params.znoise_dim:].values
    sc.pp.neighbors(ann, use_rep='spiral')
    resolution = find_optimal_resolution(ann, nclust)
    sc.tl.louvain(ann, resolution=resolution, key_added="predicted_domain")

    embed_df_aligned = embed_df.reindex(adata_combine.obs_names)
    adata_combine.obsm["X_embedding"] = embed_df_aligned.values
    adata_combine.obs["predicted_domain"] = ann.obs["predicted_domain"].reindex(adata_combine.obs_names)
    adata_combine.write(f"{output_file}/{sample}.h5ad")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run SPIRAL integration with profiling.')
    parser.add_argument('--input_file', type=str, required=True)
    parser.add_argument('--output_file', type=str, required=True)
    parser.add_argument('--model_path', type=str, required=True)
    parser.add_argument('--sample', type=str, default='SPIRAL')
    parser.add_argument('--hvgs', type=int, default=1000)
    parser.add_argument('--runNormalization', type=bool, default=True)
    parser.add_argument('--KNN', type=int, default=6)
    parser.add_argument('--data', type=str, required=True)
    parser.add_argument('--seed', type=int, default=1234)
    parser.add_argument('--lamda', type=float, default=1.0)
    parser.add_argument('--nclust', type=int, required=True)

    args = parser.parse_args()
    main(args.input_file, args.output_file, args.model_path, args.sample, args.hvgs, args.runNormalization, args.KNN, args.data, args.seed, args.lamda, args.nclust)
