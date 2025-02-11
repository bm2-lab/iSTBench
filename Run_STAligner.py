import STAligner
import os
#os.environ['R_HOME'] = '/opt/R/4.3.2/lib/R'
#os.environ['R_USER'] = "/home/dongkj/anaconda3/envs/MultiSpatial/lib/python3.10/site-packages/rpy2"
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import scipy.linalg
from scipy.sparse import csr_matrix
import torch
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse
import random

def main(input_file, input_data, output_file, batches, landmark_domain, landmark_domain_original,domain, step, runNormalization, hvgs, seed = 123):
    np.random.seed(seed)
    random.seed(seed)
    
    used_device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    
    model = os.path.splitext(os.path.basename(input_data))[0]  
    # batches = os.listdir(input_file)
    # batches = [batch.split('.')[0] for batch in batches]
    # batches = sorted(batches, key=lambda x: int(x.replace('slices', '')))
    batches = batches.split(",")
    
    output_dir = f"{output_file}/{model}"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    ###1.Load Data
    predict_data = sc.read_h5ad(input_data)
    predict_data.obs["original_domain"] = predict_data.obs["original_domain"].astype(str)
    Batch_list = []
    #adj_list = []
    for batche in batches:
        print(batche)
        adata = sc.read_h5ad(f"{input_file}/{batche}.h5ad")
        adata.obs["original_domain"] = adata.obs["original_domain"].astype(str)
        
        adata_obs = adata.obs.copy()
        if domain == "predicted_domain":
          predict_data.obs["predicted_domain"] = predict_data.obs["predicted_domain"].astype(str)
          predict_data_obs = predict_data.obs[['barcode', domain]].copy()
          merged_obs = adata_obs.merge(predict_data_obs, on='barcode', how='left')
          adata.obs = merged_obs
          if model == "GraphST" or model == "GraphSTwithPASTE":
            X_embedding = predict_data.obsm["emb"]
          else:
            X_embedding = predict_data.obsm["X_embedding"]
          barcode_to_embedding = dict(zip(predict_data.obs['barcode'], X_embedding))
          adata_embedding = np.array([barcode_to_embedding[barcode] for barcode in adata.obs['barcode']])
          adata.obsm['X_embedding'] = adata_embedding
        adata.X = csr_matrix(adata.X)
        adata.var_names_make_unique(join="++")
        adata.obs_names = adata.obs["barcode"]
        #STAligner.Cal_Spatial_Net(adata, rad_cutoff=50)
        # if runNormalization:
        #     sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=hvgs)
        #     sc.pp.normalize_total(adata, target_sum=1e4)
        #     sc.pp.log1p(adata)
        #     adata = adata[:, adata.var['highly_variable']]
        #adj_list.append(adata.uns['adj'])
        Batch_list.append(adata)
        
    adata_concat = ad.concat(Batch_list, label="slice_name", keys=batches)
    adata_concat.obs["batch_name"] = adata_concat.obs["slice_name"].astype('category')
    print('adata_concat.shape: ', adata_concat.shape)
    
    ###2.Running STAligner
    adata_concat.obsm["STAligner"] =  adata_concat.obsm["X_embedding"]
    # iter_comb1 = [(i, i+1) for i in range(0, len(batches)-1)]
    # adata_concat2 = STAligner.train_STAligner_subgraph(adata_concat, verbose=True, knn_neigh = 100, n_epochs = 600, iter_comb = iter_comb1,
    #                                                         Batch_list=Batch_list, device=used_device)
    
    ###3.Spatial domain guided 3D slices alignment
    adata_concat.obs["louvain"] = adata_concat.obs[domain].astype(str)
    for a in Batch_list:
      a.obs['louvain'] = a.obs['predicted_domain']
    
    iter_comb2 = [(i, len(batches) - 1) for i in range(len(batches) - 1)]
    for comb in iter_comb2:
        print(comb)
        i, j = comb[0], comb[1]
        adata_target = Batch_list[i]
        adata_ref = Batch_list[j]
        slice_target = batches[i]
        slice_ref = batches[j]
        aligned_coor = STAligner.ICP_align(adata_concat, adata_target, adata_ref, slice_target, slice_ref, [str(landmark_domain)])
        if np.all(np.isnan(aligned_coor)):
          print("can not correct")
        else:
          adata_target.obsm["spatial"] = aligned_coor
          
    adata_concat_correct = ad.concat(Batch_list, label="slice_name", keys=batches)
    adata_concat_correct.write_h5ad(f"{output_dir}/Spatial_correct_data.h5ad")
    
    ###4.可视化-1
    fig, ax = plt.subplots(1, len(batches), figsize=(20, 5), gridspec_kw={'wspace': 0.05, 'hspace': 0.1})
    for it in range(len(batches)):
      batch = Batch_list[it]
      spatial_coords = batch.obsm['spatial']
      original_domain = batch.obs['original_domain']
      selected_coords = spatial_coords[original_domain == landmark_domain_original]
      ax[it].scatter(spatial_coords[:, 0], spatial_coords[:, 1], c='lightgray', s=1)  # 绘制所有点
      ax[it].scatter(selected_coords[:, 0], selected_coords[:, 1], c='red', s=1)  # 绘制选定类别的点
      ax[it].set_title(batch.obs["slices"].unique().astype(str))
      ax[it].axis('off')
    plt.savefig(f"{output_dir}/Slices_correct.pdf")
    
    ###4.可视化-2
    Z_values = []
    for i, batch in enumerate(Batch_list):
        Z_values += [i * step] * batch.shape[0]
    adata_concat_correct.obs['Z'] = Z_values
    adata_concat_correct.obs["X"] = adata_concat_correct.obsm["spatial"][:,0]
    adata_concat_correct.obs["Y"] = adata_concat_correct.obsm["spatial"][:,1]
    All_coor = adata_concat_correct.obs[['X', 'Y', 'Z']].copy()
    
    fig = plt.figure(figsize=(5, 4))
    ax1 = plt.axes(projection='3d')
    cmap = plt.get_cmap('tab10') 
    color_dict = {label: cmap(i % cmap.N) for i, label in enumerate(np.unique(adata_concat.obs['original_domain']))}
    
    for it, label in enumerate(np.unique(adata_concat_correct.obs['original_domain'])):
        temp_Coor = All_coor.loc[adata_concat.obs['original_domain'] == label, :]
        temp_xd = temp_Coor['X']
        temp_yd = temp_Coor['Y']
        temp_zd = temp_Coor['Z']
        if label == landmark_domain_original:
            ax1.scatter3D(temp_xd, temp_yd, temp_zd, c = "red",
                          s=0.05, marker=".", label=label, alpha=1)
        else:
            ax1.scatter3D(temp_xd, temp_yd, temp_zd, c = color_dict[label],
                          s=0.05, marker=".", label=label, alpha=0.6)
                          
    plt.legend(bbox_to_anchor=(1.2, 0.8), markerscale=10, frameon=False)
    ax1.elev = 20
    ax1.azim = -40
    ax1.set_xlabel('')
    ax1.set_ylabel('')
    ax1.set_zlabel('')
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    ax1.set_zticklabels([])
    plt.savefig(f"{output_dir}/Slices_correct_3d.pdf")
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run STAligner with given parameters.')
    parser.add_argument('--input_file', type=str, required=True, help='Input file directory.')
    parser.add_argument('--input_data', type=str, required=True, help='Input h5ad data file.')
    parser.add_argument('--output_file', type=str, required=True, help='Output file directory.')
    parser.add_argument('--batches', type=str, required=True, help='Output file directory.')
    parser.add_argument('--landmark_domain', type=str, required=True, help='Landmark domain.')
    parser.add_argument('--landmark_domain_original', type=str, required=True, help='Landmark domain original name.')
    parser.add_argument('--domain', type=str, required=True, help='Domain column name.')
    parser.add_argument('--step', type=int, required=True, help='Step size.')
    parser.add_argument('--runNormalization', type=bool, required=True, help='Run normalization.')
    parser.add_argument('--hvgs', type=int, default=10000, help='Number of highly variable genes.')

    args = parser.parse_args()
    main(args.input_file, args.input_data, args.output_file, args.batches, args.landmark_domain, args.landmark_domain_original, args.domain, args.step, args.runNormalization, args.hvgs)
