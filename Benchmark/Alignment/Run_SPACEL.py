import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import os
import warnings
warnings.filterwarnings("ignore")
import SPACEL
from SPACEL import Scube
import argparse
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def main(input_file, input_data, output_file, batches, step, seed = 123):
    np.random.seed(seed)
    random.seed(seed)
    
    batches = batches.split(",")
    model = os.path.splitext(os.path.basename(input_data))[0]  
    
    output_dir = f"{output_file}/SPACEL_{model}"
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
    
    predict_data = sc.read_h5ad(input_data)
    predict_data.obs["predicted_domain"] = predict_data.obs["predicted_domain"].astype(str)

    Batch_list = []
    for batche in batches:
      print(batche)
      adata = sc.read_h5ad(f"{input_file}/{batche}.h5ad")
      common_barcodes = adata.obs['barcode'].isin(predict_data.obs['barcode'])
      adata = adata[common_barcodes].copy()
      adata_obs = adata.obs.copy()
      
      predict_data_obs = predict_data.obs[['barcode', "predicted_domain"]].copy()
      merged_obs = adata_obs.merge(predict_data_obs, on='barcode', how='left')
      adata.obs = merged_obs
      Batch_list.append(adata)

    Scube.align(Batch_list, cluster_key='predicted_domain', n_neighbors = 15, n_threads=10, p=2)
    for bl in Batch_list:
      bl.obsm["spatial_aligned"] = bl.obsm["spatial_aligned"].to_numpy()

    adata_concat_correct = ad.concat(Batch_list, label="slice_name", keys=batches)
    adata_concat_correct.obsm["spatial"] = adata_concat_correct.obsm["spatial_aligned"]
    adata_concat_correct.write_h5ad(f"{output_dir}/Spatial_correct_data.h5ad") 
    print("Alognment finished")

    fig, ax = plt.subplots(1, len(batches), figsize=(20, 5), gridspec_kw={'wspace': 0.05, 'hspace': 0.1})
    for it in range(len(batches)):
      batch = Batch_list[it]
      spatial_coords = batch.obsm['spatial']
      original_domain = batch.obs['original_domain']
      domain_labels = pd.Categorical(original_domain).codes
      scatter = ax[it].scatter(spatial_coords[:, 0], spatial_coords[:, 1], c=domain_labels, cmap='tab10',  s=1)
      ax[it].set_title(batch.obs["slices"].unique().astype(str))
      ax[it].axis('off')
      
    plt.savefig(f"{output_dir}/Slices_correct.pdf")
  
    Z_values = []
    for i, batch in enumerate(Batch_list):
      Z_values += [i * step] * batch.shape[0]
    
    adata_concat_correct.obs['Z'] = Z_values
    adata_concat_correct.obs['X'] = adata_concat_correct.obsm['spatial'][:, 0]
    adata_concat_correct.obs['Y'] = adata_concat_correct.obsm['spatial'][:, 1]
    All_coor = adata_concat_correct.obs[['X', 'Y', 'Z']].copy()
    labels = adata_concat_correct.obs['original_domain'].values
    cmap = plt.get_cmap('tab10')
    unique_labels = np.unique(labels)
    color_dict = {label: cmap(i % cmap.N) for i, label in enumerate(unique_labels)}
  
    fig = plt.figure(figsize=(5, 4))
    ax1 = fig.add_subplot(111, projection='3d')
  
    for label in unique_labels:
      mask = adata_concat_correct.obs['original_domain'] == label
      temp_Coor = All_coor[mask]
      color = color_dict[label]
      ax1.scatter3D(
        temp_Coor['X'], temp_Coor['Y'], temp_Coor['Z'],
        c=color,
        s=0.05, marker=".", label=label,
        alpha=1)
    ax1.view_init(elev=20, azim=-40)
    ax1.set_xlabel('')
    ax1.set_ylabel('')
    ax1.set_zlabel('')
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    ax1.set_zticklabels([])

    plt.legend(bbox_to_anchor=(1.2, 0.8), markerscale=10, frameon=False)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/Slices_correct_3d.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run STAligner with given parameters.')
    parser.add_argument('--input_file', type=str, required=True, help='Input file directory.')
    parser.add_argument('--input_data', type=str, required=True, help='Input h5ad data file.')
    parser.add_argument('--output_file', type=str, required=True, help='Output file directory.')
    parser.add_argument('--batches', type=str, required=True, help='Output file directory.')
    parser.add_argument('--step', type=int, required=True, help='Step size.')
    parser.add_argument('--seed', type=int, default=123, help='Random seed.')
    args = parser.parse_args()
    main(args.input_file, args.input_data, args.output_file, args.batches, args.step, args.seed)


