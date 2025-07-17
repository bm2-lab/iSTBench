import numpy as np
import pandas as pd
import scanpy as sc
import torch
import paste as pst
import anndata as ad
import ot
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse
import random

def process_batches(input_file, batches):
    slices = [sc.read_h5ad(f"{input_file}/{name}.h5ad") for name in batches]
    
    pis = []
    for i in range(len(slices) - 1):
        pi0 = pst.match_spots_using_spatial_heuristic(slices[i].obsm['spatial'],slices[i+1].obsm['spatial'],use_ot=True)
        pi = pst.pairwise_align(slices[i], slices[i+1],alpha=0.1,G_init=pi0,norm=True,verbose=False, use_gpu = True, backend=ot.backend.TorchBackend())
        pis.append(pi)

    new_slices = pst.stack_slices_pairwise(slices, pis)
    
    return new_slices
  


def main(input_file, output_file, batches, step, seed = 123):
  np.random.seed(seed)
  random.seed(seed)
    
  batches = batches.split(",")
  
  new_slices = process_batches(input_file, batches)
  print("Alignment finished")
  
  output_dir = f"{output_file}/PASTE"
  if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
  adata_concat_correct = ad.concat(new_slices, label="slice_name", keys=batches)
  adata_concat_correct.write_h5ad(f"{output_dir}/Spatial_correct_data.h5ad")
  
  fig, ax = plt.subplots(1, len(batches), figsize=(20, 5), gridspec_kw={'wspace': 0.05, 'hspace': 0.1})
  for it in range(len(batches)):
    batch = new_slices[it]
    spatial_coords = batch.obsm['spatial']
    original_domain = batch.obs['original_domain']
    domain_labels = pd.Categorical(original_domain).codes
    scatter = ax[it].scatter(spatial_coords[:, 0], spatial_coords[:, 1], c=domain_labels, cmap='tab10',  s=1)
    ax[it].set_title(batch.obs["slices"].unique().astype(str))
    ax[it].axis('off')
    
  plt.savefig(f"{output_dir}/Slices_correct.pdf")
  
  Z_values = []
  for i, batch in enumerate(new_slices):
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
    parser.add_argument('--output_file', type=str, required=True, help='Output file directory.')
    parser.add_argument('--batches', type=str, required=True, help='Output file directory.')
    parser.add_argument('--step', type=int, required=True, help='Step size.')
    parser.add_argument('--seed', type=int, default=123, help='Random seed.')
    args = parser.parse_args()
    main(args.input_file, args.output_file, args.batches, args.step, args.seed)


