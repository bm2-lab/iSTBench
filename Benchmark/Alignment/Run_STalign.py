import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import anndata as ad
import torch
import plotly
import requests
import scanpy as sc
from STalign import STalign
import argparse
import random
import os

def main(input_file, output_file, batches, step, dx = 30, seed = 123):
    np.random.seed(seed)
    random.seed(seed)
    #device = torch.device('cuda:2' if torch.cuda.is_available() else 'cpu')
    batches = batches.split(",")
    
    batch_target = batches[-1]
    
    adata_target = sc.read_h5ad(f"{input_file}/{batch_target}.h5ad")
    target_df = adata_target.obsm["spatial"]
    xJ = np.array(target_df[:, 0])
    yJ = np.array(target_df[:, 1])
    XJ, YJ, J, _ = STalign.rasterize(xJ, yJ, dx=dx)
    
    corrected_adatas = []
    for batch in batches[:-1]:
      adata_source = sc.read_h5ad(f"{input_file}/{batch}.h5ad")
      source_df = adata_source.obsm["spatial"]
      xI = np.array(source_df[:, 0])
      yI = np.array(source_df[:, 1])
      XI, YI, I, _ = STalign.rasterize(xI, yI, dx=dx)
      
      params = {'niter': 10000, 'device': device, 'epV': 50}
      out = STalign.LDDMM([YI, XI], I, [YJ, XJ], J, **params)
      A, v, xv = out['A'], out['v'], out['xv']
    
      dtype = A.dtype
      points_np = np.stack([yI, xI], axis=1)
      points_tensor = torch.tensor(points_np, dtype=dtype).to(device)
      tpointsI = STalign.transform_points_source_to_target(xv, v, A, points_tensor)
      if tpointsI.is_cuda:
        tpointsI = tpointsI.cpu()
     
      adata_source.obsm["spatial"] = tpointsI[:, [1, 0]].numpy()
      corrected_adatas.append(adata_source)
      print(f"{batch} finished!")

    corrected_adatas.append(adata_target)
   
    output_dir = f"{output_file}/STalign"
    if not os.path.exists(output_dir):
     os.makedirs(output_dir)
    
    adata_concat_correct = ad.concat(corrected_adatas, label="slice_name", keys=batches)
    adata_concat_correct.write_h5ad(f"{output_dir}/Spatial_correct_data.h5ad")
   
    fig, ax = plt.subplots(1, len(batches), figsize=(20, 5), gridspec_kw={'wspace': 0.05, 'hspace': 0.1})
    for it in range(len(batches)):
      batch = corrected_adatas[it]
      spatial_coords = batch.obsm['spatial']
      original_domain = batch.obs['original_domain']
      domain_labels = pd.Categorical(original_domain).codes
      scatter = ax[it].scatter(spatial_coords[:, 0], spatial_coords[:, 1], c=domain_labels, cmap='tab10',  s=1)
      ax[it].set_title(batch.obs["slices"].unique().astype(str))
      ax[it].axis('off')
    
    plt.savefig(f"{output_dir}/Slices_correct.pdf")
  
    Z_values = []
    for i, batch in enumerate(corrected_adatas):
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
    parser = argparse.ArgumentParser(description='Run STalign with given parameters.')
    parser.add_argument('--input_file', type=str, required=True, help='Input file directory.')
    parser.add_argument('--output_file', type=str, required=True, help='Output file directory.')
    parser.add_argument('--batches', type=str, required=True, help='Output file directory.')
    parser.add_argument('--step', type=int, required=True, help='Step size.')
    parser.add_argument('--dx', type=int, default=30, help='rasterize at 30um resolution')
    parser.add_argument('--seed', type=int, default=123, help='Random seed.')
    args = parser.parse_args()
    main(args.input_file, args.output_file, args.batches, args.step, args.dx, args.seed)


