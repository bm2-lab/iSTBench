import STAligner
import os
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

# Main function to process spatial transcriptomics data using STAligner.
# The function loads spatial data from different batches, applies domain-guided 3D alignment using STAligner, 
# and saves the aligned data. It also includes data visualization steps for 2D and 3D plots.
def main(input_file, input_data, output_file, batches, landmark_domain, landmark_domain_original, domain, step, runNormalization, hvgs, seed = 123):
    # Set the random seeds for reproducibility
    np.random.seed(seed)
    random.seed(seed)
    
    # Set the device for running the model (CUDA or CPU)
    used_device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    
    # Extract model name from input data file name
    model = os.path.splitext(os.path.basename(input_data))[0]  
    
    # Split batches input into a list
    batches = batches.split(",")
    
    # Prepare output directory for saving results
    output_dir = f"{output_file}/{model}"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    ### 1. Load Data
    # Load prediction data from the input file
    predict_data = sc.read_h5ad(input_data)
    predict_data.obs["original_domain"] = predict_data.obs["original_domain"].astype(str)
    
    # Initialize list to store batch data
    Batch_list = []
    for batche in batches:
        print(batche)
        
        # Load spatial data for each batch
        adata = sc.read_h5ad(f"{input_file}/{batche}.h5ad")
        adata.obs["original_domain"] = adata.obs["original_domain"].astype(str)
        
        adata_obs = adata.obs.copy()
        
        # If domain is predicted, merge with prediction data
        if domain == "predicted_domain":
          predict_data.obs["predicted_domain"] = predict_data.obs["predicted_domain"].astype(str)
          predict_data_obs = predict_data.obs[['barcode', domain]].copy()
          merged_obs = adata_obs.merge(predict_data_obs, on='barcode', how='left')
          adata.obs = merged_obs
          
          # Depending on the model, select the embedding data
          if model == "GraphST" or model == "GraphSTwithPASTE":
            X_embedding = predict_data.obsm["emb"]
          else:
            X_embedding = predict_data.obsm["X_embedding"]
          
          barcode_to_embedding = dict(zip(predict_data.obs['barcode'], X_embedding))
          adata_embedding = np.array([barcode_to_embedding[barcode] for barcode in adata.obs['barcode']])
          adata.obsm['X_embedding'] = adata_embedding
        
        # Convert sparse matrix to CSR format for efficient storage
        adata.X = csr_matrix(adata.X)
        adata.var_names_make_unique(join="++")
        adata.obs_names = adata.obs["barcode"]
        Batch_list.append(adata)
        
    # Concatenate all batches into a single AnnData object
    adata_concat = ad.concat(Batch_list, label="slice_name", keys=batches)
    adata_concat.obs["batch_name"] = adata_concat.obs["slice_name"].astype('category')
    print('adata_concat.shape: ', adata_concat.shape)
    
    ### 2. Running STAligner
    # Store the X_embedding as the initial alignment for STAligner
    adata_concat.obsm["STAligner"] =  adata_concat.obsm["X_embedding"]
    
    ### 3. Spatial domain-guided 3D slices alignment
    # Perform the alignment using ICP (Iterative Closest Point) for each pair of batches
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
        
        # Perform ICP alignment between the target and reference batches
        aligned_coor = STAligner.ICP_align(adata_concat, adata_target, adata_ref, slice_target, slice_ref, [str(landmark_domain)])
        
        # If alignment fails, print a warning; otherwise, update the spatial coordinates
        if np.all(np.isnan(aligned_coor)):
          print("can not correct")
        else:
          adata_target.obsm["spatial"] = aligned_coor
          
    # Concatenate the corrected batches and save the results
    adata_concat_correct = ad.concat(Batch_list, label="slice_name", keys=batches)
    adata_concat_correct.write_h5ad(f"{output_dir}/Spatial_correct_data.h5ad")
    
    ### 4. Visualization - 1
    # Create a 2D plot for visualizing the alignment of spatial coordinates across batches
    fig, ax = plt.subplots(1, len(batches), figsize=(20, 5), gridspec_kw={'wspace': 0.05, 'hspace': 0.1})
    for it in range(len(batches)):
      batch = Batch_list[it]
      spatial_coords = batch.obsm['spatial']
      original_domain = batch.obs['original_domain']
      selected_coords = spatial_coords[original_domain == landmark_domain_original]
      
      # Plot all points in light gray and selected points in red
      ax[it].scatter(spatial_coords[:, 0], spatial_coords[:, 1], c='lightgray', s=1)
      ax[it].scatter(selected_coords[:, 0], selected_coords[:, 1], c='red', s=1)
      ax[it].set_title(batch.obs["slices"].unique().astype(str))
      ax[it].axis('off')
    
    # Save the 2D visualization
    plt.savefig(f"{output_dir}/Slices_correct.pdf")
    
    ### 5. Visualization - 2
    # Create a 3D plot for visualizing the spatial alignment across batches
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
    
    # Plot the 3D coordinates for each domain
    for it, label in enumerate(np.unique(adata_concat_correct.obs['original_domain'])):
        temp_Coor = All_coor.loc[adata_concat.obs['original_domain'] == label, :]
        temp_xd = temp_Coor['X']
        temp_yd = temp_Coor['Y']
        temp_zd = temp_Coor['Z']
        if label == landmark_domain_original:
            ax1.scatter3D(temp_xd, temp_yd, temp_zd, c="red", s=0.05, marker=".", label=label, alpha=1)
        else:
            ax1.scatter3D(temp_xd, temp_yd, temp_zd, c=color_dict[label], s=0.05, marker=".", label=label, alpha=0.6)
    
    # Finalize the 3D plot
    plt.legend(bbox_to_anchor=(1.2, 0.8), markerscale=10, frameon=False)
    ax1.elev = 20
    ax1.azim = -40
    ax1.set_xlabel('')
    ax1.set_ylabel

('')
    ax1.set_zlabel('')
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    ax1.set_zticklabels([])
    
    # Save the 3D visualization
    plt.savefig(f"{output_dir}/Slices_correct_3d.pdf")


# Command-line interface setup for running the script
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run STAligner with given parameters.')
    parser.add_argument('--input_file', type=str, required=True, help='Input file directory.')
    parser.add_argument('--input_data', type=str, required=True, help='Input h5ad data file.')
    parser.add_argument('--output_file', type=str, required=True, help='Output file directory.')
    parser.add_argument('--batches', type=str, required=True, help='List of batches to process.')
    parser.add_argument('--landmark_domain', type=str, required=True, help='Landmark domain for alignment.')
    parser.add_argument('--landmark_domain_original', type=str, required=True, help='Original domain name for visualization.')
    parser.add_argument('--domain', type=str, required=True, help='Domain column name in the data.')
    parser.add_argument('--step', type=int, required=True, help='Step size for visualization.')
    parser.add_argument('--runNormalization', type=bool, required=True, help='Flag to indicate normalization.')
    parser.add_argument('--hvgs', type=int, default=10000, help='Number of highly variable genes (default: 10000).')

    args = parser.parse_args()
    main(args.input_file, args.input_data, args.output_file, args.batches, args.landmark_domain, args.landmark_domain_original, args.domain, args.step, args.runNormalization, args.hvgs)


