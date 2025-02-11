import os
import warnings
from datetime import datetime
import random
import gdown
import time
from memory_profiler import memory_usage
import csv

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import squidpy as sq
from matplotlib import gridspec
from sklearn.preprocessing import MinMaxScaler
import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse as sp

from nichecompass.models import NicheCompass
from nichecompass.utils import (add_gps_from_gp_dict_to_adata,
                                filter_and_combine_gp_dict_gps,
                                extract_gp_dict_from_omnipath_lr_interactions,
                                extract_gp_dict_from_nichenet_lrt_interactions,
                                extract_gp_dict_from_mebocost_es_interactions)

import argparse
import mygene

def train_model(adata, counts_key, adj_key, cat_covariates_keys, gp_names_key, active_gp_names_key,
                gp_targets_mask_key, gp_targets_categories_mask_key, gp_sources_mask_key, gp_sources_categories_mask_key,
                latent_key, conv_layer_encoder, active_gp_thresh_ratio, n_epochs, n_epochs_all_gps, lr, lambda_edge_recon,
                lambda_gene_expr_recon, lambda_l1_masked, lambda_l1_addon, edge_batch_size, n_sampled_neighbors, use_cuda_if_available,cat_covariates_embeds_nums, seed):
    
    def model_execution():
      model = NicheCompass(adata,
                      counts_key=counts_key,
                      adj_key=adj_key,
                      cat_covariates_embeds_injection=["gene_expr_decoder"],
                      cat_covariates_keys=cat_covariates_keys,
                      cat_covariates_no_edges=[True],
                      cat_covariates_embeds_nums=cat_covariates_embeds_nums,
                      gp_names_key=gp_names_key,
                      active_gp_names_key=gp_names_key,
                      gp_targets_mask_key=gp_targets_mask_key,
                      gp_targets_categories_mask_key=gp_targets_categories_mask_key,
                      gp_sources_mask_key=gp_sources_mask_key,
                      gp_sources_categories_mask_key=gp_sources_categories_mask_key,
                      latent_key=latent_key,
                      conv_layer_encoder=conv_layer_encoder,
                      active_gp_thresh_ratio=active_gp_thresh_ratio,
                      seed = seed)
      
      np.random.seed(seed)
      import torch
      torch.set_num_threads(10)
      model.train(n_epochs=n_epochs,
                n_epochs_all_gps=n_epochs_all_gps,
                lr=lr,
                lambda_edge_recon=lambda_edge_recon,
                lambda_gene_expr_recon=lambda_gene_expr_recon,
                lambda_l1_masked=lambda_l1_masked,
                edge_batch_size=edge_batch_size,
                n_sampled_neighbors=n_sampled_neighbors,
                use_cuda_if_available=False,
                verbose=False)
      return model
    
    start_time = time.time()
    mem_usage = memory_usage((model_execution,), max_usage = True, retval=True)
    end_time = time.time()

    model = mem_usage[1]
    mem_usage_values = mem_usage[0]
    return mem_usage_values, end_time - start_time, model


def find_optimal_resolution(adata, target_clusters, key_added, neighbors_key, max_iter=500, tol=1):
    lower_res = 0.0
    upper_res = 5.0
    optimal_res = (lower_res + upper_res) / 2
    iter = 0
    
    while iter < max_iter:
        iter += 1
        sc.tl.leiden(adata, resolution=optimal_res, neighbors_key=neighbors_key, key_added = key_added)
        num_clusters = adata.obs['predicted_domain'].nunique()
        
        if abs(num_clusters - target_clusters) <= tol:
            return optimal_res
        elif num_clusters < target_clusters:
            lower_res = optimal_res
        else:
            upper_res = optimal_res
        
        optimal_res = (lower_res + upper_res) / 2
    
    return optimal_res

def convert_id(adata, from_id, to_id, species):
    mg = mygene.MyGeneInfo()
    original_ids = adata.var_names.tolist()
    query_result = mg.querymany(original_ids, scopes=from_id, fields=to_id, species=species)
    id_mapping = {}
    for item in query_result:
        if to_id in item:
            id_mapping[item['query']] = item[to_id]
        else:
            id_mapping[item['query']] = item['query']
            
    new_var_names = [id_mapping.get(id, id) for id in adata.var_names]
    adata.var[to_id] = new_var_names
    adata.var_names = adata.var[to_id]
    return adata


def run_nichecompass(input_file, output_file, sample, nclust, species,convertID,from_id,to_id, seed):
    # Define other constants and configurations here
    spatial_key = "spatial"
    n_neighbors = 4
    
    counts_key = "counts"
    adj_key = "spatial_connectivities"
    cat_covariates_keys = ["slices"]
    gp_names_key = "nichecompass_gp_names"
    active_gp_names_key = "nichecompass_active_gp_names"
    gp_targets_mask_key = "nichecompass_gp_targets"
    gp_targets_categories_mask_key = "nichecompass_gp_targets_categories"
    gp_sources_mask_key = "nichecompass_gp_sources"
    gp_sources_categories_mask_key = "nichecompass_gp_sources_categories"
    latent_key = "nichecompass_latent"
    conv_layer_encoder = "gcnconv"
    
    active_gp_thresh_ratio = 0.01
    n_epochs = 400
    n_epochs_all_gps = 25
    lr = 0.001
    lambda_edge_recon = 500000.
    lambda_gene_expr_recon = 300.
    lambda_l1_masked = 0.
    lambda_l1_addon = 100.
    edge_batch_size = 1024
    n_sampled_neighbors = 4
    use_cuda_if_available = True
    
    cell_type_key = "original_domain"
    latent_cluster_key = "predicted_domain"
    spot_size = 0.2
    differential_gp_test_results_key = "nichecompass_differential_gp_test_results"

    # Define paths
    ga_data_folder_path = "/home/dongkj/home_dkj/FD_yzy/Intergrated_method/NicheCompass/nichecompass-reproducibility-main/dkj_re/data/gene_annotations"
    gp_data_folder_path = "/home/dongkj/home_dkj/FD_yzy/Intergrated_method/NicheCompass/nichecompass-reproducibility-main/dkj_re/data/gene_programs"
    so_data_folder_path = "/home/dongkj/home_dkj/FD_yzy/Intergrated_method/NicheCompass/nichecompass-reproducibility-main/dkj_re/data/spatial_omics"
    omnipath_lr_network_file_path = f"{gp_data_folder_path}/omnipath_lr_network.csv"
    collectri_tf_network_file_path = f"{gp_data_folder_path}/collectri_tf_network_{species}.csv"
    nichenet_lr_network_file_path = f"{gp_data_folder_path}/nichenet_lr_network_v2_{species}.csv"
    nichenet_ligand_target_matrix_file_path = f"{gp_data_folder_path}/nichenet_ligand_target_matrix_v2_{species}.csv"
    mebocost_enzyme_sensor_interactions_folder_path = f"{gp_data_folder_path}/metabolite_enzyme_sensor_gps"
    gene_orthologs_mapping_file_path = f"{ga_data_folder_path}/human_mouse_gene_orthologs.csv"

    # Retrieve OmniPath GPs (source: ligand genes; target: receptor genes)
    omnipath_gp_dict = extract_gp_dict_from_omnipath_lr_interactions(
        species=species,
        min_curation_effort=0,
        load_from_disk=True,
        save_to_disk=False,
        lr_network_file_path=omnipath_lr_network_file_path,
        gene_orthologs_mapping_file_path=gene_orthologs_mapping_file_path,
        plot_gp_gene_count_distributions=False)
        
    # Retrieve MEBOCOST GPs (source: enzyme genes; target: sensor genes)
    mebocost_gp_dict = extract_gp_dict_from_mebocost_es_interactions(
        dir_path=mebocost_enzyme_sensor_interactions_folder_path,
        species=species,
        plot_gp_gene_count_distributions=False)


    # Retrieve NicheNet GPs (source: ligand genes; target: receptor genes, target genes)
    nichenet_gp_dict = extract_gp_dict_from_nichenet_lrt_interactions(
        species=species,
        version="v2",
        keep_target_genes_ratio=1.,
        max_n_target_genes_per_gp=250,
        load_from_disk=True,
        save_to_disk=False,
        lr_network_file_path=nichenet_lr_network_file_path,
        ligand_target_matrix_file_path=nichenet_ligand_target_matrix_file_path,
        gene_orthologs_mapping_file_path=gene_orthologs_mapping_file_path,
        plot_gp_gene_count_distributions=False)


    # Add GPs into one combined dictionary for model training
    combined_gp_dict = dict(omnipath_gp_dict)
    combined_gp_dict.update(mebocost_gp_dict)
    combined_gp_dict.update(nichenet_gp_dict)

    # Filter and combine GPs to avoid overlaps
    combined_new_gp_dict = filter_and_combine_gp_dict_gps(
        gp_dict=combined_gp_dict,
        gp_filter_mode="subset",
        combine_overlap_gps=True,
        overlap_thresh_source_genes=0.9,
        overlap_thresh_target_genes=0.9,
        overlap_thresh_genes=0.9)

    # Load data
    batches = os.listdir(input_file)
    batches = [batch.split('.')[0] for batch in batches]
    cat_covariates_embeds_nums = [np.array(batches).shape[0]]
    adata_batch_list = []
    for batch in batches:
        print(f"Processing batch {batch}...")
        print("Loading data...")
        adata_batch = sc.read_h5ad(f"{input_file}/{batch}.h5ad")
        print("Computing spatial neighborhood graph...\n")
        # Compute (separate) spatial neighborhood graphs
        sq.gr.spatial_neighbors(adata_batch,
                                coord_type="generic",
                                spatial_key=spatial_key,
                                n_neighs=n_neighbors)
        
        # Make adjacency matrix symmetric
        adata_batch.obsp[adj_key] = (
            adata_batch.obsp[adj_key].maximum(
                adata_batch.obsp[adj_key].T))
        adata_batch_list.append(adata_batch)
    
    adata = ad.concat(adata_batch_list, join="inner")

    # Combine spatial neighborhood graphs as disconnected components
    batch_connectivities = []
    len_before_batch = 0
    for i in range(len(adata_batch_list)):
        if i == 0: # first batch
            after_batch_connectivities_extension = sp.csr_matrix(
                (adata_batch_list[0].shape[0],
                (adata.shape[0] -
                adata_batch_list[0].shape[0])))
            batch_connectivities.append(sp.hstack(
                (adata_batch_list[0].obsp[adj_key],
                after_batch_connectivities_extension)))
        elif i == (len(adata_batch_list) - 1): # last batch
            before_batch_connectivities_extension = sp.csr_matrix(
                (adata_batch_list[i].shape[0],
                (adata.shape[0] -
                adata_batch_list[i].shape[0])))
            batch_connectivities.append(sp.hstack(
                (before_batch_connectivities_extension,
                adata_batch_list[i].obsp[adj_key])))
        else: # middle batches
            before_batch_connectivities_extension = sp.csr_matrix(
                (adata_batch_list[i].shape[0], len_before_batch))
            after_batch_connectivities_extension = sp.csr_matrix(
                (adata_batch_list[i].shape[0],
                (adata.shape[0] -
                adata_batch_list[i].shape[0] -
                len_before_batch)))
            batch_connectivities.append(sp.hstack(
                (before_batch_connectivities_extension,
                adata_batch_list[i].obsp[adj_key],
                after_batch_connectivities_extension)))
        len_before_batch += adata_batch_list[i].shape[0]

    adata.obsp[adj_key] = sp.vstack(batch_connectivities)

    if convertID:
      adata = convert_id(adata, from_id, to_id, species)

    # Add the GP dictionary as binary masks to the adata
    add_gps_from_gp_dict_to_adata(
        gp_dict=combined_new_gp_dict,
        adata=adata,
        gp_targets_mask_key=gp_targets_mask_key,
        gp_targets_categories_mask_key=gp_targets_categories_mask_key,
        gp_sources_mask_key=gp_sources_mask_key,
        gp_sources_categories_mask_key=gp_sources_categories_mask_key,
        gp_names_key=gp_names_key,
        min_genes_per_gp=2,
        min_source_genes_per_gp=1,
        min_target_genes_per_gp=1,
        max_genes_per_gp=None,
        max_source_genes_per_gp=None,
        max_target_genes_per_gp=None)
        
    adata.layers[counts_key] = adata.X.copy()
    count_data = adata.X.copy()
    
    memory_usage_value, execution_time, model = train_model(
        adata, counts_key, adj_key, cat_covariates_keys, gp_names_key, active_gp_names_key,
        gp_targets_mask_key, gp_targets_categories_mask_key, gp_sources_mask_key, gp_sources_categories_mask_key,
        latent_key, conv_layer_encoder, active_gp_thresh_ratio, n_epochs, n_epochs_all_gps, lr, lambda_edge_recon,
        lambda_gene_expr_recon, lambda_l1_masked, lambda_l1_addon, edge_batch_size, n_sampled_neighbors, use_cuda_if_available,cat_covariates_embeds_nums, seed)
    
    # Compute latent neighbor graph
    sc.pp.neighbors(model.adata,
                    use_rep=latent_key,
                    key_added=latent_key)

    # Compute UMAP embedding
    sc.tl.umap(model.adata,
               neighbors_key=latent_key)

    # Perform analysis
    samples = model.adata.obs["slices"].unique().tolist()

    # Compute latent Leiden clustering
    optimal_resolution = find_optimal_resolution(model.adata, nclust, latent_cluster_key, latent_key)
    sc.tl.leiden(adata=model.adata,
                 resolution=optimal_resolution,
                 key_added=latent_cluster_key,
                 neighbors_key=latent_key)

    # Write output to file
    model.adata.obsm["X_embedding"] = model.adata.obsm["nichecompass_latent"]
    model.adata.write(f"{output_file}/{sample}.h5ad")
    
    # Save the statistics to a CSV file
    stats_file = os.path.join(output_file, f"{sample}_stats.csv")
    with open(stats_file, "w", newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Memory_Usage_MiB", "Execution_Time_s"])
        csvwriter.writerow([memory_usage_value, execution_time])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run NicheCompass")
    parser.add_argument("--input_file", type=str, help="Input file path")
    parser.add_argument("--output_file", type=str, help="Output file path")
    parser.add_argument("--sample", type=str, default=None, help="Sample size for downsampling")
    parser.add_argument("--nclust", type=int, help="Cluster number")
    parser.add_argument("--species", type=str, default="human", help="Species name (human/mouse)")
    parser.add_argument('--convertID', type=bool, default=False, help='Need convert gene id?')
    parser.add_argument('--from_id', type=str, default = "ensembl.gene", help='species of the data')
    parser.add_argument('--to_id', type=str, default = "symbol", help='species of the data')
    parser.add_argument('--seed', type=int, default = 0, help='seed')
    args = parser.parse_args()
    
    run_nichecompass(args.input_file, args.output_file, args.sample, args.nclust, args.species, args.convertID, args.from_id, args.to_id, args.seed)
