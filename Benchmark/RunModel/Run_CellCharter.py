import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import cellcharter as cc
import scvi
import squidpy as sq
import os
import time
from memory_profiler import memory_usage
import csv
import argparse
import torch
import scarches as sca

print("CUDA_VISIBLE_DEVICES:", os.getenv('CUDA_VISIBLE_DEVICES'))
print("CUDA_HOME:", os.getenv('CUDA_HOME'))


def run_model(adata, seeds, n_latent, nhood_layers, nclust, isCODEX):
    def model_execution():
        for seed in seeds:
            if isCODEX:
              source_conditions = adata.obs["slices"].unique().tolist()
              trvae = sca.models.TRVAE(
                adata=adata,
                condition_key="slices",
                conditions=source_conditions
                )
              early_stopping_kwargs = {
                "early_stopping_metric": "val_unweighted_loss",
                "threshold": 0,
                "patience": 20,
                "reduce_lr": True,
                "lr_patience": 13,
                "lr_factor": 0.1,
                }
              trvae.train(early_stopping_kwargs=early_stopping_kwargs)
              adata.obsm['X_scVI'] = trvae.get_latent()
            else:
              scvi.settings.seed = seed
              scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="slices")
              model = scvi.model.SCVI(adata, n_latent=n_latent)
              model.train(early_stopping=True)
              adata.obsm['X_scVI'] = model.get_latent_representation(adata).astype(np.float32)
  
            sq.gr.spatial_neighbors(adata, library_key='slices', coord_type='generic', delaunay=True, spatial_key='spatial')
            cc.gr.remove_long_links(adata)
            cc.gr.aggregate_neighbors(adata, nhood_layers, use_rep='X_scVI', out_key='X_cellcharter', sample_key='slices')
            cls = cc.tl.Cluster(nclust, random_state=seed)
            cls.fit(adata, use_rep='X_cellcharter')
            adata.obs["predicted_domain"] = cls.predict(adata, use_rep='X_cellcharter')
        adata.obsm["X_embedding"] = adata.obsm["X_cellcharter"]
        return adata
      
    start_time = time.time()
    mem_usage = memory_usage((model_execution,), max_usage = True, retval=True)
    end_time = time.time()
    adata = mem_usage[1]
    mem_usage_values = mem_usage[0] 
    return mem_usage_values, end_time - start_time, adata

def process_data(input_file, output_file, sample, nclust, hvgs, runNormalization, isCODEX, n_latent, nhood_layers, seed):

    adata = sc.read_h5ad(input_file)
    adata.obs["slices"] = adata.obs["slices"].astype(str)
    adata.obs["slices"] = adata.obs["slices"].astype('category')
    adata.obs["original_domain"] = adata.obs["original_domain"].astype(str)
    count_data = adata.X.copy()
    adata.layers["counts"] = adata.X.copy()

    if runNormalization:
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=hvgs,
            subset=True,
            layer="counts",
            flavor="seurat_v3"
        )

    adata.var_names_make_unique()

    #rng = np.random.default_rng(12345)
    seeds = np.array([seed])

    torch.set_num_threads(10)
    memory_stats, time_stats,adata = run_model(adata, seeds, n_latent, nhood_layers, nclust, isCODEX)

    # Save the statistics to a CSV file
    stats_file = os.path.join(output_file, f"{sample}_stats.csv")
    with open(stats_file, "w", newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Memory_Usage_MiB", "Execution_Time_s"])
        csvwriter.writerow([memory_stats, time_stats])
    adata.write(os.path.join(output_file, f"{sample}.h5ad"))
    print("success!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--input_file', type=str, help='input file path')
    parser.add_argument('--output_file', type=str, help='output file path')
    parser.add_argument('--sample', type=str, help='sample name you analysis')
    parser.add_argument('--nclust', type=int, help='number of clusters')
    parser.add_argument('--hvgs', type=int, default=5000, help='number of highly variable genes (default: 5000)')
    parser.add_argument('--runNormalization', type=bool, default=True, help='Run data normalization before run integration?')
    parser.add_argument('--isCODEX', type=bool, default=False, help='Is the data is CODEX?')
    parser.add_argument('--n_latent', type=int, help='Hyperparameter: n_latent')
    parser.add_argument('--nhood_layers', type=int, help='Hyperparameter: nhood_layers')
    parser.add_argument('--seed', type=int, default=1234, help='seed')
    args = parser.parse_args()


    process_data(args.input_file, args.output_file, args.sample, args.nclust, args.hvgs, args.runNormalization, args.isCODEX, args.n_latent, args.nhood_layers, args.seed)
