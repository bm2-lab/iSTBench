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

# Display CUDA environment variables for debugging purposes
print("CUDA_VISIBLE_DEVICES:", os.getenv('CUDA_VISIBLE_DEVICES'))
print("CUDA_HOME:", os.getenv('CUDA_HOME'))

def process_data(input_file, output_file, sample, nclust, hvgs, runNormalization, isCODEX, n_latent, nhood_layers, seed):
    """
    Preprocesses the input data and runs the model to predict domains using CellCharter or SCVI.

    Parameters:
    input_file: Path to the multi-slice data that needs to be integrated. All slices should be contained in a single h5ad file. Refer to the example data for the specific format. This is a required parameter.
    output_file: Path to the directory where the output results will be saved. This is a required parameter.
    sample: Name for the output result files. It is recommended to use the model name, for example, setting it as "CellCharter" will result in output files named "CellCharter.h5ad". This is a required parameter.
    nclust: Number of identified domains. This is a required parameter.
    hvgs: The number of highly variable genes identified if data normalization is performed. This is an optional parameter, with the default set to 5000.
    runNormalization: A boolean variable indicating whether data normalization should be performed before integration. This is an optional parameter, with the default set to True.
    isCODEX: A boolean variable indicates whether the data were obtained by CODEX sequencing. This is an optional parameter, with the default set to False.
    n_latent: Hyperparameter for the CellCharter model. This is a required parameter. The original recommendation is to conduct a grid search in the hyperparameter space [5, 10, 15] to select the best hyperparameter.
    nhood_layers: Hyperparameter for the CellCharter model. This is a required parameter. The original recommendation is to conduct a grid search in the hyperparameter space [2, 3, 4] to select the best hyperparameter.
    seed: Random seed for model execution. This is an optional parameter, with the default set to 1234.

    Returns:
    None
    """
    
    # Load the data from an h5ad file
    adata = sc.read_h5ad(input_file)
    adata.obs["slices"] = adata.obs["slices"].astype(str)
    adata.obs["slices"] = adata.obs["slices"].astype('category')
    adata.obs["original_domain"] = adata.obs["original_domain"].astype(str)
    count_data = adata.X.copy()  # Preserve the raw count data
    adata.layers["counts"] = adata.X.copy()  # Store the counts in layers

    # Normalize the data if the 'runNormalization' flag is True
    if runNormalization:
        sc.pp.normalize_total(adata, target_sum=1e6)  # Normalize to a total count of 1e6
        sc.pp.log1p(adata)  # Log transform the data
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=hvgs,  # Select top 'hvgs' highly variable genes
            subset=True,
            layer="counts",  # Apply this to the 'counts' layer
            flavor="seurat_v3"  # Seurat v3 method for HVG selection
        )

    adata.var_names_make_unique()  # Ensure unique gene names

    # Convert the random seed to a numpy array for reproducibility
    seeds = np.array([seed])

    # Set the number of threads for PyTorch
    torch.set_num_threads(10)

    # Run the model and get memory usage, time statistics, and the updated adata object
    memory_stats, time_stats, adata = run_model(adata, seeds, n_latent, nhood_layers, nclust, isCODEX)

    # Save the memory usage and execution time statistics to a CSV file
    stats_file = os.path.join(output_file, f"{sample}_stats.csv")
    with open(stats_file, "w", newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Memory_Usage_MiB", "Execution_Time_s"])
        csvwriter.writerow([memory_stats, time_stats])

    # Save the processed AnnData object to a file
    adata.write(os.path.join(output_file, f"{sample}.h5ad"))
    print("success!")


def run_model(adata, seeds, n_latent, nhood_layers, nclust, isCODEX):
    """
    Run the CellCharter model or SCVI model, depending on whether the data is from CODEX or not.

    Parameters:
    adata: AnnData object containing the dataset to process
    seeds: List of random seeds to ensure reproducibility
    n_latent: Number of latent variables for the model
    nhood_layers: Number of neighborhood layers for the CellCharter model
    nclust: Number of clusters (domains) to be identified
    isCODEX: Boolean flag indicating whether the data is from CODEX sequencing

    Returns:
    mem_usage_values: Memory usage statistics during the model execution
    time_stats: Total execution time for the model
    adata: The processed AnnData object with the predicted domains and latent representations
    """
    
    def model_execution():
        """
        Executes the model training and clustering for each seed. Either CellCharter or SCVI is used 
        depending on whether the dataset is from CODEX.
        """
        for seed in seeds:
            if isCODEX:
                # Set up and train the CellCharter TRVAE model if the data is from CODEX
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
                # Set up and train the SCVI model if the data is not from CODEX
                scvi.settings.seed = seed
                scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="slices")
                model = scvi.model.SCVI(adata, n_latent=n_latent)
                model.train(early_stopping=True)
                adata.obsm['X_scVI'] = model.get_latent_representation(adata).astype(np.float32)
  
            # Calculate spatial neighbors for downstream analysis
            sq.gr.spatial_neighbors(adata, library_key='slices', coord_type='generic', delaunay=True, spatial_key='spatial')
            cc.gr.remove_long_links(adata)  # Clean long-range links that may affect clustering
            cc.gr.aggregate_neighbors(adata, nhood_layers, use_rep='X_scVI', out_key='X_cellcharter', sample_key='slices')
            
            # Cluster the data and assign domain labels
            cls = cc.tl.Cluster(nclust, random_state=seed)
            cls.fit(adata, use_rep='X_cellcharter')
            adata.obs["predicted_domain"] = cls.predict(adata, use_rep='X_cellcharter')
        
        # Store the final embedding in the adata object
        adata.obsm["X_embedding"] = adata.obsm["X_cellcharter"]
        return adata
      
    # Measure memory usage and execution time during the model training and clustering
    start_time = time.time()
    mem_usage = memory_usage((model_execution,), max_usage=True, retval=True)
    end_time = time.time()
    
    # Return the memory usage, execution time, and the processed AnnData object
    adata = mem_usage[1]
    mem_usage_values = mem_usage[0] 
    return mem_usage_values, end_time - start_time, adata


# Main execution block to handle input arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--input_file', type=str, help='Input file path (h5ad format)')
    parser.add_argument('--output_file', type=str, help='Output directory path')
    parser.add_argument('--sample', type=str, help='Sample name for output files')
    parser.add_argument('--nclust', type=int, help='Number of clusters (domains)')
    parser.add_argument('--hvgs', type=int, default=5000, help='Number of highly variable genes (default: 5000)')
    parser.add_argument('--runNormalization', type=bool, default=True, help='Run data normalization before integration?')
    parser.add_argument('--isCODEX', type=bool, default=False, help='Is the data from CODEX sequencing?')
    parser.add_argument('--n_latent', type=int, help='Number of latent dimensions for the model')
    parser.add_argument('--nhood_layers', type=int, help='Number of neighborhood layers for the model')
    parser.add_argument('--seed', type=int, default=1234, help='Random seed for model execution')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the process_data function with the parsed arguments
    process_data(args.input_file, args.output_file, args.sample, args.nclust, args.hvgs, args.runNormalization, args.isCODEX, args.n_latent, args.nhood_layers, args.seed)
