# Import necessary libraries and modules
import os
import argparse
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp
import matplotlib.pyplot as plt
from GraphST import GraphST
import numpy as np
from GraphST.utils import clustering
import time
from memory_profiler import memory_usage
import csv


# Function to process the input data and perform normalization (if specified)
# This function reads input data, normalizes it if required, selects highly variable genes, and calls the model execution
def process_data(input_file, output_file, sample, nclust, hvgs, run_normalization, tool, device, seed):
    # input_file: Path to the multi-slice data that needs to be integrated. All slices should be contained in a single h5ad file. Refer to the example data for the specific format. This is a required parameter.
    # output_file: Path to the directory where the output results will be saved. This is a required parameter.
    # sample: Name for the output result files. It is recommended to use the model name, for example, setting it as "GraphST" will result in output files named "GraphST.h5ad". This is a required parameter.
    # nclust: Number of identified domains. This is a required parameter.
    # hvgs: The number of highly variable genes identified if data normalization is performed. This is an optional parameter, with the default set to 3000.
    # runNormalization: A boolean variable indicating whether data normalization should be performed before integration. This is an optional parameter, with the default set to True.
    # tool: Methods for spatial clustering, which can be set to "mclust", "leiden" and "louvain". This is an optional parameter, with the default set to "mclust".
    # device: The device on which the model runs can be set to "cpu" and "gpu". This is an optional parameter, with the default set to "cpu".
    # seed: Random seed for model execution. This is an optional parameter, with the default set to 0.


    # Set the environment variable for R (required for Scanpy)
    # Set R based on your environment
    # os.environ['R_HOME'] = '/opt/R/4.3.2/lib/R'

    # Read the input data file into an AnnData object
    adata = sc.read_h5ad(input_file)
    adata.var_names_make_unique()
    adata.obs['data'] = adata.obs['slices'].astype(str) 
    count_data = adata.X.copy() 

    # Perform normalization if specified
    if run_normalization:
        total_counts_per_cell = np.sum(adata.X, axis=1) 
        sc.pp.log1p(adata) 
        
        sc.pp.normalize_total(adata, target_sum=total_counts_per_cell)
        sc.pp.scale(adata, max_value=10)  
        sc.pp.highly_variable_genes(adata, n_top_genes=hvgs)
        
        highly_variable_data = adata[:, adata.var['highly_variable']] 
        highly_variable_data.obs_names = highly_variable_data.obs_names + '_' + highly_variable_data.obs['slices'].astype(str)
    else:
        # If no normalization, use the original data
        highly_variable_data = adata
        highly_variable_data.obs_names = highly_variable_data.obs_names + '_' + highly_variable_data.obs['slices'].astype(str)

    # Set the number of threads for PyTorch
    torch.set_num_threads(10)
    
    # Run the model and get the memory usage and execution time
    memory_usage_value, execution_time, adata = run_model(highly_variable_data, device, nclust, tool, seed)

    # Save performance statistics (memory usage and execution time) to a CSV file
    stats_file = os.path.join(output_file, f"{sample}_stats.csv")
    with open(stats_file, "w", newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Memory_Usage_MiB", "Execution_Time_s"])
        csvwriter.writerow([memory_usage_value, execution_time])
    
    # Store the embedding and predicted domain in the AnnData object
    adata.obsm["X_embedding"] = adata.obsm["emb"]
    adata.obs['predicted_domain'] = adata.obs['domain']

    # Save the final results to an h5ad file
    adata.write(f"{output_file}/{sample}.h5ad")


# Function to run the model and return performance metrics
# This function initializes the GraphST model, trains it, and performs clustering
# It also tracks memory usage and execution time during the training process
def run_model(highly_variable_data, device, nclust, tool, seed):
    # Parameters:
    # highly_variable_data: The input data with highly variable genes selected
    # device: The computing device to be used (e.g., 'cpu' or 'gpu')
    # nclust: The number of clusters for clustering algorithms
    # tool: The clustering method to be used (e.g., 'mclust', 'leiden', 'louvain')
    # seed: Random seed for reproducibility
    np.random.seed(seed)  # Set random seed for reproducibility

    def model_execution():
        # Initialize and train the GraphST model
        model = GraphST.GraphST(highly_variable_data, device=device, random_seed=seed)
        adata = model.train()

        # Perform clustering using the specified clustering tool
        if tool == 'mclust':
            clustering(adata, nclust, method=tool, refinement=True)
        elif tool in ['leiden', 'louvain']:
            clustering(adata, nclust, method=tool, start=0.1, end=2.0, increment=0.01)

        return adata

    # Measure execution time and memory usage during model execution
    start_time = time.time()  # Start time for tracking execution duration
    mem_usage = memory_usage((model_execution,), max_usage=True, retval=True)  # Track memory usage
    end_time = time.time()  # End time after model execution

    adata = mem_usage[1]  # The resulting annotated data object
    mem_usage_values = mem_usage[0]  # Memory usage values during execution
    return mem_usage_values, end_time - start_time, adata  # Return memory usage, execution time, and final data object

# Main function to parse arguments and run the data processing pipeline
# This function parses input arguments and calls the process_data function
if __name__ == "__main__":
    # Set up argument parser to accept command-line inputs
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--input_file', type=str, help='input file path')
    parser.add_argument('--output_file', type=str, help='output file path')
    parser.add_argument('--sample', type=str, help='sample name')
    parser.add_argument('--nclust', type=int, help='number of clusters')
    parser.add_argument('--hvgs', type=int, default=3000, help='number of highly variable genes')
    parser.add_argument('--runNormalization', type=str, default='True', help='whether to run normalization')
    parser.add_argument('--tool', type=str, default="mclust", help='clustering tool')
    parser.add_argument('--device', type=str, default="cpu", help='gpu or cpu')
    parser.add_argument('--seed', type=int, default=41, help='random seed')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Convert normalization flag to boolean
    run_normalization = args.runNormalization.lower() in ('true', 't', 'yes', 'y', '1')

    # Call the process_data function with the parsed arguments
    process_data(args.input_file, args.output_file, args.sample, args.nclust, args.hvgs, run_normalization, args.tool, args.device, args.seed)
