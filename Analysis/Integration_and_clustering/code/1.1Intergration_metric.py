# Import necessary libraries for spatial transcriptomics data analysis and metric computation
import scanpy as sc
import scib
import os
import anndata as ad
import pandas as pd
import argparse
import numpy as np
import random

# Set environment variable for R to integrate with Python
os.environ['R_HOME'] = '/home/dongkj/anaconda3/envs/MultiSpatial/lib/R'
# os.environ['R_HOME'] = '/opt/R/4.3.2/lib/R'

# Function to calculate various metrics for spatial transcriptomics integration results
def calculate_metrics(input_file, input_path, output_path, seed = 123):
   # Parameters:
   # input_file: Path to directory containing .h5ad files with integration results
   # input_path: Path to the original combined slices data in .h5ad format
   # output_path: Path where the resulting CSV file with integration metrics will be saved
   # seed: Random seed for reproducibility of results (default: 123)
  
    # Set the random seed for numpy and random to ensure reproducibility of results
    np.random.seed(seed)
    random.seed(seed)
    
    # List all .h5ad files in the input directory that do not include 'withRaw' in their names
    fl = [f for f in os.listdir(input_file) if f.endswith('.h5ad') and 'withRaw' not in f]

    # Initialize a list to store the results for each model
    results = []
    
    # Loop over each file in the directory to compute metrics
    for f in fl:
        model = os.path.splitext(f)[0]  # Extract the model name from the file name
        # Read in the integration results and the original data (AnnData format)
        inter_data = ad.read_h5ad(os.path.join(input_file, f))
        inter_data.obs["slices"] = inter_data.obs["slices"].astype(str).astype('category')
        inter_data.obs["original_domain"] = inter_data.obs["original_domain"].astype(str).astype('category')

        # Load the original combined slices data
        original_data = ad.read_h5ad(input_path)
        original_data.obs["slices"] = original_data.obs["slices"].astype(str).astype('category')
        original_data.obs["original_domain"] = original_data.obs["original_domain"].astype(str).astype('category')

        # Determine which embedding to use based on the model type
        if model == "GraphST" or model == "GraphSTwithPASTE":
            embed = "emb"  # Use "emb" embedding for these specific models
        else:
            embed = "X_embedding"  # Default to "X_embedding" for other models

        # Calculate silhouette scores for domain and batch consistency
        # dASW: Domain silhouette score (measures clustering quality within a domain)
        # bASW: Batch silhouette score (measures clustering quality within batches)
        dASW = scib.me.silhouette(inter_data, label_key="original_domain", embed = embed, scale = True)
        bASW = scib.me.silhouette_batch(inter_data, batch_key="slices", label_key="original_domain", embed=embed, scale = True)
        
        # Calculate kNN-based metrics for domain and batch consistency
        # ilF1: Isolated labels F1 score (measures how well labels are separated within batches)
        # kBET: kBET score (measures the preservation of batch structure across domains)
        sc.pp.neighbors(inter_data, use_rep=embed)  # Compute nearest neighbors using the selected embedding
        ilF1  = scib.me.isolated_labels_f1(inter_data, batch_key="slices", label_key="original_domain", embed = None)
        kBET = scib.me.kBET(inter_data, batch_key="slices", label_key="original_domain", type_="knn", scaled = True)

        # Append the results for the current model to the results list
        results.append([model, dASW, bASW, ilF1, kBET])

    # Convert the results list into a pandas DataFrame for easy manipulation
    results_df = pd.DataFrame(results, columns=["Model", "dASW", "bASW", "ilF1", "kBET"])
    
    # Save the results DataFrame to a CSV file at the specified output path
    results_df.to_csv(output_path, index=False)

# Main function to parse command-line arguments and call the metrics calculation function
if __name__ == "__main__":
    # Set up argument parser to collect input and output paths from the command line
    parser = argparse.ArgumentParser(description="Calculate integration metrics for spatial transcriptomics data.")
    
    # Define the command-line arguments
    parser.add_argument('--input_file', type=str, required=True, 
                        help="Directory containing the integration result files in .h5ad format.")
    parser.add_argument('--input_path', type=str, required=True, 
                        help="Path to the original combined slices data in .h5ad format.")
    parser.add_argument('--output_path', type=str, required=True, 
                        help="Path where the output CSV file with metrics will be saved.")
    
    # Parse the command-line arguments
    args = parser.parse_args()
    
    # Call the calculate_metrics function with the parsed arguments
    calculate_metrics(args.input_file, args.input_path, args.output_path)
