# Import necessary libraries for spatial transcriptomics data analysis and metric computation
import scanpy as sc
import scib
import os
import anndata as ad
import pandas as pd
import argparse
import numpy as np
import random
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection 

# Main function to calculate integration metrics
def calculate_metrics(input_file, input_path, output_path, seed = 123):
    # Set random seed for reproducibility
    np.random.seed(seed)
    random.seed(seed)
    
    # Get a list of .h5ad files in the input directory, excluding those containing 'withRaw'
    fl = [f for f in os.listdir(input_file) if f.endswith('.h5ad') and 'withRaw' not in f]

    results = []

    # Iterate over each integrated data file
    for f in fl:
        model = os.path.splitext(f)[0]  # Get model name from file name
        inter_data = ad.read_h5ad(os.path.join(input_file, f))  # Load integration result as AnnData

        # Ensure correct categorical format for batch and biological label
        inter_data.obs["slices"] = inter_data.obs["slices"].astype(str).astype('category')
        inter_data.obs["original_domain"] = inter_data.obs["original_domain"].astype(str).astype('category')

        # Choose embedding key depending on the model name
        if model == "GraphST" or model == "GraphSTwithPASTE":
            embed = "emb"
        else:
            embed = "X_embedding"

        # Initialize the Benchmarker with selected metrics
        bm = Benchmarker(
            inter_data,
            batch_key="slices",  # Key indicating batch (e.g., slices)
            label_key="original_domain",  # Key indicating biological identity
            bio_conservation_metrics=BioConservation(
                isolated_labels=True,
                nmi_ari_cluster_labels_leiden=False,
                nmi_ari_cluster_labels_kmeans=False,
                silhouette_label=True,
                clisi_knn=True),
            batch_correction_metrics=BatchCorrection(
                silhouette_batch=True,
                ilisi_knn=True,
                kbet_per_label=False,
                graph_connectivity=True,
                pcr_comparison=False),
            embedding_obsm_keys=[embed],  # Specify which embedding to evaluate
            n_jobs=10                     # Parallel processing
        )

        # Run the benchmark
        bm.benchmark()

        # Retrieve raw metric scores (not min-max scaled)
        df = bm.get_results(min_max_scale=False)

        # Extract selected metrics from the result DataFrame
        dASW = df.loc[embed, 'Silhouette label']         # Domain ASW
        dLISI = df.loc[embed, 'cLISI']                   # Domain LISI
        ILL = df.loc[embed, 'Isolated labels']           # Isolated label F1

        bASW = df.loc[embed, 'Silhouette batch']         # Batch ASW
        iLISI = df.loc[embed, 'iLISI']                   # Batch LISI
        GC = df.loc[embed, 'Graph connectivity']         # Graph connectivity

        # Append metric results for the current model
        results.append([model, dASW, dLISI, ILL, bASW, iLISI, GC])

    # Convert results to DataFrame and save as CSV
    results_df = pd.DataFrame(results, columns=["Model", "dASW", "dLISI", "ILL", "bASW", "iLISI", "GC"])
    results_df.to_csv(output_path, index=False)

# Entry point for command-line execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate integration metrics for spatial transcriptomics data.")
    
    # Define required input arguments
    parser.add_argument('--input_file', type=str, required=True, 
                        help="Directory containing the integration result files in .h5ad format.")
    parser.add_argument('--input_path', type=str, required=True, 
                        help="Path to the original combined slices data in .h5ad format.")
    parser.add_argument('--output_path', type=str, required=True, 
                        help="Path where the output CSV file with metrics will be saved.")
    
    # Parse arguments and run the metric calculation
    args = parser.parse_args()
    calculate_metrics(args.input_file, args.input_path, args.output_path)
