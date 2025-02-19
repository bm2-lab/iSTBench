import warnings
warnings.filterwarnings("ignore")
import MENDER
import scanpy as sc
import os
import time
from memory_profiler import memory_usage
import csv
import argparse
import pandas as pd
import numpy as np

def run_model(adata, n_cls, tech, scale, radius, seed):
    """
    Runs the MENDER model for spatial transcriptomics data integration and clustering.

    Parameters:
    input_file: Path to the multi-slice data that needs to be integrated. All slices should be contained in a single h5ad file. Refer to the example data for the specific format. This is a required parameter.
    output_file: Path to the directory where the output results will be saved. This is a required parameter.
    sample: Name for the output result files. It is recommended to use the model name, for example, setting it as "MENDER" will result in output files named "MENDER.h5ad". This is a required parameter.
    nclust: Number of identified domains. This is a required parameter.
    hvgs: The number of highly variable genes identified if data normalization is performed. This is an optional parameter, with the default set to 4000.
    runNormalization: A boolean variable indicating whether data normalization should be performed before integration. This is an optional parameter, with the default set to True.
    tech: Sequencing technology that generated the data. Can be set to "ST", "Visium" and others. This is a required parameter.
    scale: Hyperparameter for the MENDER model. This is an optional parameter, with the default set to 6.
    radius: Hyperparameter for the MENDER model. This is an optional parameter, with the default set to 15.
    n_pcs: The number of principal components set when performing PCA dimensionality reduction on the data. This is an optional parameter, with the default set to 50.
    seed: Random seed for model execution. This is an optional parameter, with the default set to 666.

    Returns:
    mem_usage_values: Memory usage statistics during model execution
    time_stats: Total execution time for the model
    adata: The processed AnnData object with the predicted domains and latent representations
    """
    
    np.random.seed(seed)

    def model_execution():
        """
        Executes the MENDER model for integration and clustering based on the sequencing technology.
        This includes setting parameters and performing clustering.

        Returns:
        adata: The updated AnnData object after model execution
        """
        # Batch and ground truth observation names
        batch_obs = 'slices'
        gt_obs = 'original_domain'
        
        # Initialize the MENDER model
        msm = MENDER.MENDER(
            adata,
            batch_obs='slices',
            ct_obs='ct',
            random_seed=seed
        )
        msm.prepare()  # Prepare the model with the dataset
        
        # Set parameters based on sequencing technology
        if tech == "ST":
            msm.set_MENDER_para(
                n_scales=scale,
                nn_mode='ring',
                nn_para=4
            )
        elif tech == "Visium":
            msm.set_MENDER_para(
                n_scales=scale,
                nn_mode='ring',
                nn_para=6
            )
        else:
            msm.set_MENDER_para(
                n_scales=scale,
                nn_mode='radius',
                nn_para=radius
            )
        
        # Run model representation and clustering
        msm.run_representation_mp(2)
        msm.run_clustering_normal(n_cls)
        
        # Reorder the embeddings and predicted domains based on barcodes
        barcodes_adata = adata.obs['barcode'].values
        barcodes_mender = msm.adata_MENDER.obs['barcode'].values
        adata_df = pd.DataFrame(index=barcodes_adata)
        
        # Reorder the embeddings according to the barcodes
        mender_X_df = pd.DataFrame(msm.adata_MENDER.X, index=barcodes_mender)
        mender_X_reordered = mender_X_df.reindex(barcodes_adata).values
        adata.obsm["X_embedding"] = mender_X_reordered
        
        # Reorder and assign predicted domains
        mender_obs_df = pd.DataFrame(msm.adata_MENDER.obs[['MENDER']].values, index=barcodes_mender, columns=['MENDER'])
        mender_obs_reordered = mender_obs_df.reindex(barcodes_adata)['MENDER'].values
        adata.obs["predicted_domain"] = mender_obs_reordered
        
        return adata

    # Measure execution time and memory usage
    start_time = time.time()
    mem_usage = memory_usage((model_execution,), max_usage=True, retval=True)
    end_time = time.time()

    adata = mem_usage[1]  # Get the final AnnData object after model execution
    mem_usage_values = mem_usage[0]  # Get memory usage during model execution
    
    return mem_usage_values, end_time - start_time, adata


def process_data(input_file, output_file, sample, n_cls, hvgs, runNormalization, tech, scale, radius, n_pcs, seed):
    """
    Preprocesses the input data, normalizes it (if specified), and runs the MENDER model for integration and clustering.

    Parameters:
    input_file: Path to the multi-slice data that needs to be integrated. All slices should be contained in a single h5ad file.
    output_file: Path to the directory where output results will be saved
    sample: Name for the output result files
    n_cls: Number of clusters (domains) to identify
    hvgs: The number of highly variable genes identified if data normalization is performed
    runNormalization: Whether to normalize the data (True/False)
    tech: Sequencing technology used ('ST', 'Visium', etc.)
    scale: Scale parameter for the MENDER model
    radius: Radius parameter for the MENDER model
    n_pcs: The number of principal components to use for PCA
    seed: Random seed for model execution

    Returns:
    None
    """
    
    # Load the data from an h5ad file
    adata_raw = sc.read_h5ad(input_file)
    adata_raw.obs['slices'] = adata_raw.obs['slices'].astype('category')
    adata_raw.obs['original_domain'] = adata_raw.obs['original_domain'].astype('category')
    batch_obs = 'slices'
    gt_obs = 'original_domain'
    adata = adata_raw.copy()
    
    # If celltype exists, add it to the observation data
    if 'celltype' in adata.obs.columns:
        adata.obs[batch_obs] = adata.obs[batch_obs].astype('category')
        adata.obs['ct'] = adata.obs['celltype'].astype('category')
    else:
        # Perform normalization if specified
        if runNormalization:
            sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=hvgs)
            sc.pp.normalize_total(adata, inplace=True)
            sc.pp.log1p(adata)
            sc.pp.pca(adata, n_comps=n_pcs)
            sc.pp.neighbors(adata)
            sc.tl.leiden(adata, resolution=2, key_added='ct', random_state=666)
            adata.obs[batch_obs] = adata.obs[batch_obs].astype('category')
            adata.obs['ct'] = adata.obs['ct'].astype('category')
        else:
            # Skip normalization if not required
            sc.pp.pca(adata, n_comps=n_pcs)
            sc.pp.neighbors(adata)
            sc.tl.leiden(adata, resolution=2, key_added='ct', random_state=666)
            adata.obs[batch_obs] = adata.obs[batch_obs].astype('category')
            adata.obs['ct'] = adata.obs['ct'].astype('category')

    # Run the MENDER model and gather performance stats
    memory_stats, time_stats, adata = run_model(adata, n_cls, tech, scale, radius, seed)
    
    # Save the memory usage and execution time statistics to a CSV file
    stats_file = os.path.join(output_file, f"{sample}_stats.csv")
    with open(stats_file, "w", newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Memory_Usage_MiB", "Execution_Time_s"])
        csvwriter.writerow([memory_stats, time_stats])

    # Save the processed AnnData object to a file
    output_path = os.path.join(output_file, f"{sample}.h5ad")
    adata.write(output_path)


if __name__ == "__main__":
    """
    Main execution block that parses command-line arguments and processes the data.

    This section provides the entry point for running the script with specified arguments, 
    including the input file, output directory, sample name, and model parameters.
    """
    
    # Parse the command-line arguments
    parser = argparse.ArgumentParser(description='Process some single-cell data with MENDER.')
    parser.add_argument('--input_file', type=str, required=True, help='Input file path')
    parser.add_argument('--output_file', type=str, required=True, help='Output file path')
    parser.add_argument('--sample', type=str, required=True, help='Sample name for analysis')
    parser.add_argument('--nclust', type=int, required=True, help='Number of clusters')
    parser.add_argument('--hvgs', type=int, default=4000, help='Number of highly variable genes (default: 4000)')
    parser.add_argument('--runNormalization', type=bool, default=True, help='Run data normalization before integration')
    parser.add_argument('--tech', type=str, required=True, help='Sequencing technology (e.g., "ST", "Visium")')
    parser.add_argument('--scale', type=int, default=6, help='Scale parameter for MENDER (default: 6)')
    parser.add_argument('--radius', type=int, default=15, help='Radius parameter for MENDER (default: 15)')
    parser.add_argument('--n_pcs', type=int, default=50, help='Number of principal components (default: 50)')
    parser.add_argument('--seed', type=int, default=666, help='Random seed')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the process_data function with the parsed arguments
    process_data(args.input_file, args.output_file, args.sample, args.nclust, args.hvgs, args.runNormalization, args.tech, args.scale, args.radius, args.n_pcs, args.seed)
    print("success!")
