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

def run_model(adata, n_cls, tech,scale, radius, seed):
    np.random.seed(seed)
    def model_execution():
        batch_obs = 'slices'
        gt_obs = 'original_domain'
        msm = MENDER.MENDER(
            adata,
            batch_obs='slices',
            ct_obs='ct',
            random_seed=seed
        )
        msm.prepare()
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
          
        msm.run_representation_mp(2)
        msm.run_clustering_normal(n_cls)
        
        barcodes_adata = adata.obs['barcode'].values
        barcodes_mender = msm.adata_MENDER.obs['barcode'].values
        adata_df = pd.DataFrame(index=barcodes_adata)
        
        mender_X_df = pd.DataFrame(msm.adata_MENDER.X, index=barcodes_mender)
        mender_X_reordered = mender_X_df.reindex(barcodes_adata).values
        adata.obsm["X_embedding"] = mender_X_reordered
        
        # 提取并重排序msm.adata_MENDER.obs["MENDER"]
        mender_obs_df =  pd.DataFrame(msm.adata_MENDER.obs[['MENDER']].values, index=barcodes_mender,columns=['MENDER'])
        mender_obs_reordered = mender_obs_df.reindex(barcodes_adata)['MENDER'].values
        adata.obs["predicted_domain"] = mender_obs_reordered
        return adata
    start_time = time.time()
    mem_usage = memory_usage((model_execution,), max_usage=True, retval=True)
    end_time = time.time()
    adata = mem_usage[1]  # 获取函数的返回值
    mem_usage_values = mem_usage[0]  # 获取内存使用情况
    return mem_usage_values, end_time - start_time, adata

def process_data(input_file, output_file, sample, n_cls, hvgs, runNormalization, tech, scale, radius, n_pcs, seed):
    adata_raw = sc.read_h5ad(input_file)
    adata_raw.obs['slices'] = adata_raw.obs['slices'].astype('category')
    adata_raw.obs['original_domain'] = adata_raw.obs['original_domain'].astype('category')
    batch_obs = 'slices'
    gt_obs = 'original_domain'
    adata = adata_raw.copy()
    
    if 'celltype' in adata.obs.columns:
      adata.obs[batch_obs] = adata.obs[batch_obs].astype('category')
      adata.obs['ct'] = adata.obs['celltype'].astype('category')
    else:
      if runNormalization:
        sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=hvgs)
        sc.pp.normalize_total(adata, inplace=True)
        sc.pp.log1p(adata)
  
        sc.pp.pca(adata, n_comps=n_pcs)
        sc.pp.neighbors(adata)
        sc.tl.leiden(adata,resolution=2,key_added='ct',random_state=666)
        adata.obs[batch_obs] = adata.obs[batch_obs].astype('category')
        adata.obs['ct'] = adata.obs['ct'].astype('category')
      else:
        sc.pp.pca(adata, n_comps=n_pcs)
        sc.pp.neighbors(adata)
        sc.tl.leiden(adata,resolution=2,key_added='ct',random_state=666)
        adata.obs[batch_obs] = adata.obs[batch_obs].astype('category')
        adata.obs['ct'] = adata.obs['ct'].astype('category')
    
    memory_stats, time_stats, adata = run_model(adata, n_cls, tech,scale, radius, seed)
    
    # Save the statistics to a CSV file
    stats_file = os.path.join(output_file, f"{sample}_stats.csv")
    with open(stats_file, "w", newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Memory_Usage_MiB", "Execution_Time_s"])
        csvwriter.writerow([memory_stats, time_stats])
        
    output_path = os.path.join(output_file, f"{sample}.h5ad")
    adata.write(output_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some single-cell data with MENDER.')
    parser.add_argument('--input_file', type=str, required=True, help='Input file path')
    parser.add_argument('--output_file', type=str, required=True, help='Output file path')
    parser.add_argument('--sample', type=str, required=True, help='Sample name you analysis')
    parser.add_argument('--nclust', type=int, required=True, help='Number of clusters')
    parser.add_argument('--hvgs', type=int, default=4000, help='Number of highly variable genes (default: 4000)')
    parser.add_argument('--runNormalization', type=bool, default=True, help='Run data normalization before integration')
    parser.add_argument('--tech', type=str, required=True, help='Sequencing techonology')
    parser.add_argument('--scale', type=int, default=6, help='Scale parameter for MENDER (default: 6)')
    parser.add_argument('--radius', type=int, default=15, help='Radius parameter for MENDER (default: 15)')
    parser.add_argument('--n_pcs', type=int, default=50, help='Number of principal components (default: 50)')
    parser.add_argument('--seed', type=int, default = 666, help='random seed')

    args = parser.parse_args()

    process_data(args.input_file, args.output_file, args.sample, args.nclust, args.hvgs, args.runNormalization, args.tech, args.scale, args.radius, args.n_pcs, args.seed)
    print("success!")
    
    
