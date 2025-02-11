import scanpy as sc
import scib
import os
import anndata as ad
import pandas as pd
import argparse
import numpy as np
import random

def calculate_metrics(input_file, input_path, output_path, seed = 123):
    # set random seed
    np.random.seed(seed)
    random.seed(seed)
    
    fl = [f for f in os.listdir(input_file) if f.endswith('.h5ad') and 'withRaw' not in f]

    results = []
    for f in fl:
        model = os.path.splitext(f)[0]
        inter_data = ad.read_h5ad(os.path.join(input_file, f))
        inter_data.obs["slices"] = inter_data.obs["slices"].astype(str).astype('category')
        inter_data.obs["original_domain"] = inter_data.obs["original_domain"].astype(str).astype('category')

        original_data = ad.read_h5ad(input_path)
        original_data.obs["slices"] = original_data.obs["slices"].astype(str).astype('category')
        original_data.obs["original_domain"] = original_data.obs["original_domain"].astype(str).astype('category')

        if model == "GraphST" or model == "GraphSTwithPASTE":
            embed = "emb"
        else:
            embed = "X_embedding"

        # 计算指标
        #embedding-domain/batch
        dASW = scib.me.silhouette(inter_data, label_key="original_domain", embed = embed, scale = True)
        bASW = scib.me.silhouette_batch(inter_data, batch_key="slices", label_key="original_domain", embed=embed, scale = True)
        
        #knn-domian/batch
        sc.pp.neighbors(inter_data, use_rep=embed)
        ilF1  = scib.me.isolated_labels_f1(inter_data, batch_key="slices", label_key="original_domain", embed = None)
        kBET = scib.me.kBET(inter_data, batch_key="slices", label_key="original_domain", type_="knn", scaled = True)
        
        # PCA = scib.me.pcr_comparison(original_data, inter_data, covariate="slices", embed=embed, n_comps = 20, scale = True)
        # GC = scib.me.graph_connectivity(inter_data, label_key="original_domain")

        # 将结果添加到列表中
        results.append([model, dASW, bASW, ilF1, kBET])
        #results.append([model, dASW, bASW, ilF1, kBET, PCA, GC])

    # 将结果转换为DataFrame并保存为CSV文件
    results_df = pd.DataFrame(results, columns=["Model", "dASW", "bASW", "ilF1", "kBET"])
    #results_df = pd.DataFrame(results, columns=["Model", "dASW", "bASW", "ilF1", "kBET", "PCA", "GC"])
    results_df.to_csv(output_path, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate integration metrics for spatial transcriptomics data.")
    
    parser.add_argument('--input_file', type=str, required=True, 
                        help="Directory containing the integration result files in .h5ad format.")
    parser.add_argument('--input_path', type=str, required=True, 
                        help="Path to the original combined slices data in .h5ad format.")
    parser.add_argument('--output_path', type=str, required=True, 
                        help="Path where the output CSV file with metrics will be saved.")
    
    args = parser.parse_args()
    
    calculate_metrics(args.input_file, args.input_path, args.output_path)
