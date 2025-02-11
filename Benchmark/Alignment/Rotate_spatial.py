import os
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import math

def rotate(coords, angle):
    # 将角度从度数转换为弧度
    theta = np.radians(angle)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array(((c, -s), (s, c)))  # 旋转矩阵

    # 计算中心点
    center = np.mean(coords, axis=0)
    
    # 平移到原点，旋转，平移回原位置
    rotated_coords = np.dot(coords - center, R) + center
    
    return rotated_coords

def process_data(input_file, output_file, batches, angles):
    if not os.path.exists(output_file):
        os.makedirs(output_file)
    
    batches = batches.split(",")
    # batches = os.listdir(input_file)
    # batches = [batch.split('.')[0] for batch in batches]
    # batches = sorted(batches, key=lambda x: int(x.replace('slices', '')))
    
    for batche, angle in zip(batches, angles):
        print(f"Processing {batche} with rotation angle {angle} degrees")
        adata = sc.read_h5ad(f"{input_file}/{batche}.h5ad")
        original_spatial = adata.obsm["spatial"]
        rotate_spatial = rotate(original_spatial, angle)
        
        adata.obsm["spatial"] = rotate_spatial
        
        adata.write_h5ad(f"{output_file}/{batche}.h5ad")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process spatial data with rotation.")
    parser.add_argument("--input_file", type=str, required=True, help="Path to the input file directory.")
    parser.add_argument("--output_file", type=str, required=True, help="Path to the output file directory.")
    parser.add_argument("--batches", type=str, required=True, help="The sequence of slices")
    parser.add_argument("--angles", type=int, nargs='+', required=True, help="List of rotation angles.")
    
    args = parser.parse_args()
    
    process_data(args.input_file, args.output_file, args.batches, args.angles)
