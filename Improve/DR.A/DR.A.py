import numpy as np
from scipy.spatial import KDTree
from scipy.spatial.distance import pdist, squareform
import scanpy as sc
from mnn_utils import create_dictionary_mnn
import math
import matplotlib.pyplot as plt
import os
import anndata as ad
import argparse

def domain_dist_cal(source_data):
    domain_mean_spots = []
    domain_names = []
    for idx, domain in enumerate(np.unique(source_data.obs['original_domain'])):
        domain_mean_spot = source_data[source_data.obs['original_domain']==domain].obsm['X_pca'].mean(axis=0)
        domain_mean_spots.append(np.array(domain_mean_spot))
        domain_names.append(domain)
    domain_mean_spots = np.array(domain_mean_spots)
    distances = squareform(pdist(domain_mean_spots, 'euclidean'))
    return domain_names, distances

def compute_weights(distances, sigma=1.0):
    return np.exp(-distances**2 / (2 * sigma**2))

def find_correspondences(source_data, target_data, source_domains, target_domains, src_idx, tgt_idx, alpha=0.5):
    source_points = source_data.obsm['X_pca']
    target_points = target_data.obsm['X_pca']
    # source_tree = KDTree(source_points)
    source_domain_names, source_domain_dist = domain_dist_cal(source_data)
    correspondences = []
    for i, src_idx_tmp in enumerate(src_idx):
        dist = np.linalg.norm(source_points[src_idx_tmp] - target_points[tgt_idx[i]])
        domain_dist = source_domain_dist[source_domain_names.index(source_domains[src_idx_tmp]), source_domain_names.index(target_domains[tgt_idx[i]])]
        combined_dist = alpha * dist + (1 - alpha) * domain_dist
        correspondences.append((src_idx_tmp, tgt_idx[i], combined_dist))
    return correspondences

def estimate_transformation(source_points, target_points, correspondences, weights):
    src_pts = np.array([source_points[i] for i, _ in correspondences])
    tgt_pts = np.array([target_points[j] for _, j in correspondences])
    
    centroid_src = np.average(src_pts, axis=0, weights=weights)
    centroid_tgt = np.average(tgt_pts, axis=0, weights=weights)
    
    src_pts_centered = src_pts - centroid_src
    tgt_pts_centered = tgt_pts - centroid_tgt
    H = np.dot(src_pts_centered.T, tgt_pts_centered)
    U, S, Vt = np.linalg.svd(H)
    R = np.dot(Vt.T, U.T)
    if np.linalg.det(R) < 0:
        Vt[1, :] *= -1
        R = np.dot(Vt.T, U.T)
    
    t = centroid_tgt - np.dot(R, centroid_src)
    return R, t

def apply_transformation(points, R, t):
    return np.dot(points, R.T) + t

def domain_icp(source, target, src_idx, tgt_idx, max_iterations=20, tolerance=1e-15, sigma = 1.0):
    source_points = np.asarray(source.obsm['X_pca'])
    target_points = np.asarray(target.obsm['X_pca'])
    source_domains = np.asarray(source.obs['original_domain'])
    target_domains = np.asarray(target.obs['original_domain'])
    
    prev_error = float('inf')
    transformation = np.eye(3)
    
    # correspondences = find_correspondences(source, target, source_domains, target_domains, src_idx, tgt_idx, alpha)
    target_domain_names, target_domain_dist = domain_dist_cal(target)
    weights_target = compute_weights(target_domain_dist, sigma)
    
    weights = []
    for i_src_idx, i_src in enumerate(source_domains[src_idx]):
      try:
        weights.append(weights_target[target_domain_names.index(i_src), target_domain_names.index(target_domains[tgt_idx][i_src_idx])])
      except:
        weights.append(0)
        
    
    correspondences = list(zip(src_idx, tgt_idx))
    
    coor_src = source.obsm['spatial']
    for i in range(max_iterations):
        print(i)
        R, t = estimate_transformation(coor_src, target.obsm['spatial'], correspondences, weights)
        
        coor_src = apply_transformation(coor_src, R, t)
        
        current_error = np.mean([math.sqrt(((coor_src[src_idx][kk, 0] - target.obsm['spatial'][tgt_idx][kk, 0]) ** 2) + ((coor_src[src_idx][kk, 1] - target.obsm['spatial'][tgt_idx][kk, 1]) ** 2))
                             for kk in range(len(coor_src[src_idx]))])
        
        # current_error = sum(dist for _, _, dist in correspondences) / len(correspondences)
        
        if np.abs(prev_error - current_error) < tolerance:
            break
        
        prev_error = current_error
        
        T = np.eye(3)
        T[:2, :2] = R
        T[:2, 2] = t
        transformation = T @ transformation
    
    return transformation

def find_mnn_idx(new_data, source_slice='151675', target_slice='151673'):
    mnn_dict = create_dictionary_mnn(new_data, use_rep='X_pca', batch_name='slices', k=1, iter_comb=None, verbose=0)
    adata_1 = new_data[new_data.obs['slices'] == source_slice]
    adata_2 = new_data[new_data.obs['slices'] == target_slice]
    
    anchor_list = []
    positive_list = []
    for batch_pair_name in mnn_dict.keys(): 
        for anchor in mnn_dict[batch_pair_name].keys():
            positive_spot = mnn_dict[batch_pair_name][anchor][0]
            ### anchor should only in the ref slice, pos only in the target slice
            if anchor in adata_1.obs_names and positive_spot in adata_2.obs_names:                 
                anchor_list.append(anchor)
                positive_list.append(positive_spot)
                
    batch_as_dict = dict(zip(list(new_data.obs_names), range(0, new_data.shape[0])))
    anchor_ind = list(map(lambda _: batch_as_dict[_], anchor_list))
    positive_ind = list(map(lambda _: batch_as_dict[_], positive_list))
    anchor_arr = new_data.obsm['X_pca'][anchor_ind, ]
    positive_arr = new_data.obsm['X_pca'][positive_ind, ]
    dist_list = [np.sqrt(np.sum(np.square(anchor_arr[ii, :] - positive_arr[ii, :]))) for ii in range(anchor_arr.shape[0])]
    
    key_points_src =  np.array(anchor_list)[dist_list < np.percentile(dist_list, 50)] ## remove remote outliers
    key_points_dst =  np.array(positive_list)[dist_list < np.percentile(dist_list, 50)]
    
    # coor_src = adata_1.obsm["spatial"] ## to_be_aligned
    # coor_dst = adata_2.obsm["spatial"] ## reference_points
    
    ## index number
    MNN_ind_src = [list(adata_1.obs_names).index(key_points_src[ii]) for ii in range(len(key_points_src))]
    MNN_ind_dst = [list(adata_2.obs_names).index(key_points_dst[ii]) for ii in range(len(key_points_dst))]
    return MNN_ind_src, MNN_ind_dst, adata_1, adata_2

def get_aligned_coor(transformation_matrix, coor):
    R = transformation_matrix[:2, :2]
    t = transformation_matrix[:2, 2]
    aligned_coor = apply_transformation(coor, R, t)
    return aligned_coor

def correct_coor(new_data, source_slice, target_slice):
    MNN_ind_src, MNN_ind_dst, adata_source, adata_target = find_mnn_idx(new_data, source_slice=source_slice, target_slice=target_slice)
    transformation = domain_icp(adata_source, adata_target, MNN_ind_src, MNN_ind_dst)
    aligned_coor = get_aligned_coor(transformation, adata_source.obsm['spatial'])
    return aligned_coor

def process_data(data_path, integrated_path, slices_file, slices, output_path, model='original'):
    slices_file = slices_file.split(",")
    slices = slices.split(",")
    
    # Construct paths based on dataset and model type
    data_paths = [os.path.join(data_path, f"{s}.h5ad") for s in slices_file]

    # Read data based on model type
    if model == 'original':
        data = [sc.read_h5ad(path) for path in data_paths]
    else:
        data = [sc.read_h5ad(path) for path in data_paths]
        
        predict_data = sc.read_h5ad(os.path.join(integrated_path, f"{model}.h5ad"))
        predict_data.obs["predicted_domain"] = predict_data.obs["predicted_domain"].astype(str)
        predict_data_obs = predict_data.obs[['barcode', "predicted_domain"]].copy()
        
        for d in data:
            d.obs["slices"] = d.obs["slices"].astype(str)
            adata_obs = d.obs
            merged_obs = adata_obs.merge(predict_data_obs, on='barcode', how='left')
            merged_obs.set_index(adata_obs.index, inplace=True)
            d.obs["original_domain"] = merged_obs["predicted_domain"]
    
    new_data = ad.concat(data, pairwise=True) 
    new_data.obs_names_make_unique()
    new_data.obs = new_data.obs.astype(str)
    # Preprocessing steps
    sc.pp.scale(new_data)
    sc.pp.pca(new_data, n_comps=30, svd_solver='arpack')
    
    # Spatial alignment
    corrected_coor = []
    for ii in range(len(slices)-1):
        aligned_coor = correct_coor(new_data, source_slice=slices[ii], target_slice=slices[-1])
        corrected_coor.append(aligned_coor)
    
    corrected_coor.append(data[-1].obsm['spatial'])
    all_coordinates = np.vstack([np.array(arr) for arr in corrected_coor])
    new_data.obsm["spatial"] = all_coordinates
    new_data.obs = new_data.obs.astype(str)
    
    # Determine output path based on model type
    if model == 'original':
        output_dir = os.path.join(output_path,'DR.A')
    else:
        output_dir = os.path.join(output_path,f'DR.A_{model}')
    
    # Create directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Write processed data
    new_data.write_h5ad(os.path.join(output_dir, 'Spatial_correct_data.h5ad'))
    
    # Plotting
    fig, axes = plt.subplots(1, len(slices), figsize=(6*len(slices), 6))
    
    for i, ax in enumerate(axes):
        cols = dict(zip(np.unique(data[i].obs['original_domain']), range(len(np.unique(data[i].obs['original_domain'])))))
        ax.scatter(corrected_coor[i][:,0], corrected_coor[i][:,1], c=[cols[j] for j in data[i].obs['original_domain']])
        ax.set_title(f'Slice {slices[i]}')
        
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'tmp.png'))
  

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--data_path', type=str, help='input file path')
    parser.add_argument('--integrated_path', type=str, default = "none", help='input integration results path')
    parser.add_argument('--slices_file', type=str, help='the file name of each slice data')
    parser.add_argument('--slices', type=str, help='the slice label')
    parser.add_argument('--output_path', type=str, help='output file path')
    parser.add_argument('--model', type=str, default="original", help='the integration model')
    args = parser.parse_args()
    
    process_data(args.data_path, args.integrated_path, args.slices_file, args.slices, args.output_path, args.model)
    
    
    
