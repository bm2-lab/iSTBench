import os
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import math

# Main data processing function.
# This function reads in spatial data, rotates the spatial coordinates for each specified slice (batch),
# and writes the rotated data to the output directory. It processes multiple batches with corresponding rotation angles.
def process_data(input_file, output_file, batches, angles):
    # input_file: Each slice needs to be saved separately to an h5ad file, and the path represents the location where these slices are stored. Refer to the example data for the specific format. This is a required parameter.
    # output_file: Path to the directory where the output results will be saved. This is a required parameter.
    # batches: The name of the file corresponding to each slice. Different slices are separated by "," and need to be arranged in adjacent order. For example: "slices1,slices2,slices3". This is a required parameter.
    # angles: The angles by which each slice needs to be rotated sequentially, with the length matching the number of slices. Different angles are separated by spaces, for example, 40 20 0.

  
    # Ensure the output directory exists
    if not os.path.exists(output_file):
        os.makedirs(output_file)
    
    # Convert the comma-separated batches string into a list
    batches = batches.split(",")
    
    # Process each batch with the corresponding rotation angle
    for batche, angle in zip(batches, angles):
        print(f"Processing {batche} with rotation angle {angle} degrees")
        
        # Read in the AnnData object for the current batch
        adata = sc.read_h5ad(f"{input_file}/{batche}.h5ad")
        
        # Extract the original spatial coordinates
        original_spatial = adata.obsm["spatial"]
        
        # Apply rotation to the spatial coordinates
        rotate_spatial = rotate(original_spatial, angle)
        
        # Update the AnnData object with the rotated coordinates
        adata.obsm["spatial"] = rotate_spatial
        
        # Write the updated AnnData object to the output directory
        adata.write_h5ad(f"{output_file}/{batche}.h5ad")
    
    print("Finished!")


# Function to rotate the given spatial coordinates by a specified angle.
# This function performs a 2D rotation on the coordinates around their mean center.
# The rotation is applied by first translating the coordinates to the origin, performing the rotation, 
# and then translating them back to their original positions.
def rotate(coords, angle):
    # Convert angle from degrees to radians
    theta = np.radians(angle)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array(((c, -s), (s, c)))  # Rotation matrix

    # Compute the center of the coordinates (mean position)
    center = np.mean(coords, axis=0)
    
    # Translate coordinates to origin, apply rotation, and translate back to original position
    rotated_coords = np.dot(coords - center, R) + center
    
    return rotated_coords


# Command-line interface (CLI) setup.
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process spatial data with rotation.")
    
    # Argument for input file directory
    parser.add_argument("--input_file", type=str, required=True, help="Path to the input file directory.")
    
    # Argument for output file directory
    parser.add_argument("--output_file", type=str, required=True, help="Path to the output file directory.")
    
    # Argument for specifying the sequence of slices (batches) to process
    parser.add_argument("--batches", type=str, required=True, help="The sequence of slices")
    
    # Argument for specifying the list of rotation angles corresponding to each batch
    parser.add_argument("--angles", type=int, nargs='+', required=True, help="List of rotation angles.")
    
    args = parser.parse_args()
    
    # Execute the main processing function with the parsed arguments
    process_data(args.input_file, args.output_file, args.batches, args.angles)
