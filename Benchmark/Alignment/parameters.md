# Rotate_spatial.py
```python
input_file: Each slice needs to be saved separately to an h5ad file, and the path represents the location where these slices are stored. Refer to the example data for the specific format. This is a required parameter.
output_file: Path to the directory where the output results will be saved. This is a required parameter.
batches: The name of the file corresponding to each slice. Different slices are separated by "," and need to be arranged in adjacent order. For example: "slices1,slices2,slices3". This is a required parameter.
angles: The angles by which each slice needs to be rotated sequentially, with the length matching the number of slices. Different angles are separated by spaces, for example, 40 20 0.
```
# Run_STAligner.py
```python
input_file: Each slice needs to be saved separately to an h5ad file, and the path represents the location where these slices are stored. Refer to the example data for the specific format. This is a required parameter.
input_data: Path to save the integration data. This is a required parameter.
output_file:  Path to the directory where the output results will be saved. This is a required parameter.
batches: The name of the file corresponding to each slice. Different slices are separated by "," and need to be arranged in adjacent order. For example: "slices1,slices2,slices3". This is a required parameter.
landmark_domain: The name of the selected landmark domain. This name is the domian name recognized based on the integration method. This is a required parameter.
landmark_domain_original: The original domain name corresponding to the selected landmark domain. This is a required parameter.
domain: The name of metadata that store identified domians. This is a required parameter.
step: Hyperparameter for the STAligner model. This is a required parameter.
runNormalization: A boolean variable indicating whether data normalization should be performed before integration. This is a required parameter.
hvgs: The number of highly variable genes identified if data normalization is performed. This is an optional parameter, with the default set to 10000.
```
