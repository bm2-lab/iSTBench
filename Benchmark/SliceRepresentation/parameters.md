# SlicesRepresentation
## Parameters:
```python
original_file: Path to the multi-slice data that needs to be integrated. All slices should be contained in a single h5ad file. Refer to the example data for the specific format. This is a required parameter.
file: The file path that stores the result of data consolidation when the number of domains is set. Refer to the example data for the specific format. This is a required parameter.
model: The model for consolidating the data is also the name of the folder next to "file". This is a required parameter.
cluster_number: The number of classes of samples. This is a required parameter.
cluster_method: The method used for sample clustering can be set to "hclust" and "kmeans". This is an optional parameter, with the default set to "hclust".
predicted_domain: The name of metadata that store identified domians. This is an optional parameter, with the default set to "predicted_domain".
dist_method: The method of calculating sample distance when doing sample clustering. Can be set to "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". This is an optional parameter, with the default set to "euclidean".
slices_class: The name of metadata that store sample class. This is an optional parameter, with the default set to "slices_class".
```
## Terminal running
```python
nohup Rscript Benchmark/SliceRepresentation/Slices_classification.R \
Data/TNBC/sample_all_data/Slices_combind_data.h5ad \
Data/TNBC/SlicesEmbedding/Banksy \
Banksy 3 hclust \
> Data/TNBC/SlicesEmbedding/Banksy/Metric/SlicesClustering.output &
```
