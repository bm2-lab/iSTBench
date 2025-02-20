# SlicesRepresentation
## Parameters:
```python
file (abbreviation: f): The file path that stores the result of data consolidation when the number of domains is set. Refer to the example data for the specific format. This is a required parameter.
model (abbreviation: m): The model for consolidating the data is also the name of the folder next to "file". This is a required parameter.
cluster_number (abbreviation: n): The number of classes of samples. This is a required parameter.
cluster_method (abbreviation: M): The method used for sample clustering can be set to "hclust" and "leiden". This is an optional parameter, with the default set to "hclust".
predicted_domain (abbreviation: p): The name of metadata that store identified domians. This is an optional parameter, with the default set to "predicted_domain".
dist_method (abbreviation: d): The method of calculating sample distance when doing sample clustering. Can be set to "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". This is an optional parameter, with the default set to "euclidean".
slices_class (abbreviation: s): The name of metadata that store sample class. This is an optional parameter, with the default set to "slices_class".
```
## Terminal running
```python
nohup /opt/R/4.3.2/lib/R/bin/Rscript Benchmark/SliceRepresentation/SliceRepresentation.R \
-f /NFS2_home/NFS2_home_3/dongkj/home_dkj/FD_yzy/Result/GitHub_test/iSTBench/Data/TNBC/SlicesEmbedding/MENDER \
-m MENDER -n 3 \
> Data/TNBC/SlicesEmbedding/MENDER/Metric/SlicesClustering.output &
```
