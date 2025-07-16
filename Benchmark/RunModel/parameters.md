# Run_Banksy
```python
input_file (abbreviation: i): Path to the multi-slice data that needs to be integrated. All slices should be contained in a single RDS or RData (the variable name after loading should be data) file. Refer to the example data for the specific format. This is a required parameter.

output_file (abbreviation: o): Path to the directory where the output results will be saved. This is a required parameter.

sample (abbreviation: s): Name for the output result files. It is recommended to use the model name, for example, setting it as "Banksy" will result in output files named "Banksy.h5ad". This is a required parameter.

nclust (abbreviation: c): Number of identified domains. This is a required parameter.

dataType (abbreviation: d): The type of input data, which should be set to "RData" or "RDS". This is an optional parameter, with the default set to "RData".

hvgs (abbreviation: h): The number of highly variable genes identified if data normalization is performed. This is an optional parameter, with the default set to 2000.

runNormalization (abbreviation: n): A boolean variable indicating whether data normalization should be performed before integration. This is an optional parameter, with the default set to TRUE.

k_geom (abbreviation: k): Hyperparameter for the Banksy model. This is an optional parameter, with the default set to c(15, 30).

lambda (abbreviation: l): Hyperparameter for the Banksy model. This is an optional parameter, with the default set to 0.8.

SEED (abbreviation: S): Random seed for model execution. This is an optional parameter, with the default set to 1234.

```
# Run_CellCharter
```python
input_file: Path to the multi-slice data that needs to be integrated. All slices should be contained in a single h5ad file. Refer to the example data for the specific format. This is a required parameter.

output_file: Path to the directory where the output results will be saved. This is a required parameter.

sample: Name for the output result files. It is recommended to use the model name, for example, setting it as "CellCharter" will result in output files named "CellCharter.h5ad". This is a required parameter.

nclust: Number of identified domains. This is a required parameter.

hvgs: The number of highly variable genes identified if data normalization is performed. This is an optional parameter, with the default set to 5000.

runNormalization: A boolean variable indicating whether data normalization should be performed before integration. This is an optional parameter, with the default set to True.

isCODEX: A boolean variable indicates whether the data were obtained by CODEX sequencing. This is an optional parameter, with the default set to False.

n_latent: Hyperparameter for the CellCharter model. This is a required parameter. The original recommendation is to conduct a grid search in the hyperparameter space [5, 10, 15] to select the best hyperparameter.

nhood_layers: Hyperparameter for the CellCharter model. This is a required parameter. The original recommendation is to conduct a grid search in the hyperparameter space [2, 3, 4] to select the best hyperparameter.

seed: Random seed for model execution. This is an optional parameter, with the default set to 1234.
```
# Run_CN
```python
input_file: Path to the multi-slice data that needs to be integrated. All slices should be contained in a single h5ad file. Refer to the example data for the specific format. This is a required parameter.

output_file: Path to the directory where the output results will be saved. This is a required parameter.

sample: Name for the output result files. It is recommended to use the model name, for example, setting it as "CN" will result in output files named "CN.h5ad". This is a required parameter.

nclust: Number of identified domains. This is a required parameter.

runNormalization: A boolean variable indicating whether data normalization should be performed before integration. This is an optional parameter, with the default set to True.

hvgs: The number of highly variable genes identified if data normalization is performed. This is an optional parameter, with the default set to 4000.

n_pcs: The number of principal components set when performing PCA dimensionality reduction on the data. This is an optional parameter, with the default set to 50.

seed: Random seed for model execution. This is an optional parameter, with the default set to 0.
```
# Run_GraphST/Run_GraphST-PASTE
```python
input_file: Path to the multi-slice data that needs to be integrated. All slices should be contained in a single h5ad file. Refer to the example data for the specific format. This is a required parameter.

output_file: Path to the directory where the output results will be saved. This is a required parameter.

sample: Name for the output result files. It is recommended to use the model name, for example, setting it as "GraphST" will result in output files named "GraphST.h5ad". This is a required parameter.

nclust: Number of identified domains. This is a required parameter.

hvgs: The number of highly variable genes identified if data normalization is performed. This is an optional parameter, with the default set to 3000.

runNormalization: A boolean variable indicating whether data normalization should be performed before integration. This is an optional parameter, with the default set to True.

tool: Methods for spatial clustering, which can be set to "mclust", "leiden" and "louvain". This is an optional parameter, with the default set to "mclust".

device: The device on which the model runs can be set to "cpu" and "gpu". This is an optional parameter, with the default set to "cpu".

seed: Random seed for model execution. This is an optional parameter, with the default set to 0.
```
# Run_MENDER
```python
input_file: Path to the multi-slice data that needs to be integrated. All slices should be contained in a single h5ad file. Refer to the example data for the specific format. This is a required parameter.

output_file: Path to the directory where the output results will be saved. This is a required parameter.

sample: Name for the output result files. It is recommended to use the model name, for example, setting it as "MENDER" will result in output files named "MENDER.h5ad". This is a required parameter.

nclust: Number of identified domains. This is a required parameter.

hvgs: The number of highly variable genes identified if data normalization is performed. This is an optional parameter, with the default set to 4000.

runNormalization: A boolean variable indicating whether data normalization should be performed before integration. This is an optional parameter, with the default set to True.

tech: Sequencing technology that generated the data. Can be set to "ST", "Visium" or others. This is a required parameter.

scale: Hyperparameter for the MENDER model. This is an optional parameter, with the default set to 6.

radius: Hyperparameter for the MENDER model. This is an optional parameter, with the default set to 15.

n_pcs: The number of principal components set when performing PCA dimensionality reduction on the data. This is an optional parameter, with the default set to 50.

seed: Random seed for model execution. This is an optional parameter, with the default set to 666.
```

# Run_NicheCompass
```python
input_file: Each slice needs to be saved separately to an h5ad file, and the path represents the location where these slices are stored. Refer to the example data for the specific format. This is a required parameter.

output_file: Path to the directory where the output results will be saved. This is a required parameter.

sample: Name for the output result files. It is recommended to use the model name, for example, setting it as "NicheCompass" will result in output files named "NicheCompass.h5ad". This is a required parameter.

nclust: Number of identified domains. This is a required parameter.

species: The detected species type can be set to "human" and "mouse". This is an optional parameter, with the default set to "human".

convertID: Whether gene id conversion is required. Genes need to be named with "symbol". This is an optional parameter, with the default set to False.

from_id/to_id: If a gene id conversion is required, enter the original id type and the target id type. They are both optional parameters. The default of to_id is "symbol". The gene id conversion is based on the "mygene" package. The type of gene id allowed to be input can be referred to "mygene".

seed: Random seed for model execution. This is an optional parameter, with the default set to 0.
```
# Run_PRECAST
```python
input_file (abbreviation: i): Each slice needs to be saved separately to an h5ad file, and the path represents the location where these slices are stored. Refer to the example data for the specific format. This is a required parameter.

output_file (abbreviation: o): Path to the directory where the output results will be saved. This is a required parameter.

sample (abbreviation: s): Name for the output result files. It is recommended to use the model name, for example, setting it as "PRECAST" will result in output files named "PRECAST.h5ad". This is a required parameter.

nclust (abbreviation: c): Number of identified domains. This is a required parameter.

hvgs (abbreviation: v): The number of highly variable genes identified if data normalization is performed. This is an optional parameter, with the default set to 2000.

runNormalization (abbreviation: n): A boolean variable indicating whether data normalization should be performed before integration. This is an optional parameter, with the default set to TRUE.

tech (abbreviation: t): Sequencing technology that generated the data. Can be set to "ST", "Visium" or "Other_SRT". This is a required parameter.

species (abbreviation: p): The detected species type can be set to "Human", "Mouse" or "Unknown". This is an optional parameter, with the default set to "human".

seed (abbreviation: S): Random seed for model execution. This is an optional parameter, with the default set to 1234.
```
# Run_SpaDo
```python
input_file (abbreviation: i): Path to the multi-slice data that needs to be integrated. All slices should be contained in a single h5ad file. Refer to the example data for the specific format. This is a required parameter.

output_file (abbreviation: o): Path to the directory where the output results will be saved. This is a required parameter.

sample (abbreviation: s): Name for the output result files. It is recommended to use the model name, for example, setting it as "SpaDo" will result in output files named "SpaDo.h5ad". This is a required parameter.

nclust (abbreviation: c): Number of identified domains. This is a required parameter.

hvgs (abbreviation: v): The number of highly variable genes identified if data normalization is performed. This is an optional parameter, with the default set to 2000.

runNormalization (abbreviation: n): A boolean variable indicating whether data normalization should be performed before integration. This is an optional parameter, with the default set to TRUE.

n_pcs: The number of principal components set when performing PCA dimensionality reduction on the data. This is an optional parameter, with the default set to 50.
```

# Run_SPIRAL
```python
input_file: Each slice needs to be saved separately to an h5ad file, and the path represents the location where these slices are stored. Refer to the example data for the specific format. This is a required parameter.

output_file: Path to the directory where the output results will be saved. This is a required parameter.

model_path: Path to store the intermediate results. This is a required parameter.

sample: Name for the output result files. It is recommended to use the model name, for example, setting it as "SPIRAL" will result in output files named "SPIRAL.h5ad". This is a required parameter.

hvgs: The number of highly variable genes identified if data normalization is performed. This is an optional parameter, with the default set to 1000.

runNormalization: A boolean variable indicating whether data normalization should be performed before integration. This is an optional parameter, with the default set to TRUE.

KNN: In constructing the k-nearest neighbor graph, the choice of k determines the number of nearest neighbors considered for each cell or spot, with the default set to 6.

data: The name of the data. This is a required parameter.

seed: Random seed for model execution. This is an optional parameter, with the default set to 1234.

nclust: Number of identified domains. This is a required parameter.
```
# Run_STAIG
```python
input_file: Each slice needs to be saved separately to an h5ad file, and the path represents the location where these slices are stored. Refer to the example data for the specific format. This is a required parameter.

output_file: Path to the directory where the output results will be saved. This is a required parameter.

model_path: Path to store the intermediate results. This is a required parameter.

sample: Name for the output result files. It is recommended to use the model name, for example, setting it as "STAIG" will result in output files named "STAIG.h5ad". This is a required parameter.

hvgs: The number of highly variable genes identified if data normalization is performed. This is an optional parameter, with the default set to 3000.

runNormalization: A boolean variable indicating whether data normalization should be performed before integration. This is an optional parameter, with the default set to TRUE.

nclust: Number of identified domains. This is a required parameter.

nPCA: Number of PCA, with the default set to 64.

num_layers: Number of GCN layers.

k_neighbor: Number of neighboring nodes for each spot.

tau: Temperature parameter.

num_epochs: Number of epochs.

#The original paper provides reference values for parameters such as num_layers, k_neighbor, tau, and num_epochs across different data types.

seed: Random seed for model execution. This is an optional parameter, with the default set to 1234.

device: The device on which the model runs can be set to "cpu" and "gpu". This is an optional parameter, with the default set to "gpu".
```
# Run_STAligner
```python
input_file: Each slice needs to be saved separately to an h5ad file, and the path represents the location where these slices are stored. Refer to the example data for the specific format. This is a required parameter.

output_file: Path to the directory where the output results will be saved. This is a required parameter.

sample: Name for the output result files. It is recommended to use the model name, for example, setting it as "STAligner" will result in output files named "STAligner.h5ad". This is a required parameter.

batches: Each slice name should be specified according to its sequential order. This is a required parameter.

hvgs: The number of highly variable genes identified if data normalization is performed. This is an optional parameter, with the default set to 5000.

runNormalization: A boolean variable indicating whether data normalization should be performed before integration. This is an optional parameter, with the default set to TRUE.

nclust: Number of identified domains. This is a required parameter.

r: Radius for identify neighbors, with the default set to 50.

seed: Random seed for model execution. This is an optional parameter, with the default set to 1234.

device: The device on which the model runs can be set to "cpu" and "gpu". This is an optional parameter, with the default set to "gpu".
```
