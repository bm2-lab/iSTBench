# Run methods on BaristaSeq dataset
## Banksy
```python
nohup Rscript Benchmark/RunModel/Run_Banksy.R \
-i Data/BaristaSeq/sample_all_data/Slices_combind_data.RDS \
-o Data/BaristaSeq/IntergrationRe \
-s Banksy -c 7 -d RDS -k 18 -l 0.2 \
> Data/BaristaSeq/IntergrationRe/Banksy.output &
```  
## CellCharter
```python
nohup python Benchmark/RunModel/Run_CellCharter.py \
--input_file Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/BaristaSeq/IntergrationRe \
--sample CellCharter --nclust 7 --hvgs 2000 --runNormalization True \
--n_latent 15 --nhood_layers 2 \
> Data/BaristaSeq/IntergrationRe/CellCharter.output &
```
## CN
```python
nohup python Benchmark/RunModel/Run_CN.py \
--input_file Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/BaristaSeq/IntergrationRe \
--sample CN --nclust 7 \
> Data/BaristaSeq/IntergrationRe/CN.output &
```
## GraphST
```python
nohup python Benchmark/RunModel/GraphST/Run_GraphST.py \
--input_file Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/BaristaSeq/IntergrationRe \
--sample GraphST --nclust 7 --device cuda \
> Data/BaristaSeq/IntergrationRe/GraphST.output  &
```
## GraphST-PASTE
```python
nohup python Benchmark/RunModel/GraphST/Run_GraphST-PASTE.py \
--input_file Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/BaristaSeq/IntergrationRe \
--sample GraphSTwithPASTE --nclust 7 --device cuda \
> Data/BaristaSeq/IntergrationRe/GraphSTwithPASTE.output  &
```
## MENDER
```python
nohup python Benchmark/RunModel/Run_MENDER.py \
--input_file Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/BaristaSeq/IntergrationRe \
--sample MENDER --nclust 7 --tech BaristaSeq \
> Data/BaristaSeq/IntergrationRe/MENDER.output &
```
## NicheCompass 
```python
nohup python Benchmark/RunModel/Run_NicheCompass.py \
--input_file Data/BaristaSeq/sample_data \
--output_file Data/BaristaSeq/IntergrationRe \
--sample NicheCompass --nclust 7 \
> Data/BaristaSeq/IntergrationRe/NicheCompass.output &
```
## Spado
```python
nohup /opt/R/4.3.2/lib/R/bin/Rscript Benchmark/RunModel/Run_Spado.R \
-i Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
-o Data/BaristaSeq/IntergrationRe \
-s Spado -v 2000 -c 7 -n TRUE \
> Data/BaristaSeq/IntergrationRe/Spado.output &    
```
# Run methods on other dataset
## Banksy
### DLPFC S1
```python
nohup Rscript Benchmark/RunModel/Run_Banksy.R \
-i Data/DLPFC_sample1/sample_all_data/Slices_combind_data.RDS \
-o Data/DLPFC_sample1/IntergrationRe \
-s Banksy -c 7 -d RDS -k 18 -l 0.2 \
> Data/DLPFC_sample1/IntergrationRe/Banksy.output &
```
### DLPFC S2
```python
nohup Rscript Benchmark/RunModel/Run_Banksy.R \
-i Data/DLPFC_sample2/sample_all_data/Slices_combind_data.RDS \
-o Data/DLPFC_sample2/IntergrationRe \
-s Banksy -c 5 -d RDS -k 18 -l 0.2 \
> Data/DLPFC_sample2/IntergrationRe/Banksy.output &
```
### DLPFC S3
```python
nohup Rscript Benchmark/RunModel/Run_Banksy.R \
-i Data/DLPFC_sample3/sample_all_data/Slices_combind_data.RDS \
-o Data/DLPFC_sample3/IntergrationRe \
-s Banksy -c 7 -d RDS -k 18 -l 0.2 \
> Data/DLPFC_sample3/IntergrationRe/Banksy.output &
```
### MERFISH Preoptic
```python
nohup Rscript Benchmark/RunModel/Run_Banksy.R \
-i Data/MERFISH/sample_all_data/Slices_combind_data.RDS \
-o Data/MERFISH/IntergrationRe \
-s Banksy -c 8 -d RDS -n FALSE \
> Data/MERFISH/IntergrationRe/Banksy.output &
```
### MERFISH Brain S2-S12
```python
nohup Rscript Benchmark/RunModel/Run_Banksy.R \
-i Data/MERFISH_Brain_S2/sample_all_data/Slices_combind_data.RDS \
-o Data/MERFISH_Brain_S2/IntergrationRe \
-s Banksy -c 8 -d RDS -h 374 \
> Data/MERFISH_Brain_S2/IntergrationRe/Banksy.output &
# You just need to change the path, like MERFISH_Brain_S3, MERFISH_Brain_S4...
```
### STARMap
```python
nohup Rscript Benchmark/RunModel/Run_Banksy.R \
-i Data/STARMap/sample_all_data/Slices_combind_data.RDS \
-o Data/STARMap/IntergrationRe \
-s Banksy -c 4 -d RDS -h 166 \
> Data/STARMap/IntergrationRe/Banksy.output &
```
### Large-scale Dataset
```python
nohup Rscript Benchmark/RunModel/Run_Banksy.R \
-i Data/Mouse/sample_all_data/Slices_combind_data.RDS \
-o Data/Mouse/IntergrationRe \
-s Banksy -c 7 -d RDS -n FALSE \
> Data/Mouse/IntergrationRe/Banksy.output &
```
### TNBC
```python
nohup Rscript Benchmark/RunModel/Run_Banksy.R \
-i Data/TNBC/sample_all_data/Slices_combind_data.RDS \
-o Data/TNBC/IntergrationRe \
-s Banksy -c 4 -d RDS -n FALSE \
> Data/TNBC/SlicesEmbedding/Banksy/Banksy4.output &
# You need to modify the parameter "c" to identify the different number of domains
```
## CellCharter
### DLPFC S1
```python
nohup python Benchmark/RunModel/Run_CellCharter.py \
--input_file Data/DLPFC_sample1/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/DLPFC_sample1/IntergrationRe \
--sample CellCharter --nclust 7 --hvgs 2000 --runNormalization True \
--n_latent 5 --nhood_layers 4 \
> Data/DLPFC_sample1/IntergrationRe/CellCharter.output &
```
### DLPFC S2
```python
nohup python Benchmark/RunModel/Run_CellCharter.py \
--input_file Data/DLPFC_sample2/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/DLPFC_sample2/IntergrationRe \
--sample CellCharter --nclust 5 --hvgs 2000 --runNormalization True \
--n_latent 10 --nhood_layers 2 \
> Data/DLPFC_sample2/IntergrationRe/CellCharter.output &
```
### DLPFC S3
```python
nohup python Benchmark/RunModel/Run_CellCharter.py \
--input_file Data/DLPFC_sample3/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/DLPFC_sample3/IntergrationRe \
--sample CellCharter --nclust 7 --hvgs 2000 --runNormalization True \
--n_latent 15 --nhood_layers 2 \
> Data/DLPFC_sample3/IntergrationRe/CellCharter.output &
```
### MERFISH Preoptic
```python
nohup python Benchmark/RunModel/Run_CellCharter.py \
--input_file Data/MERFISH/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/MERFISH/IntergrationRe \
--sample CellCharter --nclust 8 --runNormalization False \
--n_latent 10 --nhood_layers 4 \
> Data/MERFISH/IntergrationRe/CellCharter.output &
```
### MERFISH Brain S2
```python
nohup python Benchmark/RunModel/Run_CellCharter.py \
--input_file Data/MERFISH_Brain_S2/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/MERFISH_Brain_S2/IntergrationRe \
--sample CellCharter --nclust 8 --hvgs 374 --runNormalization True \
--n_latent 10 --nhood_layers 4 \
> Data/MERFISH_Brain_S2/IntergrationRe/CellCharter.output &
```
### MERFISH Brain S3
```python
nohup python Benchmark/RunModel/Run_CellCharter.py \
--input_file Data/MERFISH_Brain_S3/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/MERFISH_Brain_S3/IntergrationRe \
--sample CellCharter --nclust 8 --hvgs 374 --runNormalization True \
--n_latent 15 --nhood_layers 4 \
> Data/MERFISH_Brain_S3/IntergrationRe/CellCharter.output &
```
### MERFISH Brain S4
```python
nohup python Benchmark/RunModel/Run_CellCharter.py \
--input_file Data/MERFISH_Brain_S4/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/MERFISH_Brain_S4/IntergrationRe \
--sample CellCharter --nclust 8 --hvgs 374 --runNormalization True \
--n_latent 10 --nhood_layers 4 \
> Data/MERFISH_Brain_S4/IntergrationRe/CellCharter.output &
```
### MERFISH Brain S5
```python
nohup python Benchmark/RunModel/Run_CellCharter.py \
--input_file Data/MERFISH_Brain_S5/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/MERFISH_Brain_S5/IntergrationRe \
--sample CellCharter --nclust 8 --hvgs 374 --runNormalization True \
--n_latent 15 --nhood_layers 4 \
> Data/MERFISH_Brain_S5/IntergrationRe/CellCharter.output &
```
### MERFISH Brain S6
```python
nohup python Benchmark/RunModel/Run_CellCharter.py \
--input_file Data/MERFISH_Brain_S6/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/MERFISH_Brain_S6/IntergrationRe \
--sample CellCharter --nclust 8 --hvgs 374 --runNormalization True \
--n_latent 15 --nhood_layers 4 \
> Data/MERFISH_Brain_S6/IntergrationRe/CellCharter.output &
```
### MERFISH Brain S7
```python
nohup python Benchmark/RunModel/Run_CellCharter.py \
--input_file Data/MERFISH_Brain_S7/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/MERFISH_Brain_S7/IntergrationRe \
--sample CellCharter --nclust 8 --hvgs 374 --runNormalization True \
--n_latent 15 --nhood_layers 4 \
> Data/MERFISH_Brain_S7/IntergrationRe/CellCharter.output &
```
### MERFISH Brain S8
```python
nohup python Benchmark/RunModel/Run_CellCharter.py \
--input_file Data/MERFISH_Brain_S8/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/MERFISH_Brain_S8/IntergrationRe \
--sample CellCharter --nclust 8 --hvgs 374 --runNormalization True \
--n_latent 10 --nhood_layers 3 \
> Data/MERFISH_Brain_S8/IntergrationRe/CellCharter.output &
```
### MERFISH Brain S9
```python
nohup python Benchmark/RunModel/Run_CellCharter.py \
--input_file Data/MERFISH_Brain_S9/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/MERFISH_Brain_S9/IntergrationRe \
--sample CellCharter --nclust 8 --hvgs 374 --runNormalization True \
--n_latent 10 --nhood_layers 4 \
> Data/MERFISH_Brain_S9/IntergrationRe/CellCharter.output &
```
### MERFISH Brain S10
```python
nohup python Benchmark/RunModel/Run_CellCharter.py \
--input_file Data/MERFISH_Brain_S10/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/MERFISH_Brain_S10/IntergrationRe \
--sample CellCharter --nclust 8 --hvgs 374 --runNormalization True \
--n_latent 10 --nhood_layers 4 \
> Data/MERFISH_Brain_S10/IntergrationRe/CellCharter.output &
```
### MERFISH Brain S11
```python
nohup python Benchmark/RunModel/Run_CellCharter.py \
--input_file Data/MERFISH_Brain_S11/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/MERFISH_Brain_S11/IntergrationRe \
--sample CellCharter --nclust 8 --hvgs 374 --runNormalization True \
--n_latent 10 --nhood_layers 4 \
> Data/MERFISH_Brain_S11/IntergrationRe/CellCharter.output &
```
### MERFISH Brain S12
```python
nohup python Benchmark/RunModel/Run_CellCharter.py \
--input_file Data/MERFISH_Brain_S12/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/MERFISH_Brain_S12/IntergrationRe \
--sample CellCharter --nclust 8 --hvgs 374 --runNormalization True \
--n_latent 15 --nhood_layers 4 \
> Data/MERFISH_Brain_S12/IntergrationRe/CellCharter.output &
```
### STARMap
```python
nohup python Benchmark/RunModel/Run_CellCharter.py \
--input_file Data/STARMap/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/STARMap/IntergrationRe \
--sample CellCharter --nclust 4 --hvgs 166 --runNormalization True \
--n_latent 5 --nhood_layers 4 \
> Data/STARMap/IntergrationRe/CellCharter.output &
```
### Large-scale Dataset
```python
nohup python Benchmark/RunModel/Run_CellCharter.py \
--input_file Data/Mouse/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/Mouse/IntergrationRe \
--sample CellCharter --nclust 6 --runNormalization False \
--n_latent 15 --nhood_layers 4 \
> Data/Mouse/IntergrationRe/CellCharter.output &
```
### TNBC
```python
nohup python Benchmark/RunModel/Run_CellCharter.py \
--input_file Data/TNBC/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/TNBC/IntergrationRe \
--sample CellCharter --nclust 4 --hvgs 2000 --runNormalization False \
--n_latent 15 --nhood_layers 3 \
> Data/TNBC/IntergrationRe/CellCharter.output &
#  You need to modify the parameter "nclust" to identify the different number of domains
```
