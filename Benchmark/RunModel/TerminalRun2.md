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
-o Data/TNBC/SlicesEmbedding \
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
--output_file Data/TNBC/SlicesEmbedding \
--sample CellCharter --nclust 4 --hvgs 2000 --runNormalization False \
--n_latent 15 --nhood_layers 3 \
> Data/TNBC/SlicesEmbedding/CellCharter/CellCharter4.output &
#  You need to modify the parameter "nclust" to identify the different number of domains
```
## CN
### DLPFC S1
```python
nohup python Benchmark/RunModel/Run_CN.py \
--input_file Data/DLPFC_sample1/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/DLPFC_sample1/IntergrationRe \
--sample CN --nclust 7 \
> Data/DLPFC_sample1/IntergrationRe/CN.output &
```
### DLPFC S2
```python
nohup python Benchmark/RunModel/Run_CN.py \
--input_file Data/DLPFC_sample2/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/DLPFC_sample2/IntergrationRe \
--sample CN --nclust 5 \
> Data/DLPFC_sample2/IntergrationRe/CN.output &
```
### DLPFC S3
```python
nohup python Benchmark/RunModel/Run_CN.py \
--input_file Data/DLPFC_sample3/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/DLPFC_sample3/IntergrationRe \
--sample CN --nclust 7 \
> Data/DLPFC_sample3/IntergrationRe/CN.output &
```
### MERFISH Preoptic
```python
nohup python Benchmark/RunModel/Run_CN.py \
--input_file Data/MERFISH/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/MERFISH/IntergrationRe \
--sample CN --nclust 8 --runNormalization False \
> Data/MERFISH/IntergrationRe/CN.output &
```
### MERFISH Brain S2-S12
```python
nohup python Benchmark/RunModel/Run_CN.py \
--input_file Data/MERFISH_Brain_S2/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/MERFISH_Brain_S2/IntergrationRe \
--sample CN --nclust 8 --runNormalization True --hvg 374 \
> Data/MERFISH_Brain_S2/IntergrationRe/CN.output &
# You just need to change the path, like MERFISH_Brain_S3, MERFISH_Brain_S4...
```
### STARMap
```python
nohup python Benchmark/RunModel/Run_CN.py \
--input_file Data/STARMap/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/STARMap/IntergrationRe \
--sample CN --nclust 4 --runNormalization True --hvgs 166 \
> Data/STARMap/IntergrationRe/CN.output &
```
### Large-scale Dataset
```python
nohup python Benchmark/RunModel/Run_CN.py \
--input_file Data/Mouse/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/Mouse/IntergrationRe \
--sample CN --nclust 7 --runNormalization False \
> Data/Mouse/IntergrationRe/CN.output &
```
### TNBC
```python
nohup python Benchmark/RunModel/Run_CN.py \
--input_file Data/TNBC/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/TNBC/SlicesEmbedding \
--sample CN --nclust 4 \
> Data/TNBC/SlicesEmbedding/CN/CN4.output &
#  You need to modify the parameter "nclust" to identify the different number of domains
```
## GraphST
### DLPFC S1
```python
nohup python Benchmark/RunModel/GraphST/Run_GraphST.py \
--input_file Data/DLPFC_sample1/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/DLPFC_sample1/IntergrationRe \
--sample GraphST --nclust 7 --device cuda \
> Data/DLPFC_sample1/IntergrationRe/GraphST.output  &
```
### DLPFC S2
```python
nohup python Benchmark/RunModel/GraphST/Run_GraphST.py \
--input_file Data/DLPFC_sample2/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/DLPFC_sample2/IntergrationRe \
--sample GraphST --nclust 5 --device cuda \
> Data/DLPFC_sample2/IntergrationRe/GraphST.output  &
```
### DLPFC S3
```python
nohup python Benchmark/RunModel/GraphST/Run_GraphST.py \
--input_file Data/DLPFC_sample3/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/DLPFC_sample3/IntergrationRe \
--sample GraphST --nclust 7 --device cuda \
> Data/DLPFC_sample3/IntergrationRe/GraphST.output  &
```
### MERFISH Preoptic
```python
nohup python Benchmark/RunModel/GraphST/Run_GraphST.py \
--input_file Data/MERFISH/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/MERFISH/IntergrationRe \
--sample GraphST --nclust 8 --runNormalization False --device cuda \
> Data/MERFISH/IntergrationRe/GraphST.output  &
```
### MERFISH Brain S2-S12
```python
nohup python Benchmark/RunModel/GraphST/Run_GraphST.py \
--input_file Data/MERFISH_Brain_S2/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/MERFISH_Brain_S2/IntergrationRe \
--sample GraphST --nclust 8 --hvgs 374 --device cuda \
> Data/MERFISH_Brain_S2/IntergrationRe/GraphST.output  &
# You just need to change the path, like MERFISH_Brain_S3, MERFISH_Brain_S4...
```
### STARMap
```python
nohup python Benchmark/RunModel/GraphST/Run_GraphST.py \
--input_file Data/STARMap/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/STARMap/IntergrationRe \
--sample GraphST --nclust 4 --hvgs 166 --device cuda \
> Data/STARMap/IntergrationRe/GraphST.output  &
```
### Large-scale Dataset
```python
nohup python Benchmark/RunModel/GraphST/Run_GraphST.py \
--input_file Data/Mouse/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/Mouse/IntergrationRe \
--sample GraphST --nclust 7 --runNormalization False --device cuda \
> Data/Mouse/IntergrationRe/GraphST.output  &
```
## GraphST-PASTE
### DLPFC S1
```python
nohup python Benchmark/RunModel/GraphST/Run_GraphST-PASTE.py \
--input_file Data/DLPFC_sample1/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/DLPFC_sample1/IntergrationRe \
--sample GraphSTwithPASTE --nclust 7 --device cuda \
> Data/DLPFC_sample1/IntergrationRe/GraphST.output  &
```
### DLPFC S2
```python
nohup python Benchmark/RunModel/GraphST/Run_GraphST-PASTE.py \
--input_file Data/DLPFC_sample2/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/DLPFC_sample2/IntergrationRe \
--sample GraphSTwithPASTE --nclust 5 --device cuda \
> Data/DLPFC_sample2/IntergrationRe/GraphST.output  &
```
### DLPFC S3
```python
nohup python Benchmark/RunModel/GraphST/Run_GraphST-PASTE.py \
--input_file Data/DLPFC_sample3/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/DLPFC_sample3/IntergrationRe \
--sample GraphSTwithPASTE --nclust 7 --device cuda \
> Data/DLPFC_sample3/IntergrationRe/GraphST.output  &
```
### MERFISH Preoptic
```python
nohup python Benchmark/RunModel/GraphST/Run_GraphST-PASTE.py \
--input_file Data/MERFISH/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/MERFISH/IntergrationRe \
--sample GraphSTwithPASTE --nclust 8 --runNormalization False --device cuda \
> Data/MERFISH/IntergrationRe/GraphST.output  &
```
### MERFISH Brain S2-S12
```python
nohup python Benchmark/RunModel/GraphST/Run_GraphST-PASTE.py \
--input_file Data/MERFISH_Brain_S2/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/MERFISH_Brain_S2/IntergrationRe \
--sample GraphSTwithPASTE --nclust 8 --hvgs 374 --device cuda \
> Data/MERFISH_Brain_S2/IntergrationRe/GraphST.output  &
# You just need to change the path, like MERFISH_Brain_S3, MERFISH_Brain_S4...
```
### STARMap
```python
nohup python Benchmark/RunModel/GraphST/Run_GraphST-PASTE.py \
--input_file Data/STARMap/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/STARMap/IntergrationRe \
--sample GraphSTwithPASTE --nclust 4 --hvgs 166 --device cuda \
> Data/STARMap/IntergrationRe/GraphST.output  &
```
### Large-scale Dataset
```python
nohup python Benchmark/RunModel/GraphST/Run_GraphST-PASTE.py \
--input_file Data/Mouse/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/Mouse/IntergrationRe \
--sample GraphSTwithPASTE --nclust 7 --runNormalization False --device cuda \
> Data/Mouse/IntergrationRe/GraphST.output  &
```
## MENDER
### DLPFC S1
```python
nohup python Benchmark/RunModel/Run_MENDER.py \
--input_file Data/DLPFC_sample1/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/DLPFC_sample1/IntergrationRe \
--sample MENDER --nclust 7 --tech Visium \
> Data/DLPFC_sample1/IntergrationRe/MENDER.output &
```
### DLPFC S2
```python
nohup python Benchmark/RunModel/Run_MENDER.py \
--input_file Data/DLPFC_sample2/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/DLPFC_sample2/IntergrationRe \
--sample MENDER --nclust 5 --tech Visium \
> Data/DLPFC_sample2/IntergrationRe/MENDER.output &
```
### DLPFC S3
```python
nohup python Benchmark/RunModel/Run_MENDER.py \
--input_file Data/DLPFC_sample3/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/DLPFC_sample3/IntergrationRe \
--sample MENDER --nclust 7 --tech Visium \
> Data/DLPFC_sample3/IntergrationRe/MENDER.output &
```
### MERFISH Preoptic
```python
nohup python Benchmark/RunModel/Run_MENDER.py \
--input_file Data/MERFISH/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/MERFISH/IntergrationRe \
--sample MENDER --nclust 6 --runNormalization False --tech MERFISH \
> Data/MERFISH/IntergrationRe/MENDER.output &
```
### MERFISH Brain S2-S12
```python
nohup python Benchmark/RunModel/Run_MENDER.py \
--input_file Data/MERFISH_Brain_S2/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/MERFISH_Brain_S2/IntergrationRe \
--sample MENDER --nclust 8 --tech MERFISH --hvgs 374 \
> Data/MERFISH_Brain_S2/IntergrationRe/MENDER.output &
# You just need to change the path, like MERFISH_Brain_S3, MERFISH_Brain_S4...
```
### STARMap
```python
nohup python Benchmark/RunModel/Run_MENDER.py \
--input_file Data/STARMap/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/STARMap/IntergrationRe \
--sample MENDER --nclust 4 --hvgs 166 --tech STARMap \
> Data/STARMap/IntergrationRe/MENDER.output &
```
### Large-scale Dataset
```python
nohup python Benchmark/RunModel/Run_MENDER.py \
--input_file Data/Mouse/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/Mouse/IntergrationRe \
--sample MENDER --nclust 7 --runNormalization False --tech Unknown \
> Data/Mouse/IntergrationRe/MENDER.output &
```
### TNBC
```python
nohup python Benchmark/RunModel/Run_MENDER.py \
--input_file Data/TNBC/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/TNBC/SlicesEmbedding \
--sample MENDER --nclust 4 --runNormalization False --tech MIBI \
> Data/TNBC/SlicesEmbedding/MENDER/MENDER4.output &
# You need to modify the parameter "nclust" to identify the different number of domains
```
## NicheCompass
### DLPFC S1
```python
nohup python Benchmark/RunModel/Run_NicheCompass.py \
--input_file Data/DLPFC_sample1/sample_data \
--output_file Data/DLPFC_sample1/IntergrationRe \
--sample NicheCompass --nclust 7 \
> Data/DLPFC_sample1/IntergrationRe/NicheCompass.output &
```
### DLPFC S2
```python
nohup python Benchmark/RunModel/Run_NicheCompass.py \
--input_file Data/DLPFC_sample2/sample_data \
--output_file Data/DLPFC_sample2/IntergrationRe \
--sample NicheCompass --nclust 5 \
> Data/DLPFC_sample2/IntergrationRe/NicheCompass.output &
```
### DLPFC S3
```python
nohup python Benchmark/RunModel/Run_NicheCompass.py \
--input_file Data/DLPFC_sample3/sample_data \
--output_file Data/DLPFC_sample3/IntergrationRe \
--sample NicheCompass --nclust 7 \
> Data/DLPFC_sample3/IntergrationRe/NicheCompass.output &
```
### MERFISH Preoptic
```python
nohup python Benchmark/RunModel/Run_NicheCompass.py \
--input_file Data/MERFISH/sample_data \
--output_file Data/MERFISH/IntergrationRe \
--sample NicheCompass --nclust 8 --species mouse \
> Data/MERFISH/IntergrationRe/NicheCompass.output &
```
### MERFISH Brain S2-S12
```python
nohup python Benchmark/RunModel/Run_NicheCompass.py \
--input_file Data/MERFISH_Brain_S2/sample_data \
--output_file Data/MERFISH_Brain_S2/IntergrationRe \
--sample NicheCompass --nclust 8 --species mouse --convertID True \
> Data/MERFISH_Brain_S2/IntergrationRe/NicheCompass.output &
# You just need to change the path, like MERFISH_Brain_S3, MERFISH_Brain_S4...
```
### STARMap
```python
nohup python Benchmark/RunModel/Run_NicheCompass.py \
--input_file Data/STARMap/sample_data \
--output_file Data/STARMap/IntergrationRe \
--sample NicheCompass --nclust 4 --species mouse \
> Data/STARMap/IntergrationRe/NicheCompass.output &
```
### Large-scale Dataset
```python
nohup python Benchmark/RunModel/Run_NicheCompass.py \
--input_file Data/Mouse/sample_data \
--output_file Data/Mouse/IntergrationRe \
--sample NicheCompass --nclust 7 --species mouse \
> Data/Mouse/IntergrationRe/NicheCompass.output &
```
### TNBC
```python
nohup python Benchmark/RunModel/Run_NicheCompass.py \
--input_file Data/TNBC/sample_data \
--output_file Data/TNBC/SlicesEmbedding \
--sample NicheCompass --nclust 4 \
> Data/TNBC/SlicesEmbedding/NicheCompass/NicheCompass4.output &
# You need to modify the parameter "nclust" to identify the different number of domains
```
## Spado
### DLPFC S1
```python
nohup Rscript Benchmark/RunModel/Run_Spado.R \
-i Data/DLPFC_sample1/sample_all_data/Slices_combind_data.h5ad \
-o Data/DLPFC_sample1/IntergrationRe \
-s Spado -v 2000 -c 7 -n TRUE \
> Data/DLPFC_sample1/IntergrationRe/Spado.output &
```
### DLPFC S2
```python
nohup Rscript Benchmark/RunModel/Run_Spado.R \
-i Data/DLPFC_sample2/sample_all_data/Slices_combind_data.h5ad \
-o Data/DLPFC_sample2/IntergrationRe \
-s Spado -v 2000 -c 5 -n TRUE \
> Data/DLPFC_sample2/IntergrationRe/Spado.output &
```
### DLPFC S3
```python
nohup Rscript Benchmark/RunModel/Run_Spado.R \
-i Data/DLPFC_sample3/sample_all_data/Slices_combind_data.h5ad \
-o Data/DLPFC_sample3/IntergrationRe \
-s Spado -v 2000 -c 7 -n TRUE \
> Data/DLPFC_sample3/IntergrationRe/Spado.output &
``` 
### MERFISH Preoptic
```python
nohup Rscript Benchmark/RunModel/Run_Spado.R \
-i Data/MERFISH/sample_all_data/Slices_combind_data.h5ad \
-o Data/MERFISH/IntergrationRe \
-s Spado -c 8 -n FALSE \
> Data/MERFISH/IntergrationRe/Spado.output &
``` 
### MERFISH Brain S2-S12
```python
nohup Rscript Benchmark/RunModel/Run_Spado.R \
-i Data/MERFISH_Brain_S2/sample_all_data/Slices_combind_data.h5ad \
-o Data/MERFISH_Brain_S2/IntergrationRe \
-s Spado -v 374 -c 8 -n TRUE \
> Data/MERFISH_Brain_S2/IntergrationRe/Spado.output &
# You just need to change the path, like MERFISH_Brain_S3, MERFISH_Brain_S4...
```
### STARMap
```python
nohup Rscript Benchmark/RunModel/Run_Spado.R \
-i Data/STARMap/sample_all_data/Slices_combind_data.h5ad \
-o Data/STARMap/IntergrationRe \
-s Spado -v 166 -c 4 -n TRUE \
> Data/STARMap/IntergrationRe/Spado.output &
```
### Large-scale Dataset
```python
nohup Rscript Benchmark/RunModel/Run_Spado.R \
-i Data/Mouse/sample_all_data/Slices_combind_data.h5ad \
-o Data/Mouse/IntergrationRe \
-s Spado -c 7 -n FALSE \
> Data/Mouse/IntergrationRe/Spado.output &
``` 

