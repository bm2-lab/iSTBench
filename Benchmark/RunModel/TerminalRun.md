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
