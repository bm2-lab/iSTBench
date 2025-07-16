# Run methods on BaristaSeq dataset
## Banksy
```python
nohup Rscript Benchmark/RunModel/Run_Banksy.R \
-i Data/BaristaSeq/sample_all_data/Slices_combind_data.RDS \
-o Data/BaristaSeq/IntergrationRe \
-s Banksy -c 6 -d RDS -k 18 -l 0.2 \
> Data/BaristaSeq/IntergrationRe/Banksy.output &
```  
## CellCharter
```python
nohup python Benchmark/RunModel/Run_CellCharter.py \
--input_file Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/BaristaSeq/IntergrationRe \
--sample CellCharter --nclust 6 --runNormalization False \
--n_latent 15 --nhood_layers 4 \
> Data/BaristaSeq/IntergrationRe/CellCharter.output &
```
## CN
```python
nohup python Benchmark/RunModel/Run_CN.py \
--input_file Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/BaristaSeq/IntergrationRe \
--sample CN --nclust 6 --runNormalization False \
> Data/BaristaSeq/IntergrationRe/CN.output &
```
## GraphST
```python
nohup python Benchmark/RunModel/GraphST/Run_GraphST.py \
--input_file Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/BaristaSeq/IntergrationRe \
--sample GraphST --nclust 6 --runNormalization False --device cuda \
> Data/BaristaSeq/IntergrationRe/GraphST.output  &
```
## GraphST-PASTE
```python
nohup python Benchmark/RunModel/GraphST/Run_GraphST-PASTE.py \
--input_file Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/BaristaSeq/IntergrationRe \
--sample GraphSTwithPASTE --nclust 6 --runNormalization False --device cuda \
> Data/BaristaSeq/IntergrationRe/GraphSTwithPASTE.output  &
```
## MENDER
```python
nohup python Benchmark/RunModel/Run_MENDER.py \
--input_file Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/BaristaSeq/IntergrationRe \
--sample MENDER --nclust 6 --runNormalization False --tech BaristaSeq \
> Data/BaristaSeq/IntergrationRe/MENDER.output &
```
## NicheCompass 
```python
nohup python Benchmark/RunModel/Run_NicheCompass.py \
--input_file Data/BaristaSeq/sample_data \
--output_file Data/BaristaSeq/IntergrationRe \
--sample NicheCompass --nclust 6 \
> Data/BaristaSeq/IntergrationRe/NicheCompass.output &
```
## PRECAST 
```python
nohup Rscript Benchmark/RunModel/Run_PRECAST.R \
-i Data/BaristaSeq/sample_data  \
-o  Data/BaristaSeq/IntergrationRe \
-s PRECAST -c 6 -n FALSE -t Other_SRT  -p Mouse  \
> Data/BaristaSeq/IntergrationRe/PRECAST.output &    
```
## SpaDo
```python
nohup /opt/R/4.3.2/lib/R/bin/Rscript Benchmark/RunModel/Run_Spado.R \
-i Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
-o Data/BaristaSeq/IntergrationRe \
-s Spado -c 6 -n FALSE \
> Data/BaristaSeq/IntergrationRe/Spado.output &    
```
## SPIRAL
```python
nohup python Benchmark/RunModel/Run_SPIRAL.py \
--input_file Data/BaristaSeq/sample_data \
--output_file Data/BaristaSeq/IntergrationRe \
--model_path  Data/BaristaSeq/IntergrationRe/SPIRAL_intermediate_result \
--sample SPIRAL --runNormalization False --data BaristaSeq --nclust 6 \
> Data/BaristaSeq/IntergrationRe/SPIRAL.output & 
```
## STAIG
```python
nohup python Benchmark/RunModel/Run_STAIG.py \
--input_file Data/BaristaSeq/sample_data \
--output_file Data/BaristaSeq/IntergrationRe \
--sample STAIG --runNormalization False --nclust 6 \
--nPCA 50 --num_layers 1 --k_neighbor 5 --tau 10 --num_epochs 400 \
> Data/BaristaSeq/IntergrationRe/STAIG.output &
```
## STAligner
```python
nohup python Benchmark/RunModel/Run_STAligner.py \
--input_file Data/BaristaSeq/sample_data \
--output_file Data/BaristaSeq/IntergrationRe \
--sample STAligner --batches "slices1,slices2,slices3" \
--runNormalization False --nclust 6 --r 50 \
> Data/BaristaSeq/IntergrationRe/STAligner.output &
```


