# Alignment of BaristaSeq
## Rotation
```python
nohup python Benchmark/Alignment/Rotate_spatial.py \
--input_file Data/BaristaSeq/sample_data \
--output_file Data/BaristaSeq/sample_data_rotate \
--batches "slices1,slices2,slices3" \
--angles 40 20 0 \
> Data/BaristaSeq/rotate.output &
```
## Alignment
```python
# Banksy
nohup python Benchmark/Alignment/Run_STAligner.py \
--input_file Data/BaristaSeq/sample_data_rotate \
--input_data Data/BaristaSeq/IntergrationRe/Banksy.h5ad \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" \
--landmark_domain 5 --landmark_domain_original VISp_wm --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/BaristaSeq/Banksy.output &

# CellCharter
nohup python Benchmark/Alignment/Run_STAligner.py \
--input_file Data/BaristaSeq/sample_data_rotate \
--input_data Data/BaristaSeq/IntergrationRe/CellCharter.h5ad \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" \
--landmark_domain 5 --landmark_domain_original VISp_wm --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/BaristaSeq/CellCharter.output &

# CN
nohup python Benchmark/Alignment/Run_STAligner.py \
--input_file Data/BaristaSeq/sample_data_rotate \
--input_data Data/BaristaSeq/IntergrationRe/CN.h5ad \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" \
--landmark_domain 2 --landmark_domain_original VISp_wm --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/BaristaSeq/CN.output &

# GraphST
nohup python Benchmark/Alignment/Run_STAligner.py \
--input_file Data/BaristaSeq/sample_data_rotate \
--input_data Data/BaristaSeq/IntergrationRe/GraphST.h5ad \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" \
--landmark_domain 6 --landmark_domain_original VISp_wm --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/BaristaSeq/GraphST.output &

# GraphST-PASTE
nohup python Benchmark/Alignment/Run_STAligner.py \
--input_file Data/BaristaSeq/sample_data_rotate \
--input_data Data/BaristaSeq/IntergrationRe/GraphSTwithPASTE.h5ad \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" \
--landmark_domain 4 --landmark_domain_original VISp_wm --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/BaristaSeq/GraphSTwithPASTE.output &

# NicheCompass
nohup python Benchmark/Alignment/Run_STAligner.py \
--input_file Data/BaristaSeq/sample_data_rotate \
--input_data Data/BaristaSeq/IntergrationRe/NicheCompass.h5ad \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" \
--landmark_domain 3 --landmark_domain_original VISp_wm --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/BaristaSeq/NicheCompass.output &

# MENDER
nohup python Benchmark/Alignment/Run_STAligner.py \
--input_file Data/BaristaSeq/sample_data_rotate \
--input_data Data/BaristaSeq/IntergrationRe/MENDER.h5ad \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" \
--landmark_domain 4 --landmark_domain_original VISp_wm --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/BaristaSeq/MENDER.output &

# Spado
nohup python Benchmark/Alignment/Run_STAligner.py \
--input_file Data/BaristaSeq/sample_data_rotate \
--input_data Data/BaristaSeq/IntergrationRe/Spado.h5ad \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" \
--landmark_domain "Domain_6_4" --landmark_domain_original VISp_wm --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/BaristaSeq/Spado.output &
```
# Alignment of other datasets
## Rotation
```python
# DLPFC S1
nohup python Benchmark/Alignment/Rotate_spatial.py \
--input_file Data/DLPFC_sample1/sample_data \
--output_file Data/DLPFC_sample1/sample_data_rotate \
--batches "slices2,slices3,slices1,slices4" \
--angles 60 40 20 0 \
> Data/DLPFC_sample1/rotate.output &

# DLPFC S2
nohup python Benchmark/Alignment/Rotate_spatial.py \
--input_file Data/DLPFC_sample2/sample_data \
--output_file Data/DLPFC_sample2/sample_data_rotate \
--batches "slices2,slices1,slices3,slices4" \
--angles 60 40 20 0 \
> Data/DLPFC_sample2/rotate.output &

# DLPFC S3
nohup python Benchmark/Alignment/Rotate_spatial.py \
--input_file Data/DLPFC_sample3/sample_data \
--output_file Data/DLPFC_sample3/sample_data_rotate \
--batches "slices2,slices4,slices1,slices3" \
--angles 60 40 20 0 \
> Data/DLPFC_sample3/rotate.output &

# MERFISH Preoptic
nohup python Benchmark/Alignment/Rotate_spatial.py \
--input_file Data/MERFISH/sample_data \
--output_file Data/MERFISH/sample_data_rotate \
--batches "slices1,slices4,slices2,slices5,slices3" \
--angles 80 60 40 20 0 \
> Data/MERFISH/rotate.output &

# MERFISH Brain S3
nohup python Benchmark/Alignment/Rotate_spatial.py \
--input_file Data/MERFISH_Brain_S3/sample_data \
--output_file Data/MERFISH_Brain_S3/sample_data_rotate \
--batches "slices2,slices1" \
--angles 20 0 \
> Data/MERFISH_Brain_S3/rotate.output &

# STARMap
nohup python Benchmark/Alignment/Rotate_spatial.py \
--input_file Data/STARMap/sample_data \
--output_file Data/STARMap/sample_data_rotate \
--batches "slices1,slices3,slices2" \
--angles 0 20 40 \
> Data/STARMap/rotate.output &
```
## Alignment
See TTerminalRun2.txt





