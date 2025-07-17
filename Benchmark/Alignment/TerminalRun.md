# PASTE
## DLPFC S1
```python
nohup python Benchmark/Alignment/Run_PASTE.py \
--input_file Data/DLPFC_sample1/sample_data \
--output_file Benchmark/Alignment/Result/DLPFC_sample1 \
--batches "slices2,slices3,slices1,slices4" --step 10 \
> Benchmark/Alignment/Result/DLPFC_sample1/PASTE.output &
```
## DLPFC S2
```python
nohup python Benchmark/Alignment/Run_PASTE.py \
--input_file Data/DLPFC_sample2/sample_data \
--output_file Benchmark/Alignment/Result/DLPFC_sample2 \
--batches "slices2,slices3,slices1,slices4" --step 10 \
> Benchmark/Alignment/Result/DLPFC_sample2/PASTE.output &
```
## DLPFC S3
```python
nohup python Benchmark/Alignment/Run_PASTE.py \
--input_file Data/DLPFC_sample3/sample_data \
--output_file Benchmark/Alignment/Result/DLPFC_sample3 \
--batches "slices2,slices4,slices1,slices3" --step 10 \
> Benchmark/Alignment/Result/DLPFC_sample3/PASTE.output &
```
## MERFISH Preoptic
```python
nohup python Benchmark/Alignment/Run_PASTE.py \
--input_file Data/MERFISH/sample_data \
--output_file Benchmark/Alignment/Result/MERFISH \
--batches "slices1,slices4,slices2,slices5,slices3" --step 10 \
> Benchmark/Alignment/Result/MERFISH/PASTE.output &
```
## MERFISH Brain S3
```python
nohup python Benchmark/Alignment/Run_PASTE.py \
--input_file Data/MERFISH_Brain_S3/sample_data \
--output_file Benchmark/Alignment/Result/MERFISH_Brain_S3 \
--batches "slices2,slices1" --step 10 \
> Benchmark/Alignment/Result/MERFISH_Brain_S3/PASTE.output &
```
## STARMap
```python
nohup python Benchmark/Alignment/Run_PASTE.py \
--input_file Data/STARMap/sample_data \
--output_file Benchmark/Alignment/Result/STARMap \
--batches "slices1,slices3,slices2" --step 10 \
> Benchmark/Alignment/Result/STARMap/PASTE.output &
```
# STalign
## MERFISH Preoptic
```python
nohup python Benchmark/Alignment/Run_STalign.py \
--input_file Data/MERFISH/sample_data \
--output_file Benchmark/Alignment/Result/MERFISH \
--batches "slices1,slices4,slices2,slices5,slices3" --step 10 \
> Benchmark/Alignment/Result/MERFISH/STalign.output &
```
## MERFISH Brain S3
```python
nohup python Benchmark/Alignment/Run_STalign.py \
--input_file Data/MERFISH_Brain_S3/sample_data \
--output_file Benchmark/Alignment/Result/MERFISH_Brain_S3 \
--batches "slices2,slices1" --step 10 \
> Benchmark/Alignment/Result/MERFISH_Brain_S3/STalign.output &
```
## STARMap
```python
nohup python Benchmark/Alignment/Run_STalign.py \
--input_file Data/STARMap/sample_data \
--output_file Benchmark/Alignment/Result/STARMap \
--batches "slices1,slices3,slices2" --step 10 \
> Benchmark/Alignment/Result/STARMap/STalign.output &
```
# SPACEL
## DLPFC S1
```python
nohup python Benchmark/Alignment/Run_SPACEL.py \
--input_file Data/DLPFC_sample1/sample_data \
--input_data Data/DLPFC_sample1/IntergrationRe/MENDER.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample1 \
--batches "slices2,slices3,slices1,slices4" --step 10 \
> Benchmark/Alignment/Result/DLPFC_sample1/SPACEL_MENDER.output &
# You just need to change imput_path to use the results of different methods
```
## DLPFC S2
```python
nohup python Benchmark/Alignment/Run_SPACEL.py \
--input_file Data/DLPFC_sample2/sample_data \
--input_data Data/DLPFC_sample2/IntergrationRe/MENDER.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample2 \
--batches "slices2,slices3,slices1,slices4" --step 10 \
> Benchmark/Alignment/Result/DLPFC_sample2/SPACEL_MENDER.output &
# You just need to change imput_path to use the results of different methods
```
## DLPFC S3
```python
nohup python Benchmark/Alignment/Run_SPACEL.py \
--input_file Data/DLPFC_sample3/sample_data \
--input_data Data/DLPFC_sample3/IntergrationRe/MENDER.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample3 \
--batches "slices2,slices4,slices1,slices3" --step 10 \
> Benchmark/Alignment/Result/DLPFC_sample3/SPACEL_MENDER.output &
# You just need to change imput_path to use the results of different methods
```
## MERFISH Preoptic
```python
nohup python Benchmark/Alignment/Run_SPACEL.py \
--input_file Data/MERFISH/sample_data \
--input_data Data/MERFISH/IntergrationRe/MENDER.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH \
--batches "slices1,slices4,slices2,slices5,slices3" --step 10 \
> Benchmark/Alignment/Result/MERFISH/SPACEL_MENDER.output &
# You just need to change imput_path to use the results of different methods
```
## MERFISH Brain S3
```python
nohup python Benchmark/Alignment/Run_SPACEL.py \
--input_file Data/MERFISH_Brain_S3/sample_data \
--input_data Data/MERFISH_Brain_S3/IntergrationRe/MENDER.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH_Brain_S3 \
--batches "slices2,slices1" --step 10 \
> Benchmark/Alignment/Result/MERFISH_Brain_S3/SPACEL_MENDER.output &
# You just need to change imput_path to use the results of different methods
```
## STARMap
```python
nohup python Benchmark/Alignment/Run_SPACEL.py \
--input_file Data/STARMap/sample_data \
--input_data Data/STARMap/IntergrationRe/MENDER.h5ad \
--output_file Benchmark/Alignment/Result/STARMap \
--batches "slices1,slices3,slices2" --step 10 \
> Benchmark/Alignment/Result/STARMap/SPACEL_MENDER.output &
# You just need to change imput_path to use the results of different methods
```
# STAligner
## BaristaSeq
### Banksy
```python
nohup python Run_STAligner.py \
--input_file Data/BaristaSeq/sample_data \
--input_data Data/BaristaSeq/IntergrationRe/Banksy.h5ad \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" \
--landmark_domain 5 --landmark_domain_original VISp_wm --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/BaristaSeq/STAligner_Banksy.output &
```
### CellCharter
```python
nohup python Run_STAligner.py \
--input_file Data/BaristaSeq/sample_data \
--input_data Data/BaristaSeq/IntergrationRe/CellCharter.h5ad \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" \
--landmark_domain 5 --landmark_domain_original VISp_wm --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/BaristaSeq/STAligner_CellCharter.output &
```
### CN
```python
nohup python Run_STAligner.py \
--input_file Data/BaristaSeq/sample_data \
--input_data Data/BaristaSeq/IntergrationRe/CN.h5ad \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" \
--landmark_domain 3 --landmark_domain_original VISp_wm --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/BaristaSeq/STAligner_CN.output &
```
### GraphST
```python
nohup python Run_STAligner.py \
--input_file Data/BaristaSeq/sample_data \
--input_data Data/BaristaSeq/IntergrationRe/GraphST.h5ad \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" \
--landmark_domain 6 --landmark_domain_original VISp_wm --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/BaristaSeq/STAligner_GraphST.output &
```
### MENDER
```python
nohup python Run_STAligner.py \
--input_file Data/BaristaSeq/sample_data \
--input_data Data/BaristaSeq/IntergrationRe/MENDER.h5ad \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" \
--landmark_domain 4 --landmark_domain_original VISp_wm --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/BaristaSeq/STAligner_MENDER.output &
```
### NicheCompass
```python
nohup python Run_STAligner.py \
--input_file Data/BaristaSeq/sample_data \
--input_data Data/BaristaSeq/IntergrationRe/NicheCompass.h5ad \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" \
--landmark_domain 4 --landmark_domain_original VISp_wm --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/BaristaSeq/STAligner_NicheCompass.output &
```
### PRECAST
```python
nohup python Run_STAligner.py \
--input_file Data/BaristaSeq/sample_data \
--input_data Data/BaristaSeq/IntergrationRe/PRECAST.h5ad \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" \
--landmark_domain 4 --landmark_domain_original VISp_wm --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/BaristaSeq/STAligner_PRECAST.output &
```
### SpaDo
```python
nohup python Run_STAligner.py \
--input_file Data/BaristaSeq/sample_data \
--input_data Data/BaristaSeq/IntergrationRe/Spado.h5ad \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" \
--landmark_domain "Domain_6_4" --landmark_domain_original VISp_wm --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/BaristaSeq/STAligner_Spado.output &
```
### SPIRAL
```python
nohup python Run_STAligner.py \
--input_file Data/BaristaSeq/sample_data \
--input_data Data/BaristaSeq/IntergrationRe/SPIRAL.h5ad \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" \
--landmark_domain 4 --landmark_domain_original VISp_wm --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/BaristaSeq/STAligner_SPIRAL.output &
```
### STAIG
```python
nohup python Run_STAligner.py \
--input_file Data/BaristaSeq/sample_data \
--input_data Data/BaristaSeq/IntergrationRe/STAIG.h5ad \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" \
--landmark_domain 5 --landmark_domain_original VISp_wm --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/BaristaSeq/STAligner_STAIG.output &
```
### STAligner
```python
nohup python Run_STAligner.py \
--input_file Data/BaristaSeq/sample_data \
--input_data Data/BaristaSeq/IntergrationRe/STAligner.h5ad \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" \
--landmark_domain 1 --landmark_domain_original VISp_wm --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/BaristaSeq/STAligner_STAligner.output &
```
## DLPFC S1
### Banksy
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample1/sample_data \
--input_data Data/DLPFC_sample1/IntergrationRe/Banksy.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample1 \
--batches "slices2,slices3,slices1,slices4" \
--landmark_domain 6 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample1/STAligner_Banksy.output &
```
### CellCharter
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample1/sample_data \
--input_data Data/DLPFC_sample1/IntergrationRe/CellCharter.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample1 \
--batches "slices2,slices3,slices1,slices4" \
--landmark_domain 6 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample1/STAligner_CellCharter.output &
```
### CN
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample1/sample_data \
--input_data Data/DLPFC_sample1/IntergrationRe/CN.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample1 \
--batches "slices2,slices3,slices1,slices4" \
--landmark_domain 3 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample1/STAligner_CN.output &
```
### MENDER
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample1/sample_data \
--input_data Data/DLPFC_sample1/IntergrationRe/MENDER.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample1 \
--batches "slices2,slices3,slices1,slices4" \
--landmark_domain 5 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample1/STAligner_MENDER.output &
```
### PRECAST
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample1/sample_data \
--input_data Data/DLPFC_sample1/IntergrationRe/PRECAST.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample1 \
--batches "slices2,slices3,slices1,slices4" \
--landmark_domain 1 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample1/STAligner_PRECAST.output &
```
### Spado
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample1/sample_data \
--input_data Data/DLPFC_sample1/IntergrationRe/Spado.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample1 \
--batches "slices2,slices3,slices1,slices4" \
--landmark_domain "Domain_7_7" --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample1/STAligner_Spado.output &
```
### SPIRAL
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample1/sample_data \
--input_data Data/DLPFC_sample1/IntergrationRe/SPIRAL.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample1 \
--batches "slices2,slices3,slices1,slices4" \
--landmark_domain 6 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample1/STAligner_SPIRAL.output &
```
### STAIG
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample1/sample_data \
--input_data Data/DLPFC_sample1/IntergrationRe/STAIG.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample1 \
--batches "slices2,slices3,slices1,slices4" \
--landmark_domain 2 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample1/STAligner_STAIG.output &
```
### STAligner
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample1/sample_data \
--input_data Data/DLPFC_sample1/IntergrationRe/STAligner.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample1 \
--batches "slices2,slices3,slices1,slices4" \
--landmark_domain 6 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample1/STAligner_STAligner.output &
```
## DLPFC S2
### Banksy
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample2/sample_data \
--input_data Data/DLPFC_sample2/IntergrationRe/Banksy.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample2 \
--batches "slices2,slices3,slices1,slices4" \
--landmark_domain 5 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample2/STAligner_Banksy.output &
```
### CellCharter
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample2/sample_data \
--input_data Data/DLPFC_sample2/IntergrationRe/CellCharter.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample2 \
--batches "slices2,slices3,slices1,slices4" \
--landmark_domain 2 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample2/STAligner_CellCharter.output &
```
### CN
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample2/sample_data \
--input_data Data/DLPFC_sample2/IntergrationRe/CN.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample2 \
--batches "slices2,slices3,slices1,slices4" \
--landmark_domain 4 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample2/STAligner_CN.output &
```
### GraphST
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample2/sample_data \
--input_data Data/DLPFC_sample2/IntergrationRe/GraphST.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample2 \
--batches "slices2,slices3,slices1,slices4" \
--landmark_domain 1 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample2/STAligner_GraphST.output &
```
### MENDER
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample2/sample_data \
--input_data Data/DLPFC_sample2/IntergrationRe/MENDER.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample2 \
--batches "slices2,slices3,slices1,slices4" \
--landmark_domain 3 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample2/STAligner_MENDER.output &
```
### PRECAST
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample2/sample_data \
--input_data Data/DLPFC_sample2/IntergrationRe/PRECAST.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample2 \
--batches "slices2,slices3,slices1,slices4" \
--landmark_domain 2 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample2/STAligner_PRECAST.output &
```
### Spado
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample2/sample_data \
--input_data Data/DLPFC_sample2/IntergrationRe/Spado.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample2 \
--batches "slices2,slices3,slices1,slices4" \
--landmark_domain "Domain_5_5" --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample2/STAligner_Spado.output &
```
### SPIRAL
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample2/sample_data \
--input_data Data/DLPFC_sample2/IntergrationRe/SPIRAL.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample2 \
--batches "slices2,slices3,slices1,slices4" \
--landmark_domain 4 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample2/STAligner_SPIRAL.output &
```
### STAIG
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample2/sample_data \
--input_data Data/DLPFC_sample2/IntergrationRe/STAIG.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample2 \
--batches "slices2,slices3,slices1,slices4" \
--landmark_domain 1 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample2/STAligner_STAIG.output &
```
### STAligner
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample2/sample_data \
--input_data Data/DLPFC_sample2/IntergrationRe/STAligner.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample2 \
--batches "slices2,slices3,slices1,slices4" \
--landmark_domain 5 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample2/STAligner_STAligner.output &
```
## DLPFC S3
### Banksy
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample3/sample_data \
--input_data Data/DLPFC_sample3/IntergrationRe/Banksy.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample3 \
--batches "slices2,slices4,slices1,slices3" \
--landmark_domain 2 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample3/STAligner_Banksy.output &
```
### CellCharter
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample3/sample_data \
--input_data Data/DLPFC_sample3/IntergrationRe/CellCharter.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample3 \
--batches "slices2,slices4,slices1,slices3" \
--landmark_domain 0 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample3/STAligner_CellCharter.output &
```
### CN
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample3/sample_data \
--input_data Data/DLPFC_sample3/IntergrationRe/CN.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample3 \
--batches "slices2,slices4,slices1,slices3" \
--landmark_domain 6 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample3/STAligner_CN.output &
```
### GraphST
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample3/sample_data \
--input_data Data/DLPFC_sample3/IntergrationRe/GraphST.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample3 \
--batches "slices2,slices4,slices1,slices3" \
--landmark_domain 3 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample3/STAligner_GraphST.output &
```
### MENDER
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample3/sample_data \
--input_data Data/DLPFC_sample3/IntergrationRe/MENDER.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample3 \
--batches "slices2,slices4,slices1,slices3" \
--landmark_domain 4 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample3/STAligner_MENDER.output &
```
### PRECAST
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample3/sample_data \
--input_data Data/DLPFC_sample3/IntergrationRe/PRECAST.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample3 \
--batches "slices2,slices4,slices1,slices3" \
--landmark_domain 7 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample3/STAligner_PRECAST.output &
```
### Spado
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample3/sample_data \
--input_data Data/DLPFC_sample3/IntergrationRe/Spado.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample3 \
--batches "slices2,slices4,slices1,slices3" \
--landmark_domain "Domain_7_6" --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample3/STAligner_Spado.output &
```
### SPIRAL
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample3/sample_data \
--input_data Data/DLPFC_sample3/IntergrationRe/SPIRAL.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample3 \
--batches "slices2,slices4,slices1,slices3" \
--landmark_domain 3 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample3/STAligner_SPIRAL.output &
```
### STAIG
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample3/sample_data \
--input_data Data/DLPFC_sample3/IntergrationRe/STAIG.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample3 \
--batches "slices2,slices4,slices1,slices3" \
--landmark_domain 3 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample3/STAligner_STAIG.output &
```
### STAligner
```python
nohup python Run_STAligner.py \
--input_file Data/DLPFC_sample3/sample_data \
--input_data Data/DLPFC_sample3/IntergrationRe/STAligner.h5ad \
--output_file Benchmark/Alignment/Result/DLPFC_sample3 \
--batches "slices2,slices4,slices1,slices3" \
--landmark_domain 2 --landmark_domain_original WM --domain predicted_domain \
--step 10 --runNormalization True --hvgs 10000 \
>Benchmark/Alignment/Result/DLPFC_sample3/STAligner_STAligner.output &
```
## MERFISH Preoptic
### Banksy
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH/sample_data \
--input_data Data/MERFISH/IntergrationRe/Banksy.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH \
--batches "slices1,slices4,slices2,slices5,slices3" \
--landmark_domain 2 --landmark_domain_original BST --domain predicted_domain \
--step 10 --runNormalization False  \
>Benchmark/Alignment/Result/MERFISH/STAligner_Banksy.output &
```
### CellCharter
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH/sample_data \
--input_data Data/MERFISH/IntergrationRe/CellCharter.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH \
--batches "slices1,slices4,slices2,slices5,slices3" \
--landmark_domain 1 --landmark_domain_original BST --domain predicted_domain \
--step 10 --runNormalization False  \
>Benchmark/Alignment/Result/MERFISH/STAligner_CellCharter.output &
```
### CN
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH/sample_data \
--input_data Data/MERFISH/IntergrationRe/CN.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH \
--batches "slices1,slices4,slices2,slices5,slices3" \
--landmark_domain 1 --landmark_domain_original BST --domain predicted_domain \
--step 10 --runNormalization False  \
>Benchmark/Alignment/Result/MERFISH/STAligner_CN.output &
```
### GraphST
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH/sample_data \
--input_data Data/MERFISH/IntergrationRe/GraphST.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH \
--batches "slices1,slices4,slices2,slices5,slices3" \
--landmark_domain 4 --landmark_domain_original BST --domain predicted_domain \
--step 10 --runNormalization False  \
>Benchmark/Alignment/Result/MERFISH/STAligner_GraphST.output &
```
### MENDER
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH/sample_data \
--input_data Data/MERFISH/IntergrationRe/MENDER.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH \
--batches "slices1,slices4,slices2,slices5,slices3" \
--landmark_domain 1 --landmark_domain_original BST --domain predicted_domain \
--step 10 --runNormalization False  \
>Benchmark/Alignment/Result/MERFISH/STAligner_MENDER.output &
```
### NicheCompass
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH/sample_data \
--input_data Data/MERFISH/IntergrationRe/NicheCompass.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH \
--batches "slices1,slices4,slices2,slices5,slices3" \
--landmark_domain 0 --landmark_domain_original BST --domain predicted_domain \
--step 10 --runNormalization False  \
>Benchmark/Alignment/Result/MERFISH/STAligner_NicheCompass.output &
```
### PRECAST
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH/sample_data \
--input_data Data/MERFISH/IntergrationRe/PRECAST.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH \
--batches "slices1,slices4,slices2,slices5,slices3" \
--landmark_domain 8 --landmark_domain_original BST --domain predicted_domain \
--step 10 --runNormalization False  \
>Benchmark/Alignment/Result/MERFISH/STAligner_PRECAST.output &
```
### Spado
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH/sample_data \
--input_data Data/MERFISH/IntergrationRe/Spado.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH \
--batches "slices1,slices4,slices2,slices5,slices3" \
--landmark_domain "Domain_8_2" --landmark_domain_original BST --domain predicted_domain \
--step 10 --runNormalization False  \
>Benchmark/Alignment/Result/MERFISH/STAligner_Spado.output &
```
### SPIRAL
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH/sample_data \
--input_data Data/MERFISH/IntergrationRe/SPIRAL.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH \
--batches "slices1,slices4,slices2,slices5,slices3" \
--landmark_domain 1 --landmark_domain_original BST --domain predicted_domain \
--step 10 --runNormalization False  \
>Benchmark/Alignment/Result/MERFISH/STAligner_SPIRAL.output &
```
### STAIG
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH/sample_data \
--input_data Data/MERFISH/IntergrationRe/STAIG.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH \
--batches "slices1,slices4,slices2,slices5,slices3" \
--landmark_domain 6 --landmark_domain_original BST --domain predicted_domain \
--step 10 --runNormalization False  \
>Benchmark/Alignment/Result/MERFISH/STAligner_STAIG.output &
```
### STAligner
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH/sample_data \
--input_data Data/MERFISH/IntergrationRe/STAligner.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH \
--batches "slices1,slices4,slices2,slices5,slices3" \
--landmark_domain 4 --landmark_domain_original BST --domain predicted_domain \
--step 10 --runNormalization False  \
>Benchmark/Alignment/Result/MERFISH/STAligner_STAligner.output &
```
## MERFISH Brain S3
### Banksy
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH_Brain_S3/sample_data \
--input_data Data/MERFISH_Brain_S3/IntergrationRe/Banksy.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH_Brain_S3 \
--batches "slices2,slices1" \
--landmark_domain 3 --landmark_domain_original "corpus callosum" --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/MERFISH_Brain_S3/STAligner_Banksy.output &
```
### CellCharter
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH_Brain_S3/sample_data \
--input_data Data/MERFISH_Brain_S3/IntergrationRe/CellCharter.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH_Brain_S3 \
--batches "slices2,slices1" \
--landmark_domain 2 --landmark_domain_original "corpus callosum" --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/MERFISH_Brain_S3/STAligner_CellCharter.output &
```
### CN
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH_Brain_S3/sample_data \
--input_data Data/MERFISH_Brain_S3/IntergrationRe/CN.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH_Brain_S3 \
--batches "slices2,slices1" \
--landmark_domain 4 --landmark_domain_original "corpus callosum" --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/MERFISH_Brain_S3/STAligner_CN.output &
```
### GraphST
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH_Brain_S3/sample_data \
--input_data Data/MERFISH_Brain_S3/IntergrationRe/GraphST.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH_Brain_S3 \
--batches "slices2,slices1" \
--landmark_domain 6 --landmark_domain_original "corpus callosum" --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/MERFISH_Brain_S3/STAligner_GraphST.output &
```
### MENDER
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH_Brain_S3/sample_data \
--input_data Data/MERFISH_Brain_S3/IntergrationRe/MENDER.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH_Brain_S3 \
--batches "slices2,slices1" \
--landmark_domain 3 --landmark_domain_original "corpus callosum" --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/MERFISH_Brain_S3/STAligner_MENDER.output &
```
### NicheCompass
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH_Brain_S3/sample_data \
--input_data Data/MERFISH_Brain_S3/IntergrationRe/NicheCompass.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH_Brain_S3 \
--batches "slices2,slices1" \
--landmark_domain 1 --landmark_domain_original "corpus callosum" --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/MERFISH_Brain_S3/STAligner_NicheCompass.output &
```
### PRECAST
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH_Brain_S3/sample_data \
--input_data Data/MERFISH_Brain_S3/IntergrationRe/PRECAST.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH_Brain_S3 \
--batches "slices2,slices1" \
--landmark_domain 5 --landmark_domain_original "corpus callosum" --domain predicted_domain \
--step 5 --runNormalization False \
>Benchmark/Alignment/Result/MERFISH_Brain_S3/STAligner_PRECAST.output &
```
### Spado
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH_Brain_S3/sample_data \
--input_data Data/MERFISH_Brain_S3/IntergrationRe/Spado.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH_Brain_S3 \
--batches "slices2,slices1" \
--landmark_domain "Domain_8_4" --landmark_domain_original "corpus callosum" --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/MERFISH_Brain_S3/STAligner_Spado.output &
```
### SPIRAL
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH_Brain_S3/sample_data \
--input_data Data/MERFISH_Brain_S3/IntergrationRe/SPIRAL.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH_Brain_S3 \
--batches "slices2,slices1" \
--landmark_domain 0 --landmark_domain_original "corpus callosum" --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/MERFISH_Brain_S3/STAligner_SPIRAL.output &
```
### STAIG
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH_Brain_S3/sample_data \
--input_data Data/MERFISH_Brain_S3/IntergrationRe/STAIG.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH_Brain_S3 \
--batches "slices2,slices1" \
--landmark_domain 4 --landmark_domain_original "corpus callosum" --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/MERFISH_Brain_S3/STAligner_STAIG.output &
```
### STAligner
```python
nohup python Run_STAligner.py \
--input_file Data/MERFISH_Brain_S3/sample_data \
--input_data Data/MERFISH_Brain_S3/IntergrationRe/STAligner.h5ad \
--output_file Benchmark/Alignment/Result/MERFISH_Brain_S3 \
--batches "slices2,slices1" \
--landmark_domain 2 --landmark_domain_original "corpus callosum" --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/MERFISH_Brain_S3/STAligner_STAligner.output &
```
## STARMap
### Banksy
```python
nohup python Run_STAligner.py \
--input_file Data/STARMap/sample_data \
--input_data Data/STARMap/IntergrationRe/Banksy.h5ad \
--output_file Benchmark/Alignment/Result/STARMap \
--batches "slices1,slices3,slices2" \
--landmark_domain 2 --landmark_domain_original "4" --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/STARMap/STAligner_Banksy.output &
```
### CellCharter
```python
nohup python Run_STAligner.py \
--input_file Data/STARMap/sample_data \
--input_data Data/STARMap/IntergrationRe/CellCharter.h5ad \
--output_file Benchmark/Alignment/Result/STARMap \
--batches "slices1,slices3,slices2" \
--landmark_domain 0 --landmark_domain_original "4" --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/STARMap/STAligner_CellCharter.output &
```
### CN
```python
nohup python Run_STAligner.py \
--input_file Data/STARMap/sample_data \
--input_data Data/STARMap/IntergrationRe/CN.h5ad \
--output_file Benchmark/Alignment/Result/STARMap \
--batches "slices1,slices3,slices2" \
--landmark_domain 1 --landmark_domain_original "4" --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/STARMap/STAligner_CN.output &
```
### GraphST
```python
nohup python Run_STAligner.py \
--input_file Data/STARMap/sample_data \
--input_data Data/STARMap/IntergrationRe/GraphST.h5ad \
--output_file Benchmark/Alignment/Result/STARMap \
--batches "slices1,slices3,slices2" \
--landmark_domain 1 --landmark_domain_original "4" --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/STARMap/STAligner_GraphST.output &
```
### MENDER
```python
nohup python Run_STAligner.py \
--input_file Data/STARMap/sample_data \
--input_data Data/STARMap/IntergrationRe/MENDER.h5ad \
--output_file Benchmark/Alignment/Result/STARMap \
--batches "slices1,slices3,slices2" \
--landmark_domain 0 --landmark_domain_original "4" --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/STARMap/STAligner_MENDER.output &
```
### NicheCompass
```python
nohup python Run_STAligner.py \
--input_file Data/STARMap/sample_data \
--input_data Data/STARMap/IntergrationRe/NicheCompass.h5ad \
--output_file Benchmark/Alignment/Result/STARMap \
--batches "slices1,slices3,slices2" \
--landmark_domain 1 --landmark_domain_original "4" --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/STARMap/STAligner_NicheCompass.output &
```
### PRECAST
```python
nohup python Run_STAligner.py \
--input_file Data/STARMap/sample_data \
--input_data Data/STARMap/IntergrationRe/PRECAST.h5ad \
--output_file Benchmark/Alignment/Result/STARMap \
--batches "slices1,slices3,slices2" \
--landmark_domain 2 --landmark_domain_original "4" --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/STARMap/STAligner_PRECAST.output &
```
### Spado
```python
nohup python Run_STAligner.py \
--input_file Data/STARMap/sample_data \
--input_data Data/STARMap/IntergrationRe/Spado.h5ad \
--output_file Benchmark/Alignment/Result/STARMap \
--batches "slices1,slices3,slices2" \
--landmark_domain "Domain_4_4" --landmark_domain_original "4" --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/STARMap/STAligner_Spado.output &
```
### SPIRAL
```python
nohup python Run_STAligner.py \
--input_file Data/STARMap/sample_data \
--input_data Data/STARMap/IntergrationRe/SPIRAL.h5ad \
--output_file Benchmark/Alignment/Result/STARMap \
--batches "slices1,slices3,slices2" \
--landmark_domain 1 --landmark_domain_original "4" --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/STARMap/STAligner_SPIRAL.output &
```
### STAIG
```python
nohup python Run_STAligner.py \
--input_file Data/STARMap/sample_data \
--input_data Data/STARMap/IntergrationRe/STAIG.h5ad \
--output_file Benchmark/Alignment/Result/STARMap \
--batches "slices1,slices3,slices2" \
--landmark_domain 1 --landmark_domain_original "4" --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/STARMap/STAligner_STAIG.output &
```
### STAligner
```python
nohup python Run_STAligner.py \
--input_file Data/STARMap/sample_data \
--input_data Data/STARMap/IntergrationRe/STAligner.h5ad \
--output_file Benchmark/Alignment/Result/STARMap \
--batches "slices1,slices3,slices2" \
--landmark_domain 3 --landmark_domain_original "4" --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/STARMap/STAligner_STAligner.output &
```















