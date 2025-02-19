# Analysis the integration and domian identification results on BaristaSeq
## 1.1 Integration_metric.py
```python
nohup python Analysis/Integration_and_clustering/code/1.1Intergration_metric.py \
--input_file Data/BaristaSeq/IntergrationRe \
--input_path Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
--output_path Data/BaristaSeq/IntergrationRe/Metric/Intergration_value.csv \
> Data/BaristaSeq/IntergrationRe/Metric/Intergration_metric.output &
```
## 1.2 Domain_plotting.R
```python
nohup /opt/R/4.3.2/lib/R/bin/Rscript Analysis/Integration_and_clustering/code/1.2Domain_plotting.R \
Data/BaristaSeq/IntergrationRe/Metric/domain_plot.txt \
Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
Data/BaristaSeq/IntergrationRe/Metric \
> Data/BaristaSeq/IntergrationRe/Metric/domain_plot.output&
```
