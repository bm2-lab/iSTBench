## 2.1Alignment_metric.R
The script reads the results of the spatial alignment of different methods for the specified dataset and calculates the relevant indicators for the spatial alignment evaluation. You need to run the script in terminal submission.
```python
Parameters:
input_file (abbreviation: i): Path to directory containing alignemnt result file of each method.
slices (abbreviation: s): The file name of each slice. Different slices are separated by "," and need to be arranged in adjacent order. For example: "slices1,slices2,slices3". Same as the "batches" parameter when running "Rotate_spatial.py".
Terminal run:
```
```python
nohup Rscript Analysis/Spatial_alignment/code/2.1Alignment_metric.R \
-i "/NFS2_home/NFS2_home_3/dongkj/home_dkj/FD_yzy/Result/GitHub_test/iSTBench/Benchmark/Alignment/Result/BaristaSeq" \
-s "slices1,slices2,slices3" \
> /NFS2_home/NFS2_home_3/dongkj/home_dkj/FD_yzy/Result/GitHub_test/iSTBench/Benchmark/Alignment/Result/BaristaSeq/Metric/metric.output &
```
## 2.2Analysis_each_data.R
The script generates spatial alignment result plots for each dataset based on different methods. It includes both 2D and 3D illustrations, as well as result plots for the corresponding metrics. The example code demonstrates how to plot results for the BaristaSeq dataset. The script needs to be run from the command line.
## 2.3Analyze_all_alignment.R
The script reads the alignment results for each dataset and generates an overall result plot, corresponding to Extended Fig. 4. The alignment results for different datasets are stored in "Analysis/Spatial_alignment/result/integration_re.csv" and "Analysis/Spatial_alignment/result/integration_sd.csv". The script needs to be run from the command line.


