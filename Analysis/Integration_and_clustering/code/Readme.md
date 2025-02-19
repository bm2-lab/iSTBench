## 1.1Intergration_metric.py
The script reads the results of the integration of different methods for the specified dataset and calculates the relevant indicators for the integration evaluation. You need to run the script in terminal submission
```python
Parameters:
input_file: Path to directory containing .h5ad files with integration results
input_path: Path to the original combined slices data in .h5ad format
output_path: Path where the resulting CSV file with integration metrics will be saved
seed: Random seed for reproducibility of results (default: 123)
```


## 1.2Domain_plotting.R
The script draws the results of domains identification based on the integrated results of different methods for the specified dataset, and computes the domians identification evaluation metrics: ARI and NMI. You need to run the script in terminal submission
Enter three parameters:

The first parameter: domain_plot.txt, is the path to store the results of the different methods

Second parameter: Path to directory containing.h5ad files with integration results

Third parameter: the path to save the results

**domain_plot.txt**: You need to store the integration results for each method and their corresponding paths in the "domain_plot.txt" file, in the format of "method:path". The "domain_plot.txt" file should be saved in the "IntegrationRe/Metric" folder of the respective dataset. You can refer to "Data/BaristaSeq/IntegrationRe/Metric/domain_plot.txt" as an example.


## 1.3Anylyze_each_data.R
The script further calculated the doamin identification evaluation metrics CHAOS and PAS, and drew the metric result graph of different methods domain identification and the domain visualization result graph. At the same time, the cell UMAP obtained from the integration results of different methods was drawn.
The script needs to be run from the command line.

## 1.4 Analyze_all_integration.R
This script performs an integrated analysis of the performance of all methods across different datasets. The script is divided into several parts:

1. **Read data**: This section reads the performance data of various models on multi-slice integration across all datasets. We provide sample results stored in "Analysis/Integration_and_clustering/result/integration_re_all.csv". You can also run the results for different methods on different datasets, read the corresponding results, and proceed with further analysis and visualization.

2. **Integration result plotting**: This part generates performance plots for integration for each dataset and visualizes the corresponding results as shown in Extended Fig. 2.

3. **Integration funcky heatmap**: This section plots heatmaps of integration performance results for different methods across all datasets, corresponding to the results in Extended Fig. 2. The plots shown in the text are layout-optimized using Illustrator.

4. **Clustering result plotting**: This part generates plots for domain identification performance, corresponding to Extended Fig. 3.

5. **Clustering funcky heatmap**: This section plots heatmaps of domain identification performance for different methods across all datasets, corresponding to the results shown in Extended Fig. 3. The plots in the text are layout-optimized using Illustrator.

6. **Time plot**: This part plots the results for runtime and memory usage for different methods.
