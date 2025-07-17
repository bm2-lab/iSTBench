# iSTBench
We developed a comprehensive benchmarking pipeline to evaluate state-of-the-art multi-slice integration methods across diverse technologies. Our evaluation framework includes both multi-slice integration performance and three critical downstream applications that utilize the integrated embeddings, including spatial clustering, spatial alignment, and slice representation. To support this, we curated 19 spatial transcriptomics datasets from seven sources, encompassing multiple technologies such as 10X Visium, BaristaSeq, MERFISH, and STARMap, to benchmark 12 multi-slice integration methods. For each task, we perform detailed analyses of the methods and provide actionable recommendations. Our results reveal substantial data-dependent variation in performance across tasks. We further investigate the relationships between upstream and downstream tasks, showing that downstream performance often depends on upstream quality.

**12** methods are included:
| Method                                                                         | Article                                                                     | Title                                                                  | Time |
|--------------------------------------------------------------------------------|-----------------------------------------------------------------------------|------------------------------------------------------------------------|------|
| [BANKSY](https://github.com/prabhakarlab/Banksy)                               | [Nature Genetics](https://www.nature.com/articles/s41588-024-01664-3)       |BANKSY unifies cell typing and tissue domain segmentation for scalable spatial omics data analysis          | 2024 |
| [CellCharter](https://github.com/CSOgroup/cellcharter)                         | [Nature Genetics](https://www.nature.com/articles/s41588-023-01588-4)       |CellCharter reveals spatial cell niches associated with tissue remodeling and cell plasticity               | 2023 |
| [CN](https://github.com/nolanlab/NeighborhoodCoordination)                     | [Cell](https://www.cell.com/cell/fulltext/S0092-8674(20)31385-4)            |Coordinated Cellular Neighborhoods Orchestrate Antitumoral Immunity at the Colorectal Cancer Invasive Front | 2020 |
| [GraphST, GraphST-PASTE](https://github.com/JinmiaoChenLab/GraphST)            | [Nature Communications](https://www.nature.com/articles/s41467-023-36796-3) |Spatially informed clustering, integration, and deconvolution of spatial transcriptomics with GraphST       | 2023 |
| [MENDER](https://github.com/yuanzhiyuan/MENDER)                                | [Nature Communications](https://www.nature.com/articles/s41467-023-44367-9) |MENDER: fast and scalable tissue structure identification in spatial omics data                             | 2024 |
| [NicheCompass](https://github.com/Lotfollahi-lab/nichecompass)                 | [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.02.21.581428v1)      |Large-scale characterization of cell niches in spatial atlases using bio-inspired graph learning            | 2024 |
| [PRECAST](https://github.com/feiyoung/PRECAST)                                 | [Nature Communications](https://www.nature.com/articles/s41467-023-35947-w) |Probabilistic embedding, clustering, and alignment for integrating spatial transcriptomics data with PRECAST| 2023 |
| [SpaDo](https://github.com/bm2-lab/SpaDo)                                      | [Geonome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03213-x#:~:text=To%20this%20end%2C%20we%20propose%20SpaDo%20%28multi-slice%20spatial,transcriptome%20analysis%20at%20both%20single-cell%20and%20spot%20resolution.)      |Multi‐slice spatial transcriptome domain analysis with SpaDo           | 2024 |
| [SPIRAL](https://github.com/guott15/SPIRAL)                                    | [Geonome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03078-6) |SPIRAL: integrating and aligning spatially resolved transcriptomics data across different experiments, conditions, and technologies| 2023 |
| [STAIG](https://github.com/y-itao/STAIG)                                       | [Nature Communications](https://www.nature.com/articles/s41467-025-56276-0) |STAIG: Spatial transcriptomics analysis via image-aided graph contrastive learning for domain exploration and alignment-free integration| 2025 |
| [STAligner](https://github.com/zhoux85/STAligner)                              | [Nature Computational Science](https://www.nature.com/articles/s43588-023-00543-x) |STAligner enables the integration and alignment of multiple spatial transcriptomics datasets          | 2023 |


All benchmark data and corresponding results have been uploaded to Zenodo and can be accessed [here](https://zenodo.org/records/14906156). Please download all files into the ./Data directory.

# Benchmark framework
To rigorously assess the performance of multi-slice integration methods, we propose a comprehensive evaluation framework covering five key areas: multi-slice integration, spatial clustering, spatial alignment, slice representation, and method scalability.

To replicate the results or evaluate the methods with your own data, you can download the relevant code and sample data and set up a Python environment by entering the following command:
```python
# download iSTBench
git clone https://github.com/bm2-lab/iSTBench.git

# set dir to folder
cd iSTBench

# create the conda environment
conda env create -f environment.yaml
```
It is important to note that both Banksy and SpaDo are built using the R language. Therefore, to run the code and download the necessary R packages and dependencies, it is strongly recommended to use R version 4.3.2.
## 1. Multi-slice integration and spatial clustering
In this section, we use multi-slice data as input, applying different methods to integrate the data and generate the corresponding embeddings. Based on these integrated embeddings, we perform clustering to identify spatial domains. Each method has specific requirements for the input data format. The data can be found in the "Data/sample_all_data" and "data/sample_data" directories. The "sample_all_data" folder contains the merged multi-slice data, while the "sample_data" folder includes individual slice data. The embeddings and predicted domain information from the integration are stored in the metadata of the corresponding files. The specific format can be referenced in the result files located in "Data/IntegrationRe."

As an example, here is the relevant code for GraphST and MENDER on BaristaSeq dataset:
```python
# GraphST
nohup python Benchmark/RunModel/GraphST/Run_GraphST.py \
--input_file Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/BaristaSeq/IntergrationRe \
--sample GraphST --nclust 6 --device cuda \
> Data/BaristaSeq/IntergrationRe/GraphST.output  &

# MENDER
nohup python Benchmark/RunModel/Run_MENDER.py \
--input_file Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/BaristaSeq/IntergrationRe \
--sample MENDER --nclust 6 --tech BaristaSeq \
> Data/BaristaSeq/IntergrationRe/MENDER.output &
```
The "input_file" and "output_file" should be set to the exact paths of the input and output files, depending on the actual setup. The complete code for running other methods is available in the "Benchmark/RunModel/TerminalRun.md" file. Each method has specific parameter settings, and the details of these parameters can be found in the "Benchmark/RunModel/parameters.md" file.
## 2. Spatial alignment
In this section, we evaluate the performance of different methods in spatial alignment.  Current spatial alignment methods can be divided into two types: integration-based and non-integration-based methods. To comprehensively evaluate performance across these categories, we selected PASTE and STalign as representatives of non-integration-based methods, and STAligner and SPACEL as representatives of integration-based methods. For integration-based alignment, we applied the results of different multi-slice integration mthods as inputs to the STAligner and SPACEL pipelines. 

As an example, the following code demonstrates spatial alignment using PASTE and STalign, as well as SPACEL and STAligner based on the embeddings and domain labels generated by MENDER on BaristaSeq:
```python
# PASTE
nohup python Benchmark/Alignment/Run_PASTE.py \
--input_file Data/BaristaSeq/sample_data \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" --step 10 \
> Benchmark/Alignment/Result/BaristaSeq/PASTE.output &

# STalign
nohup python Benchmark/Alignment/Run_STalign.py \
--input_file Data/BaristaSeq/sample_data \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" --step 10 \
>  Benchmark/Alignment/Result/BaristaSeq/STalign.output &

# SPACEL based MENDER result
nohup python Benchmark/Alignment/Run_SPACEL.py \
--input_file Data/BaristaSeq/sample_data \
--input_data Data/BaristaSeq/IntergrationRe/MENDER.h5ad \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" --step 10 \
> Benchmark/Alignment/Result/BaristaSeq/SPACEL_MENDER.output &
# You just need to change imput_path to use the results of different methods

# STAligner based MENDER result
nohup python Benchmark/Alignment/Run_STAligner.py \
--input_file Data/BaristaSeq/sample_data \
--input_data Data/BaristaSeq/IntergrationRe/MENDER.h5ad \
--output_file Benchmark/Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" \
--landmark_domain 4 --landmark_domain_original VISp_wm --domain predicted_domain \
--step 10 --runNormalization False \
>Benchmark/Alignment/Result/BaristaSeq/STAlignerMENDER.output &
# You just need to change imput_path to use the results of different methods

```
The "input_file", "input_data" and "output_file" should be set to the exact paths of the input and output files, depending on the actual setup. The alignment results are stored in the file "Benchmark/Alignment/Result/BaristaSeq". The specific meanings of the parameters can be referenced in the "Benchmark/Alignment/parameters.md" file. The relevant code for spatial alignment based on other datasets can be found in "Benchmark/Alignment/TerminalRun.md".
## 3. Slice representation
In this section, we use the abundance of identified spatial domains in each slice as the representation. To do this, domain information must first be obtained using the integration method, and then slice representations are generated based on domain abundance. Taking MENDER as an example, the relevant code is as follows:
```python
# Firstly, domains are identified based on MENDER, where the number of domains is set to 6 or others
nohup python Benchmark/RunModel/Run_MENDER.py \
--input_file Data/TNBC/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/TNBC/SlicesEmbedding/MENDER/6 \
--sample MENDER --nclust 6 --runNormalization False --tech MIBI \
> Data/TNBC/SlicesEmbedding/MENDER/MENDER6.output &

nohup python Benchmark/RunModel/Run_MENDER.py \
--input_file Data/TNBC/sample_all_data/Slices_combind_data.h5ad \
--output_file Data/TNBC/SlicesEmbedding/MENDER/14 \
--sample MENDER --nclust 14 --runNormalization False --tech MIBI \
> Data/TNBC/SlicesEmbedding/MENDER/MENDER14.output &

# Slices are represented and clustered based on domain abundance.
nohup Rscript Benchmark/SliceRepresentation/SliceRepresentation.R \
-f Data/TNBC/SlicesEmbedding/MENDER \
-m MENDER -n 3 \
> Data/TNBC/SlicesEmbedding/MENDER/Metric/SlicesClustering.output &
```
The specific meanings of the parameters can be referenced in the "Benchmark/SliceRepresentation/parameters.md" file. 
# Modular solution
## Domain-Relationship-Aware Alignment Method (Dr.A)
The code for spatial alignment using DR.A is as follows:
```python
nohup python Improve/DR.A/DR.A.py \
--data_path Data/BaristaSeq/sample_data_rotate \
--integrated_path Data/BaristaSeq/IntergrationRe \
--output_path Improve/DR.A/Result/BaristaSeq/ \
--slices_file 'slices1,slices2,slices3' --slices '1,2,3' --model MENDER \
> Improve/DR.A/Result/MENDER.output &
```
The "data_path", "integrated_path" and "output_path" should be set to the exact paths of the input and output files, depending on the actual setup. The specific meanings of the parameters can be found in the "Improve/DR.A/parameters.md" file.
## Domain-Relationship-Aware Slice Representation Method (Dr.S)
The code for slice representation using DR.S is as follows:
```python
nohup Rscript Improve/DR.S/DR.S.R \
-f Data/TNBC/SlicesEmbedding/MENDER/14/MENDER.h5ad \
-o Data/TNBC/SlicesEmbedding/MENDER/14 \
-c 3 \
> Improve/DR.S/MENDER.output &
```
The "f" is the file path that stores the result of integration with the data, and in the obs of this data there needs to be a column "predicted_domain" to store domian information for each cell. “o” should be set to the exact paths of the output files, depending on the actual setup. "c" is the number of clusters to generate during the clustering step. The specific meanings of the parameters can be found in the "Improve/DR.S/parameters.md" file.

# Analysis
The relevant code for analyzing and visualizing the results is stored in the "Analysis" folder.
# Citation
# Contacts
bm2-lab@tongji.edu.cn

zhiyuan@fudan.edu.cn

2231451@tongji.edu.cn

