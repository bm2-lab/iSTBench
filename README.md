# iSTBench
We developed a comprehensive benchmarking pipeline to evaluate state-of-the-art multi-slice integration methods across diverse technologies. Our evaluation framework includes both multi-slice integration performance and three critical downstream applications that utilize the integrated embeddings, including spatial clustering, spatial alignment, and slice representation. Beyond benchmarking, we examine the relationship between upstream integration and downstream performance, identify task-specific limitations, and propose modular solutions to overcome them. To support this, we curated 19 spatial transcriptomics datasets from seven sources, encompassing multiple technologies such as 10X Visium, BaristaSeq, MERFISH, and STARMap, to benchmark eight multi-slice integration methods.

**8** methods are included:
- **BANKSY**:《BANKSY unifies cell typing and tissue domain segmentation for scalable spatial omics data analysis》
- **CellCharter**:《CellCharter reveals spatial cell niches associated with tissue remodeling and cell plasticity》
- **CN**:《Coordinated Cellular Neighborhoods Orchestrate Antitumoral Immunity at the Colorectal Cancer Invasive Front》
- **GraphST, GraphST-PASTE**:《Spatially informed clustering, integration, and deconvolution of spatial transcriptomics with GraphST》
- **MENDER**:《MENDER: fast and scalable tissue structure identification in spatial omics data》
- **NicheCompass**:《Large-scale characterization of cell niches in spatial atlases using bio-inspired graph learning》
- **SpaDo**:《Multi‐slice spatial transcriptome domain analysis with SpaDo》

# Benchmark framework
To rigorously assess the performance of multi-slice integration methods, we propose a comprehensive evaluation framework that examines five key areas: multi-slice integration, spatial clustering, spatial alignment, slice representation as well as the scalability of methods. 
In order to replicate the results obtained or to evaluate the methods using one's own data, it is possible to create a Python environment by entering the following command:
In order to replicate the results obtained or to evaluate the methods using your own data, it is possible to download the relevant code and sample data, and to create a Python environment by entering the following command:
```python
# download iSTBench
git clone https://github.com/bm2-lab/iSTBench.git

# set dir to folder
cd SCMMI_benchmark

# create the conda environment
conda env create -f environment.yml
```

It is imperative to note that both Banksy and DpaDo are constructed on the R language. Consequently, in order to execute the code and download the relevant R packages and dependency packages, it is strongly recommended to utilise R 4.3.2.
## 1. Multi-slice integration and Spatial clustering
In this section, we use multi-slice data as input, applying different methods to integrate the data and obtain the corresponding embeddings. Based on these integrated embeddings, we perform clustering to identify spatial domains. Different methods have varying requirements for the input data format. The data can be found in the directories "Data/sample_all_data" and "data/sample_data". The "sample_all_data" folder stores the merged multi-slice data, while the "sample_data" folder contains individual slice data. The embeddings and predicted domain information resulting from the integration are stored in the metadata of the corresponding files. The specific format can be referenced in the result files located in "Data/IntegrationRe."
Taking GraphST and MENDER as an example, the relevant code is as follows:
Command:
```python
# GraphST
nohup python ./RunModel/GraphST.py \
--input_file ./Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
--output_file ./Data/BaristaSeq/IntergrationRe \
--sample GraphST --nclust 7 --device cuda \
> ./Data/BaristaSeq/GraphST.output  &

# MENDER
nohup python MENDER.py \
--input_file ./Data/BaristaSeq/sample_all_data/Slices_combind_data.h5ad \
--output_file ./Data/BaristaSeq/IntergrationRe \
--sample MENDER --nclust 7 --tech BaristaSeq \
> ./Data/BaristaSeq/MENDER.output &
```
The complete code for running other methods is stored in the "./RunModel/TerminalRun.py" file. Different methods require distinct parameter settings. The specific meanings of the parameters can be referenced in the "./RunModel/parameters.txt" file.
## 2. Spatial alignment
In this section,  we evaluate the performance of different methods in spatial alignment. We used the STAligner framework to correct  coordinates between slices based on the integration embeddings and domains identification from each method. To demonstrate the effectiveness of spatial alignment, we applied rotations to slices within the same data. This approach allowed us to clearly illustrate the impact of spatial alignment. Specifically, for each slice in a data, rotations of 0°, 20°, 40°, and so on were applied. The code for performing slices rotation on the specified data is as follows:
```python
nohup python ./Alignment/Rotate_spatial.py \
--input_file ./Data/BaristaSeq/sample_data \
--output_file ./Data/BaristaSeq/sample_data_rotate \
--batches "slices1,slices2,slices3" \
--angles 40 20 0 \
> ./Alignment/rotate.output &
```
The specific meanings of the parameters can be referenced in the "./Alignment/Rotate_parameters.txt" file. 

As an example, the relevant code for spatial alignment based on the embeddings and domains generated by MENDER is as follows:
```python
nohup python ./Alignment/STAligner.py \
--input_file ./Data/BaristaSeq/sample_data_rotate \
--input_data ./Data/BaristaSeq/IntergrationRe/MENDER.h5ad \
--output_file ./Alignment/Result/BaristaSeq \
--batches "slices1,slices2,slices3" \
--landmark_domain 4 --landmark_domain_original VISp_wm --domain predicted_domain \
--step 10 --runNormalization False \
>./Alignment/Result/BaristaSeq/MENDER.output &
```
The specific meanings of the parameters can be referenced in the "./Alignment/STAligner_parameters.txt" file. 
