# Load necessary libraries for spatial transcriptomics analysis
library(SpaDo)
library(Seurat)
library(magrittr)
library(SeuratDisk)
library(anndata)
library(dplyr)
library(tidyr)
library(fastcluster)
library(pryr)  
library(SeuratData)
library(reticulate)

# Specify the Python environment to use (required by SeuratDisk for Python-R integration)
# use_python("/home/dongkj/anaconda3/envs/MultiSpatial/bin/python", required=TRUE)

# Function to run the spatial analysis pipeline
# This function processes input spatial transcriptomics data, performs data normalization, clustering, and dimensionality reduction.
# It also computes cell type distribution and applies hierarchical clustering to define spatial domains.
run_spatial_analysis <- function(input_file, output_file, sample, hvgs, n_pcs , nclust, runNormalization ) {
  # Parameters:
  # input_file (abbreviation: i): Path to the multi-slice data that needs to be integrated. All slices should be contained in a single h5ad file. Refer to the example data for the specific format. This is a required parameter.
  # output_file (abbreviation: o): Path to the directory where the output results will be saved. This is a required parameter.
  # sample (abbreviation: s): Name for the output result files. It is recommended to use the model name, for example, setting it as "SpaDo" will result in output files named "SpaDo.h5ad". This is a required parameter.
  # nclust (abbreviation: c): Number of identified domains. This is a required parameter.
  # hvgs (abbreviation: v): The number of highly variable genes identified if data normalization is performed. This is an optional parameter, with the default set to 2000.
  # runNormalization (abbreviation: n): A boolean variable indicating whether data normalization should be performed before integration. This is an optional parameter, with the default set to TRUE.
  # n_pcs: The number of principal components set when performing PCA dimensionality reduction on the data. This is an optional parameter, with the default set to 50.
  
  # Read the input .h5ad file and extract the expression matrix and metadata
  data <- read_h5ad(input_file) 
  exp <- data$X  
  metadata <- data$obs
  rownames(exp) <- metadata$barcode 
  
  # Create a Seurat object from the expression matrix and metadata
  data <- CreateSeuratObject(counts = t(exp), min.cells = 0, min.features = 0, meta.data = metadata)
  data$slices <- as.character(data$slices)
  
  # If normalization is requested, perform normalization, PCA, and clustering
  if (runNormalization) {
    data <- FindVariableFeatures(data, nfeatures = hvgs)  
    data <- NormalizeData(data) 
    data <- ScaleData(data)
    data <- RunPCA(data, npcs = n_pcs)  
    data <- FindNeighbors(data, dims = 1:n_pcs) 
    data <- FindClusters(data, resolution = 2, random.seed = 666) 
  } else {
    # If normalization is not requested, still normalize and scale, but use the raw count matrix for the layers
    data <- FindVariableFeatures(data, nfeatures = ncol(exp))  
    data <- NormalizeData(data) 
    data <- ScaleData(data) 
    data@assays[["RNA"]]@layers[["scale.data"]] <- t(exp) 
    data <- RunPCA(data, npcs = n_pcs) 
    data <- FindNeighbors(data, dims = 1:n_pcs) 
    data <- FindClusters(data, resolution = 2, random.seed = 666) 
  }
  
  # Initialize lists to store coordinates and cell type information for each slice
  coor_list <- list()
  cellAn_list <- list()
  slices <- unique(data$slices) 
  
  # Loop over each slice to extract coordinates and cell type information
  for (s in slices) {
    slices_data_seurat <- subset(data, slices == s)  
    coor_data <- cbind(slices_data_seurat$X, slices_data_seurat$Y) 
    rownames(coor_data) <- slices_data_seurat$barcode  
    colnames(coor_data) <- c("X", "Y") 
    
    # If 'celltype' is available in the metadata, use it; otherwise, use the cluster information
    if ("celltype" %in% colnames(slices_data_seurat@meta.data)) {
      celltype <- as.character(slices_data_seurat$celltype) 
      names(celltype) <- slices_data_seurat$barcode 
    } else {
      celltype <- as.character(slices_data_seurat$seurat_clusters)  
      names(celltype) <- slices_data_seurat$barcode 
    }
    
    # Store the coordinates and cell types in respective lists
    coor_list[[s]] <- as.data.frame(coor_data)
    cellAn_list[[s]] <- celltype
  }
  
  # Record the start time and memory usage for performance tracking
  start_time <- Sys.time()
  start_mem <- pryr::mem_used()
  
  # Perform spatial cell type distribution analysis across multiple slices
  cell_type_distribution_multiple <- SpatialCellTypeDistribution_multiple(
    sample_information_coordinate_list = coor_list,
    sequence_resolution = 'single_cell',
    sample_information_cellType_list = cellAn_list
  )
  
  # Compute the distribution distance using Jensen-Shannon Divergence (JSD) for spatial domain identification
  # "no_cores" can be set to a higher value for faster operation. 
  # However, an error message is displayed if no_cores is set to a value greater than 1 when the terminal is running.
  cell_type_distribution = cell_type_distribution_multiple$cell_type_distribution_combine
  DD <- DistributionDistance(
    cell_type_distribution,
    distance = "JSD", no_cores = 1
  )
  
  # Perform hierarchical clustering to identify spatial domains (clusters)
  domain_hclust <- DomainHclust(distribution_distance = DD, autoselection = FALSE, domain_num = nclust)
  
  # Record the end time and memory usage
  end_time <- Sys.time()
  end_mem <- pryr::mem_used()
  
  # Calculate the time taken and memory used during the analysis
  time_taken <- end_time - start_time
  mem_used <- end_mem - start_mem
  
  # Extract and format the clustering results
  re <- domain_hclust[["hclust_result_df"]]
  rownames(re) <- sapply(rownames(re), function(x) {
    unlist(strsplit(x, split = "_"))[2]
  })
  
  # Match the cluster results to the original Seurat object barcodes
  index <- match(as.character(data$barcode), rownames(re))
  re <- re[index, ]
  predicted_domain <- re[, ncol(re)]
  
  # Update the Seurat object with the predicted domain and other results
  data_seurat <- data
  data_seurat$predicted_domain <- predicted_domain
  data_seurat$original_domain <- as.character(data_seurat$original_domain)
  data_seurat$slices <- as.character(data_seurat$slices)
  
  # Rename the RNA assay and set up the embedding for visualization
  data_seurat[["RNA3"]] <- as(object = data_seurat[["RNA"]], Class = "Assay")
  DefaultAssay(data_seurat) <- "RNA3"
  data_seurat[["RNA"]] <- NULL
  data_seurat <- RenameAssays(object = data_seurat, RNA3 = 'RNA')
  
  # Create the embedding object and add it to the Seurat object
  embedding_exp <- cell_type_distribution
  rowname_embedding_exp <- sapply(rownames(embedding_exp), function(x) {
    unlist(strsplit(x, split = "_"))[2]
  })
  rownames(embedding_exp) <- rowname_embedding_exp
  colnames(embedding_exp) <- paste0("embedding_", 1:ncol(embedding_exp))
  embedding <- CreateDimReducObject(embeddings = embedding_exp, key = "embedding_", assay = "RNA")
  data_seurat[["embedding"]] <- embedding
  
  # Save the results to an H5Seurat file and convert it to an H5AD file
  h5Seurat_file_1 <- paste(output_file, paste(sample, "h5Seurat", sep = "."), sep = "/")
  SaveH5Seurat(data_seurat, filename = h5Seurat_file_1, overwrite = TRUE)
  Convert(h5Seurat_file_1, dest = "h5ad", overwrite = TRUE)
  
  # Print the execution time and memory usage
  cat("Time taken: ", time_taken, "\n")
  cat("Memory used: ", mem_used, "\n")
  
  # Log the execution statistics to a CSV file
  log_df <- data.frame(
    Memory_Usage_MiB = as.numeric(mem_used) / (1024 ^ 2), 
    Execution_Time_s = as.numeric(time_taken, units = "secs")
  )
  write.table(log_df, paste(output_file, paste(sample, "_stats.csv", sep = ""), sep = "/"), col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
}

# Command line argument parsing using getopt for script execution
library(getopt)
spec <- matrix(
  c(
    'input_file',     'i', 1, 'character',  'Path to input file',
    'output_file',    'o', 1, 'character',  'Path to output directory',
    'sample',         's', 1, 'character',  'Sample name',
    'hvgs',           'v', 2, 'integer',    'Number of highly variable genes (default: 2000)',
    'n_pcs',          'p', 2, 'integer',    'Number of principal components (default: 50)',
    'nclust',         'c', 1, 'integer',    'Number of clusters',
    'runNormalization', 'n', 2, 'logical',  'Flag to run normalization (default: TRUE)'
  ),
  byrow = TRUE,
  ncol = 5
)

opt <- getopt(spec)

# Set default values for optional parameters if not provided by user
hvgs <- ifelse(is.null(opt$hvgs), 2000, opt$hvgs)
n_pcs <- ifelse(is.null(opt$n_pcs), 50, opt$n_pcs)
runNormalization <- ifelse(is.null(opt$runNormalization), TRUE, opt$runNormalization)

# Run the spatial analysis function with the parsed arguments
run_spatial_analysis(opt$input_file, opt$output_file, opt$sample, hvgs, n_pcs, opt$nclust, runNormalization)
