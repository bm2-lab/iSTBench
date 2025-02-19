# Load necessary libraries
library(Banksy)
library(harmony)
library(SummarizedExperiment)
library(SpatialExperiment)
library(Seurat)
library(SeuratDisk)
library(scran)
library(data.table)
library(scater)
library(cowplot)
library(ggplot2)
library(pryr)
library(getopt)

# Main function to run the Banksy pipeline
run_banksy <- function(input_file, output_file, sample, nclust, dataType, hvgs, runNormalization, k_geom, lambda, SEED) {
  # Parameters:
  # input_file (abbreviation: i): Path to the multi-slice data that needs to be integrated. All slices should be contained in a single RDS or RData (the variable name after loading should be data) file. Refer to the example data for the specific format. This is a required parameter.
  # output_file (abbreviation: o): Path to the directory where the output results will be saved. This is a required parameter.
  # sample (abbreviation: s): Name for the output result files. It is recommended to use the model name, for example, setting it as "Banksy" will result in output files named "Banksy.h5ad". This is a required parameter.
  # nclust (abbreviation: c): Number of identified domains. This is a required parameter.
  # dataType (abbreviation: d): The type of input data, which should be set to "RData" or "RDS". This is an optional parameter, with the default set to "RData".
  # hvgs (abbreviation: h): The number of highly variable genes identified if data normalization is performed. This is an optional parameter, with the default set to 2000.
  # runNormalization (abbreviation: n): A boolean variable indicating whether data normalization should be performed before integration. This is an optional parameter, with the default set to TRUE.
  # k_geom (abbreviation: k): Hyperparameter for the Banksy model. This is an optional parameter, with the default set to c(15, 30).
  # lambda (abbreviation: l): Hyperparameter for the Banksy model. This is an optional parameter, with the default set to 0.8.
  # SEED (abbreviation: S): Random seed for model execution. This is an optional parameter, with the default set to 1000.
  
  # Step 1: Load input data
  # Depending on the data type (RData or RDS), load the data accordingly
  if (dataType == "RData") {
    load(input_file)
  } else if (dataType == "RDS") {
    data <- readRDS(input_file)
  } else {
    print("You should input data as RData or RDS format")
  }
  
  # Step 2: Stagger spatial coordinates for data
  # Adjust the spatial coordinates by normalizing and adjusting based on group
  locs <- spatialCoords(data)
  locs <- cbind(locs, sample_id = factor(data$slices))
  locs_dt <- data.table(locs)
  colnames(locs_dt) <- c("sdimx", "sdimy", "group")
  locs_dt[, sdimx := sdimx - min(sdimx), by = group]
  global_max <- max(locs_dt$sdimx) * 1.5
  locs_dt[, sdimx := sdimx + group * global_max]
  locs <- as.matrix(locs_dt[, 1:2])
  rownames(locs) <- colnames(data)
  spatialCoords(data) <- locs
  
  # Step 3: Data normalization (optional)
  # If normalization is required, normalize the data using Seurat
  if (runNormalization) {
    seu <- as.Seurat(data, data = NULL)
    seu <- FindVariableFeatures(seu, nfeatures = hvgs)  # Identify HVGs
    scale_factor <- median(colSums(assay(data, "counts")))  # Calculate scaling factor
    seu <- NormalizeData(seu, scale.factor = scale_factor, normalization.method = "RC")  # Normalize data
    
    # Add normalized data to SpatialExperiment
    aname <- "normcounts"
    assay(data, aname, withDimnames = FALSE) <- GetAssayData(seu)
    data <- data[VariableFeatures(seu), ]  # Subset to HVGs
  } else {
    # If no normalization, just copy the data to the assay
    seu <- as.Seurat(data, data = NULL)
    aname <- "normcounts"
    assay(data, aname, withDimnames = FALSE) <- GetAssayData(seu)
  }
  
  # Step 4: Record the time and memory usage for Banksy computation
  start_time <- Sys.time()
  start_mem <- pryr::mem_used()
  
  # Step 5: Run the Banksy method
  compute_agf <- TRUE
  data <- computeBanksy(data, assay_name = aname, compute_agf = compute_agf, k_geom = k_geom)  # Apply Banksy
  npcs <- 20
  use_agf <- TRUE
  data <- runBanksyPCA(data, use_agf = use_agf, lambda = lambda, npcs = npcs, seed = SEED)  # Perform PCA
  
  # Step 6: Record the end time and memory usage
  end_time <- Sys.time()
  end_mem <- pryr::mem_used()
  
  # Log the time and memory usage
  time_taken <- end_time - start_time
  mem_used <- end_mem - start_mem
  log_df <- data.frame(
    Memory_Usage_MiB = as.numeric(mem_used) / (1024 ^ 2),
    Execution_Time_s = as.numeric(time_taken, units = "secs")
  )
  cat("Time taken for computeBanksy: ", time_taken, "\n")
  cat("Memory used for computeBanksy: ", mem_used, "\n")
  
  # Step 7: Run Harmony batch correction on PCA
  harmony_embedding <- HarmonyMatrix(
    data_mat = reducedDim(data, paste("PCA_M1_lam", lambda, sep = "")),
    meta_data = colData(data),
    vars_use = c("slices"),
    do_pca = FALSE,
    max.iter.harmony = 20,
    verbose = FALSE
  )
  reducedDim(data, "Harmony_BANKSY") <- harmony_embedding
  
  # Step 8: Run UMAP for dimensionality reduction and visualization
  data <- runBanksyUMAP(data, use_agf = TRUE, lambda = lambda, npcs = npcs)
  data <- runBanksyUMAP(data, dimred = "Harmony_BANKSY")
  
  # Step 9: Cluster the data based on the corrected Harmony embedding
  optimal_resolution = find_optimal_resolution(data, nclust)  # Find the optimal resolution for clustering
  data <- clusterBanksy(data, dimred = "Harmony_BANKSY", resolution = optimal_resolution, seed = SEED)
  
  # Step 10: Add spatial coordinates to the data object
  data@colData@listData[["X"]] <- data@int_colData@listData[["spatialCoords"]][,1]
  data@colData@listData[["Y"]] <- data@int_colData@listData[["spatialCoords"]][,2]
  
  # Step 11: Convert data to Seurat format and save as HDF5 file
  data_seurat = as.Seurat(data, data = NULL)
  i <- paste("clust_Harmony_BANKSY_k50_res", optimal_resolution, sep = "")
  data_seurat$predicted_domain = data_seurat@meta.data[, i]
  
  # Create embedding object and add it to Seurat object
  embedding_exp <- data@assays@data@listData[["H1"]]
  embedding_exp <- t(embedding_exp)
  colnames(embedding_exp) <- paste0("embedding_", 1:ncol(embedding_exp))
  embedding <- CreateDimReducObject(embeddings = embedding_exp, key = "embedding_", assay = "originalexp")
  data_seurat[["embedding"]] <- embedding
  
  # Save the Seurat object as an H5Seurat file and convert to h5ad format
  h5Seurat_file_1 <- paste(output_file, paste(sample, "h5Seurat", sep = "."), sep = "/")
  SaveH5Seurat(data_seurat, filename = h5Seurat_file_1)
  Convert(h5Seurat_file_1, dest = "h5ad")
  
  # Save execution stats to a CSV file
  write.table(log_df, paste(output_file, paste(sample, "_stats.csv", sep = ""), sep = "/"), col.names = T, row.names = F, sep = ",", quote = F)
}

# Function to find the optimal clustering resolution based on the target number of clusters
# This function adjusts the resolution in a binary search manner until the number of clusters is close to the target
find_optimal_resolution <- function(data, target_clusters, tol = 1e-3, max_iter = 100) {
  lower_res <- 0.0
  upper_res <- 5.0
  optimal_res <- (lower_res + upper_res) / 2
  iter <- 0
  
  # Loop to adjust the resolution and check if the number of clusters matches the target
  while (iter < max_iter) {
    iter <- iter + 1
    data <- clusterBanksy(data, dimred = "Harmony_BANKSY", resolution = optimal_res, seed = 1000)
    data_seurat = as.Seurat(data, data = NULL)
    i <- paste("clust_Harmony_BANKSY_k50_res", optimal_res, sep = "")
    num_clusters = length(unique(data_seurat@meta.data[, i]))
    
    # Check if the number of clusters is within the tolerance of the target
    if (abs(num_clusters - target_clusters) <= tol) {
      return(optimal_res)
    } else if (num_clusters < target_clusters) {
      lower_res <- optimal_res
    } else {
      upper_res <- optimal_res
    }
    
    optimal_res <- (lower_res + upper_res) / 2
  }
  
  return(optimal_res)
}

# Command-line argument parsing
# Define the options that can be passed to the script
opt = matrix(c("input_file", "i", 1, "character", "Input file path",
               "output_file", "o", 1, "character", "Output file path",
               "sample", "s", 1, "character", "Sample name",
               "nclust", "c", 1, "numeric", "Cluster number",
               "dataType", "d", 2, "character", "the type of data storage(RData or RDS)",
               "hvgs", "h", 2, "numeric", "Number of HVGs (default: 2000)",
               "runNormalization", "n", 2, "logical", "Run data normalization (default: TRUE)",
               "k_geom", "k", 2, "numeric", "hiperparameter k",
               "lambda", "l", 2, "numeric", "hiperparameter lambda",
               "SEED", "S", 2, "numeric", "Random seed. default is 1000"),
             byrow = TRUE, ncol = 5)

args = getopt(opt)

if (is.null(args$dataType)) {
  args$dataType = "RData"
}

if (is.null(args$hvgs)) {
  args$hvgs = 2000
}

if (is.null(args$runNormalization)) {
  args$runNormalization = TRUE
}

if (is.null(args$k_geom)) {
  args$k_geom = c(15,30)
}

if (is.null(args$lambda)) {
  args$lambda = 0.8
}

if (is.null(args$SEED)) {
  args$SEED = 1000
}

input_file <- args$input_file
output_file <- args$output_file
sample <- args$sample
nclust <- as.numeric(args$nclust)
dataType <- args$dataType
hvgs <- args$hvgs
runNormalization <- args$runNormalization
k_geom <- args$k_geom
lambda <- args$lambda
SEED <- args$SEED

run_banksy(input_file, output_file, sample, nclust, dataType, hvgs, runNormalization, k_geom, lambda, SEED)
